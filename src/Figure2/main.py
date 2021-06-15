import gzip
import json
import math

import matplotlib
import matplotlib.pyplot as plt
import upsetplot


class Mutation:
    def __init__(self, id, position, aaBefore, aaAfter):
        self.id = id
        self.position = position
        self.aaBefore = aaBefore
        self.aaAfter = aaAfter

    @staticmethod
    def parseString(line):
        spl = line.split(':')
        return Mutation(spl[0], int(spl[1][1:-1]), spl[1][0], spl[1][-1])


class Protein:
    def __init__(self, id, sequence, mutations):
        self.id = id
        self.sequence = sequence
        self.mutations = mutations

    @staticmethod
    def readFastaFile(filename):
        records = {}
        with gzip.open(filename, "rt") as fs:
            recordId = ""
            for line in fs:
                if line[0] == '>':
                    spl = line.rstrip('\n').split('\t')
                    recordId = spl[0][1:]
                    mutations = {}
                    if len(spl[1]) != 0:
                        for s in spl[1].split(';'):
                            mutation = Mutation.parseString(s)
                            mutations[mutation.id] = mutation
                    records[recordId] = Protein(
                        recordId,
                        "",
                        mutations
                    )
                else:
                    records[recordId].sequence += line.rstrip('\n')
        return records


class Peptide:
    def __init__(self, sequence, mutationIds, proteins, experiment2intensity):
        self.sequence = sequence
        self.mutationIds = mutationIds
        self.proteins = proteins
        self.experiment2intensity = experiment2intensity

    @staticmethod
    def readMaxQuantPeptideFile(filename, columns, experimentNames):
        peptides = []
        with gzip.open(filename, "rt") as fs:
            header = fs.readline().rstrip().split("\t")
            sequenceIndex = header.index(columns["sequence"])
            proteinsIndex = header.index(columns["proteins"])
            mutatedIndex = header.index(columns["mutated"])
            mutationNamesIndex = header.index(columns["mutatedNames"])
            filterIndexes = [header.index(columns[i]) for i in ["reverse", "contaminant"]]
            experimentIndexes = {experimentName: header.index(f"Intensity {experimentName}") for experimentName in
                                 experimentNames}
            for line in fs:
                spl = line.rstrip().split('\t')
                if sum([spl[i] == '+' for i in filterIndexes]) > 0 or spl[mutatedIndex] == "Mixed":
                    continue
                proteins = [protein for protein in spl[proteinsIndex].split(';') if protein.startswith("ENSP")]
                experiments = {experimentName: int(spl[index]) for experimentName, index in experimentIndexes.items()}
                mutationIds = []
                if spl[mutatedIndex] == "Yes":
                    mutationIds = spl[mutationNamesIndex].split(';')
                if len(proteins) == 0 or sum([value for experimentName, value in experiments.items()]) == 0:
                    continue
                peptide = Peptide(
                    spl[sequenceIndex],
                    mutationIds,
                    proteins,
                    experiments
                )
                peptides.append(peptide)
        return peptides

    @staticmethod
    def getUniqueMutatedPeptidesCount(peptides, experimentNames):
        cnt = {experimentName: 0 for experimentName in experimentNames}
        for peptide in peptides:
            if len(peptide.mutationIds) == 0:
                continue
            for experimentName in experimentNames:
                if peptide.experiment2intensity[experimentName] == 0:
                    continue
                cnt[experimentName] += 1
        return cnt

    @staticmethod
    def getUniqueExpressedMutationsCount(peptides, experimentNames):
        cnt = {experimentName: set() for experimentName in experimentNames}
        for peptide in peptides:
            if len(peptide.mutationIds) == 0:
                continue
            for experimentName in experimentNames:
                if peptide.experiment2intensity[experimentName] == 0:
                    continue
                for mutationId in peptide.mutationIds:
                    cnt[experimentName].add(mutationId)
        return {experimentName: len(cnt[experimentName]) for experimentName in experimentNames}


def plotUpSet(peptides, experimentNames, filename):
    exp2num = {experimentNames[i]: 2 ** i for i in range(len(experimentNames))}
    cnt = [0 for i in range(2 ** len(experimentNames))]
    names = []
    for i in range(2 ** len(experimentNames)):
        s = []
        for exp in experimentNames:
            if exp2num[exp] & i > 0:
                s.append(exp)
        names.append(s)
    mutId2ind = {}
    for peptide in peptides:
        if len(peptide.mutationIds) == 0:
            continue
        ind = 0
        for exp in experimentNames:
            if peptide.experiment2intensity[exp] == 0:
                continue
            ind += exp2num[exp]
        for mutationId in peptide.mutationIds:
            if mutationId in mutId2ind:
                mutId2ind[mutationId] |= ind
            else:
                mutId2ind[mutationId] = ind
    for mutId, ind in mutId2ind.items():
        cnt[ind] += 1
    d = upsetplot.from_memberships(
        names[1:],
        data=cnt[1:]
    )
    upsetplot.plot(d, show_counts='%d')
    plt.savefig(filename)


def plotA(params, peptides, filename):
    plotUpSet(
        peptides,
        params["Bekker-Jensen.2017"]["input"]["proteomics"]["peptides"]["experiments"],
        filename
    )
    print("len(peptides)", len(peptides))
    print("getUniqueMutatedPeptidesCount:", Peptide.getUniqueMutatedPeptidesCount(
        peptides,
        params["Bekker-Jensen.2017"]["input"]["proteomics"]["peptides"]["experiments"]
    ))
    print("getUniqueExpressedMutationsCount:", Peptide.getUniqueExpressedMutationsCount(
        peptides,
        params["Bekker-Jensen.2017"]["input"]["proteomics"]["peptides"]["experiments"]
    ))


def plotB(params, peptides, proteins, filename, cmapFilename, circleFilename):
    aas = "ACDEFGHIKLMNPQRSTVWY"
    id2mutations = {}
    for peptide in peptides:
        if len(peptide.mutationIds) == 0:
            continue
        for mutationId in peptide.mutationIds:
            for proteinId in peptide.proteins:
                if mutationId not in proteins[proteinId].mutations:
                    continue  # TODO what is going on here?
                mut = proteins[proteinId].mutations[mutationId]
                if mutationId not in id2mutations:
                    id2mutations[mutationId] = set()
                id2mutations[mutationId].add(f"{mut.aaBefore}{mut.aaAfter}")
    id2mutations_all = {}
    for protein, fastaRecord in proteins.items():
        for mutationId, mutationRecord in fastaRecord.mutations.items():
            if mutationId not in id2mutations_all:
                id2mutations_all[mutationId] = set()
            id2mutations_all[mutationId].add(f"{mutationRecord.aaBefore}{mutationRecord.aaAfter}")
    before2after = {}
    for i in range(len(aas)):
        for j in range(len(aas)):
            if i == j:
                continue
            before2after[f"{aas[i]}{aas[j]}"] = [0, 0]
    for mutationId, mutationSet in id2mutations.items():
        for mutation in mutationSet:
            before2after[mutation][0] += 1
    for mutationId, mutationSet in id2mutations_all.items():
        for mutation in mutationSet:
            before2after[mutation][1] += 1
    before2after['IL'] = [0, 0]
    before2after['LI'] = [0, 0]

    fig, ax = plt.subplots(1, figsize=(5, 5))
    ax.set(xlim=(-1, len(aas)), ylim=(-1, len(aas)))
    maxv = max([v[1] for k, v in before2after.items()])
    for i in range(len(aas)):
        ax.hlines(y=i, xmin=0, xmax=len(aas) - 1, colors="grey", lw=0.8)
        ax.vlines(x=i, ymin=0, ymax=len(aas) - 1, colors="grey", lw=0.8)
    norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
    cmap = matplotlib.cm.get_cmap('jet')
    m = matplotlib.cm.ScalarMappable(
        norm=norm,
        cmap=cmap
    )
    for i in range(len(aas)):
        for j in range(len(aas)):
            if i == j:
                continue
            aa = f"{aas[i]}{aas[j]}"
            r = (math.log2(before2after[aa][1] + 1) / math.log2(maxv + 1)) * 0.5
            if r == 0:
                continue
            c = plt.Circle(
                (i, j),
                r,
                lw=0,
                color=m.to_rgba(before2after[aa][0] / before2after[aa][1])
            )
            ax.add_artist(c)
    plt.xticks([i for i in range(len(aas))], [c for c in aas])
    ax.set_xlabel("Reference Amino Acid")
    plt.yticks([i for i in range(len(aas))], [c for c in aas])
    ax.set_ylabel("Alternative Amino Acid")
    plt.savefig(filename)

    fig, ax = plt.subplots(1)
    matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm,
                                     orientation='horizontal')
    plt.savefig(cmapFilename)

    circles = [1]
    ii = 1
    while ii * 25 < maxv:
        circles.append(ii * 25)
        ii += 1
    fig, ax = plt.subplots(1, figsize=(5, 5))
    ax.set(xlim=(-1, len(circles)), ylim=(-1, 1), aspect="equal")
    for i in range(len(circles)):
        r = (math.log2(circles[i] + 1) / math.log2(maxv + 1)) * 0.5
        c = plt.Circle(
            (i, 0),
            r,
            lw=0,
            color="grey"
        )
        ax.add_artist(c)
    plt.savefig(circleFilename)




def plot(params):
    peptides = Peptide.readMaxQuantPeptideFile(
        params["Bekker-Jensen.2017"]["input"]["proteomics"]["peptides"]["file"],
        params["Bekker-Jensen.2017"]["input"]["proteomics"]["peptides"]["columns"],
        params["Bekker-Jensen.2017"]["input"]["proteomics"]["peptides"]["experiments"]
    )
    proteins = Protein.readFastaFile(
        params["Bekker-Jensen.2017"]["input"]["transcriptomics"]["proteins"]
    )
    plotA(
        params,
        peptides,
        params["Bekker-Jensen.2017"]["output"]["figureA"]
    )
    plotB(
        params,
        peptides,
        proteins,
        params["Bekker-Jensen.2017"]["output"]["figureB"],
        params["Bekker-Jensen.2017"]["output"]["figureB_cmap"],
        params["Bekker-Jensen.2017"]["output"]["figureB_circle"]
    )


if __name__ == "__main__":
    with open('parameters.json', 'r') as fs:
        plot(json.loads(fs.read()))
