import gzip
import json

# import numpy as np
import pandas as pd
import logomaker as lm
import matplotlib.pyplot as plt


#
# print(lm.list_font_names())
# #plt.ion()
# # provided (for na)
# crp_df = -lm.get_example_matrix('crp_energy_matrix')
# crp_df.head()
# # print(crp_df)
# logo = lm.Logo(crp_df,font_name="Arial")
# plt.show()
# print(25+23)


# # for peptides
# twentyaa = sorted([i for i in "LITSFANMPGKQYVHWDERC"])
# random_ppm = np.random.dirichlet(np.ones(20), size=9)
# df=pd.DataFrame(random_ppm, columns= twentyaa)
# print(df)
# ​
# logo = lm.Logo(df,font_name="Arial")
# plt.show()

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
            # header.index(f"Intensity {experimentName}")
            experimentIndexes = {experimentName: [i for i in range(len(header))
                                                  if "Intensity" in header[i] and experimentName in header[i]]
                                 for experimentName in experimentNames}
            for line in fs:
                spl = line.rstrip().split('\t')
                if sum([spl[i] == '+' for i in filterIndexes]) > 0 or spl[mutatedIndex] == "Mixed":
                    continue
                proteins = [protein for protein in spl[proteinsIndex].split(';') if protein.startswith("ENSP")]
                experiments = {experimentName: sum([int(spl[index]) for index in indexes])
                               for experimentName, indexes in experimentIndexes.items()}
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


def statistics(peptides, experimentNames):
    cnt = {experimentName: [0, 0] for experimentName in experimentNames}
    for peptide in peptides:
        if len(peptide.sequence) != 9:
            continue
        for experimentName in experimentNames:
            if peptide.experiment2intensity[experimentName] == 0:
                continue
            cnt[experimentName][1] += 1
            if len(peptide.mutationIds) == 0:
                continue
            cnt[experimentName][0] += 1
    return cnt


def plotImpl(peptides, experimentNames, filename, n=9):
    sequences = {experimentName: [[], []] for experimentName in experimentNames}
    for peptide in peptides:
        if len(peptide.sequence) != n:
            continue
        for experimentName in experimentNames:
            if peptide.experiment2intensity[experimentName] == 0:
                continue
            if len(peptide.mutationIds) == 0:
                sequences[experimentName][0].append(peptide.sequence)
            else:
                sequences[experimentName][1].append(peptide.sequence)
    aas = "ACDEFGHIKLMNPQRSTVWY"
    aa2index = {aas[i]: i for i in range(len(aas))}
    m = len(aas)
    fig, ax = plt.subplots(len(experimentNames), 2, figsize=(5, 5))

    for experimentIndex in range(len(experimentNames)):
        for mutatedIndex in range(2):
            data = [[0 for j in range(m)] for i in range(n)]
            print(experimentNames[experimentIndex], mutatedIndex, len(sequences[experimentNames[experimentIndex]][mutatedIndex]))
            for sequence in sequences[experimentNames[experimentIndex]][mutatedIndex]:
                for i in range(len(sequence)):
                    data[i][aa2index[sequence[i]]] += 1
            df = pd.DataFrame(data, columns=list(aas))
            lm.Logo(
                df,
                color_scheme="chemistry",
                font_name="Arial",
                fade_probabilities=True,
                stack_order='small_on_top',
                ax=ax[experimentIndex, mutatedIndex]
            )
    plt.savefig(filename)
    # plt.show()


def plot(params):
    ps = params["Sarkizova.2019"]["input"]["proteomics"]["peptides"]
    peptides = Peptide.readMaxQuantPeptideFile(
        ps["file"],
        ps["columns"],
        ps["experiments"]
    )
    print(statistics(peptides, ps["experiments"]))
    plotImpl(
        peptides,
        ps["selectedExperiments"],
        params["Sarkizova.2019"]["output"]["figureAB"]
    )


if __name__ == "__main__":
    with open('parameters.json', 'r') as fs:
        plot(json.loads(fs.read()))
