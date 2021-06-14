import numpy as np
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3


#cave: not all mhc filtered yet

# not needed for now
# class proteaseForMutID:
#     def __init__(self, chymo,
#                  gluc, lysc, tryp):
#         self.chymo = chymo
#         self.gluc = gluc
#         self.lysc = lysc
#         self.tryp = tryp


def peptidesparser(filename):
    with open(filename) as pep:
        header = pep.readline().split("\t")
        # set up indexes
        mutid_idx = header.index("Mutation names")
        mutated_idx = header.index("Mutated")
        rev_idx = header.index("Reverse")
        cont_idx = header.index("Potential contaminant")
        seq_idx = header.index("Sequence")
        # idx for the different proteases
        chymo_idx = header.index("Experiment Chymotrypsin")
        gluc_idx = header.index("Experiment GluC")
        lysc_idx = header.index("Experiment LysC")
        tryp_idx = header.index("Experiment Trypsin")
        # initialise counter for mutated peptides found per protease
        cnt_chymo = 0
        cnt_gluc = 0
        cnt_lysc = 0
        cnt_tryp = 0
        # initialise counter for peptides found per protease
        cnt_chymoall = 0
        cnt_glucall = 0
        cnt_lyscall = 0
        cnt_trypall = 0
        # initialise lists of unique mutations
        lchymo2, lgluc2, llysc2, ltryp2 = [], [], [], []

        for line in pep:
            l = line.split("\t")
            if "+" != (l[rev_idx] or l[cont_idx]):  # filter for not reverse and not contamination
                if "Yes" == l[mutated_idx] and "" != l[mutid_idx]:  # filter for mutated
                    # adds 1 to the count if there is evidence for the peptide by the respective protease
                    if l[chymo_idx] != "":
                        cnt_chymo += 1
                    if l[gluc_idx] != "":
                        cnt_gluc += 1
                    if l[lysc_idx] != "":
                        cnt_lysc += 1
                    if l[tryp_idx] != "":
                        cnt_tryp += 1
                    # now we look at individual mutations and fill the list of each protease with mutids they account for
                    for mutid in l[mutid_idx].split(";"):
                        if l[chymo_idx] != "":
                            lchymo2.append(mutid)
                        if l[gluc_idx] != "":
                            lgluc2.append(mutid)
                        if l[lysc_idx] != "":
                            llysc2.append(mutid)
                        if l[tryp_idx] != "":
                            ltryp2.append(mutid)
                else:  # do not filter for mutated, i.e. all peptides are counted and which protease gives evidence for them
                    if l[chymo_idx] != "":
                        cnt_chymoall += 1
                    if l[gluc_idx] != "":
                        cnt_glucall += 1
                    if l[lysc_idx] != "":
                        cnt_lyscall += 1
                    if l[tryp_idx] != "":
                        cnt_trypall += 1
        # generate output below
        mutperprotease = {"Chymotrypsin": cnt_chymo, "GluC": cnt_gluc, "LysC": cnt_lysc,
                          "Trypsin": cnt_tryp}  # mutated peptides found per protease
        mutperproteaseall = {"Chymotrypsin": cnt_chymoall, "GluC": cnt_glucall, "LysC": cnt_lyscall,
                             "Trypsin": cnt_trypall}  # peptides found per protease
        uniquemutdict = {"Chymotrypsin": lchymo2, "GluC": lgluc2, "LysC": llysc2, "Trypsin": ltryp2} # unique mutations found per protease
        alluniquemut = uniquemutdict["Chymotrypsin"] + uniquemutdict["GluC"] + uniquemutdict["LysC"] + uniquemutdict[
            "Trypsin"] # list of all unique mutations found per protease

    return mutperprotease, mutperproteaseall, uniquemutdict, alluniquemut


mutperprotease, mutperproteaseall, unmutdict, alluniquemut = peptidesparser(
    "C:\\Users\\Max Gerwien\\VariantProtocol\\fastrun_allfiles\\peptides_gutchymo.txt")
#check out the output below

# print("All Peptides\n", mutperproteaseall)
# print("Mutated Peptides\n", mutperprotease)
# # total number of unique mutations
print("All Unique Mutations", len(set([item for list in unmutdict.values() for item in list])))
# # number of unique mutations per protease
# for k, v in unmutdict.items():
#     print(k, len(set(v)))

# barplot of number of unique mutations found for combinations of proteases
# this below is very unelegant

chymo = set(unmutdict["Chymotrypsin"])
gluc = set(unmutdict["GluC"])
lysc = set(unmutdict["LysC"])
tryp = set(unmutdict["Trypsin"])
chymo_gluc = set(unmutdict["Chymotrypsin"] + unmutdict["GluC"])
chymo_lysc = set(unmutdict["Chymotrypsin"] + unmutdict["LysC"])
chymo_tryp = set(unmutdict["Chymotrypsin"] + unmutdict["Trypsin"])
gluc_lysc = set(unmutdict["GluC"] + unmutdict["LysC"])
gluc_tryp = set(unmutdict["GluC"] + unmutdict["Trypsin"])
lysc_tryp = set(unmutdict["LysC"] + unmutdict["Trypsin"])
chymo_gluc_lysc = set(unmutdict["Chymotrypsin"] + unmutdict["GluC"] + unmutdict["LysC"])
chymo_gluc_tryp = set(unmutdict["Chymotrypsin"] + unmutdict["GluC"] + unmutdict["Trypsin"])
chymo_lysc_tryp = set(unmutdict["Chymotrypsin"] + unmutdict["LysC"] + unmutdict["Trypsin"])
gluc_lysc_tryp = set(unmutdict["GluC"] + unmutdict["LysC"] + unmutdict["Trypsin"])
chymo_gluc_lysc_tryp = set(unmutdict["Chymotrypsin"] + unmutdict["GluC"] + unmutdict["LysC"] + unmutdict["Trypsin"])
x_idx = ["C", "G", "L", "T", "CG", "CL", "CT", "GL", "GT", "LT", "CGL", "CGT", "CLT", "GLT", "ALL"]
y_bar = [len(chymo), len(gluc), len(lysc), len(tryp), len(chymo_gluc), len(chymo_lysc), len(chymo_tryp), len(gluc_lysc),
         len(gluc_tryp), len(lysc_tryp), len(chymo_gluc_lysc), len(chymo_gluc_tryp), len(chymo_lysc_tryp),
         len(gluc_lysc_tryp), len(chymo_gluc_lysc_tryp)]

plt.bar(x_idx, y_bar, width=0.3)
plt.xlabel("Protease")
plt.ylabel("# unique mutations identified")
plt.title("Unique mutations identified for protease combinations")
plt.show()



def variantsparser(filename):
    with open(filename) as var:
        header = var.readline().split("\t")
        mutid_idx = header.index("Id")
        transcr_idx = header.index("Transcripts")
        allmut = []  # all unique mutations (mutids)
        allnonsynmut = []  # all unique non-synonymous mutations (mutids)
        allmuttranscripts = [] # all unique non-synonymous mutations gene ids of transcripts
        for line in var:
            l = line.split("\t")
            allmut.append(l[mutid_idx])
            if "NonSynonymous" in l[transcr_idx]:
                allnonsynmut.append(l[mutid_idx])
                for ids in l[transcr_idx].split(";"):
                    for el in ids.split(":"):
                        if el.startswith("ENSG"):
                            allmuttranscripts.append(el)
        allmuttranscripts = set(allmuttranscripts)
    return allmut, allnonsynmut, allmuttranscripts


#
allmut, allnonsynmut, allmuttranscripts = variantsparser(
    "C:\\Users\\Max Gerwien\\VariantProtocol\\fastrun_allfiles\\variants.txt")
#check out the output below

# print(len(allmut))
# print(len(allnonsynmut))
# print("allmuttranscripts", len(allmuttranscripts))

x = ["All Mutations", "Non-Synonymous Mutations", "Mutations Peptide Level"]
y_bar = [len(set(allmut)), len(set(allnonsynmut)), len(set([item for list in unmutdict.values() for item in list]))]
x_idx = np.arange(len(x))
barwidth = 0.3
plt.bar(x_idx, y_bar, width=barwidth)
plt.xticks(ticks=x_idx, labels=x)
plt.show()
# plt.figure()
y_venn = [set(allmut), set(allnonsynmut), set([item for list in unmutdict.values() for item in list])]
venn3(y_venn,
      set_labels=(x))
plt.title("Mutations found")
plt.show()


def readFastaFile(proteinsfastafile,
                  mhcinref):  # filename is the proteins.fa and we take mhc annotations from the reference genome
    with open(proteinsfastafile) as fs, open(mhcinref) as mhc:
        transitionsdict = {}
        mhc_list = [line.split(".")[0][1:] for line in mhc if "histocompatibility" in line]  # read mhc genes into array
        for line in fs:
            line = line.rstrip()
            if line.startswith('>'):  # fasta header found
                splline = line[1:].split('\t')  # a list with ens.id and the mutation column as splline
                if splline[0] not in mhc_list:  # only if not a mhc gene
                    if len(splline) == 2:  # this is only done if we have a mutation, i.e. there are two columns
                        for mutation in splline[1].split(';'):  # looks at all the mutations
                            splmut = mutation.split(":")
                            if splmut[0] in alluniquemut: # i.e. mutid found on peptide level
                                tmptransition = splmut[1][0] + splmut[1][-1] # concatenate AA transition
                                if tmptransition in transitionsdict.keys(): #fill dict with transitions and count them
                                    transitionsdict[tmptransition] += 1
                                else:
                                    transitionsdict[tmptransition] = 1

    return transitionsdict


transitionsdict = readFastaFile(
    "C:\\Users\\Max Gerwien\\VariantProtocol\\fastrun_allfiles\\proteins - Copy.fa",
    "C:\\Users\\Max Gerwien\\VariantProtocol\\fastrun_allfiles\\Homo_sapiens.GRCh38.pep.all - Copy.fa")

# barplot of transition frequencies; better: heatmap!

x = list(transitionsdict.keys())
y_bar = [float(i) for i in transitionsdict.values()]
x_idx = np.arange(len(x))
barwidth = 0.3
plt.bar(x_idx, y_bar, width=barwidth)
plt.xticks(ticks=x_idx, labels=x)
plt.show()

# heatmap!
# set up dict of aa:idx
idxlist = list(range(0, 20))
twentyaa = sorted([i for i in "LITSFANMPGKQYVHWDERC"])
aa2idx = {twentyaa[i]: list(range(0, 20))[i] for i in idxlist}
# add 0s (= no transitions) to the transitionsdict for AA combinations for which no transitions were found;
# also add 0s for same aa, eg AA or EE (since this is sth we cannot measure, make it maybe gray or sth in the final fig)
for rowidx in aa2idx.keys():
    for colidx in aa2idx.keys():
        if rowidx + colidx not in transitionsdict.keys():
            transitionsdict[rowidx + colidx] = 0
# construct 20x20 data matrix from transitions dict
mat = [[transitionsdict[rowidx + colidx] for rowidx in aa2idx.keys()] for colidx in aa2idx.keys()]

# # below heatmap attempt

x = range(0, 20)
y = range(0, 20)
x, y = np.meshgrid(x, y)

intensity = np.array(mat)
plt.pcolor(x, y, mat, shading="auto")
plt.colorbar()
x_idx = list(range(0, 20))
y_idx = list(range(0, 20))

plt.xticks(ticks=x_idx, labels=twentyaa)
plt.yticks(ticks=y_idx, labels=twentyaa)

plt.xlabel("aa before")
plt.ylabel("aa after")
plt.title("ab")
plt.show()
print(mat)


def expression_rpkm_parser(filename):
    with open(filename) as expr:
        header = expr.readline().split("\t")
        pei200r1_idx = header.index("Caltech_HeLa-S3_cell_PE_i200_r1_log2(x*10^9/y)_Max transcript(cds) length")
        pei200r2_idx = header.index("Caltech_HeLa-S3_cell_PE_i200_r2_log2(x*10^9/y)_Max transcript(cds) length")
        sei0r1_idx = header.index("Caltech_HeLa-S3_cell_SE_i0_r1_log2(x*10^9/y)_Max transcript(cds) length")
        sei0r2_idx = header.index("Caltech_HeLa-S3_cell_SE_i0_r2_log2(x*10^9/y)_Max transcript(cds) length")
        transcr_idx = header.index("Transcript(Protein) ids")
        next(expr)
        transcr_lvl_all = [] # list of median transcription levels of all genes
        transcr_lvl_var = [] # list of median transcription levels of genes carrying mutation
        for line in expr:
            l = line.rstrip().split("\t")
            transcr_lvl_all.append(
                np.median(list(map(float, [l[pei200r1_idx], l[pei200r2_idx], l[sei0r1_idx], l[sei0r2_idx]]))))
            for vartranscr in allmuttranscripts: # list of mutated ensembl gene names from variants.txt
                if vartranscr in line:
                    transcr_lvl_var.append(
                        np.median(list(map(float, [l[pei200r1_idx], l[pei200r2_idx], l[sei0r1_idx], l[sei0r2_idx]]))))

    return transcr_lvl_all, transcr_lvl_var


transcr_lvl_all, transcr_lvl_var = expression_rpkm_parser(
    "C:\\Users\\Max Gerwien\\VariantProtocol\\fastrun_allfiles\\expression.rpkm.txt")

# check out output below

# print("no. transcr_lvl_all", len(transcr_lvl_all))
# print("no. transcr_lvl_var", len(transcr_lvl_var))
# print("max", max(transcr_lvl_all))

# generate histogram "Transcription Distribution of All Genes vs Mutated Genes"

xrange = [min(transcr_lvl_all), max(transcr_lvl_all)]
plt.xticks(ticks=xrange, labels=xrange)

plt.hist(transcr_lvl_all, bins="auto", range=xrange, label="All")
plt.hist(transcr_lvl_var, bins="auto", range=xrange, label="With Mutations")
plt.legend(loc='upper right')
plt.xlabel("RPKM")
plt.ylabel("# of Genes")
plt.title("Transcription Distribution of All Genes vs Mutated Genes")

plt.show()

