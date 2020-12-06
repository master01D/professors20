from Bio import Align
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
aligner = Align.PairwiseAligner()
aligner.mode = 'local'
epitopes = 'C:/E-proteins epitope.txt'
topro = open(epitopes, "r")
cont = topro.read()
part = ''
slices = []
sorter1 = []
sorter2 = []
score1 = []
score2 = []
restore = []
tick = 0
for i in range(len(cont)):
    if tick == 1:
        if cont[i].isupper():
            part = ''
            part = part + cont[i]
            tick = tick - 1
    elif cont[i].isupper():
        part = part + cont[i]
    elif tick == 0:
        slices.append(part)
        tick = tick + 1
seq0 = slices[0]
slices.pop(0)
def recounting(slices, seq0, sorter1, sorter2, score1, score2, restore):
    for i in range(len(slices)):
        if slices[i] not in seq0:
            alignments1 = aligner.align(seq0[0:len(slices[i])-1], slices[i])
            for alignment in alignments1:
                if "-" not in str(alignment) and str(alignment).count('|') > 1 and '.' not in str(alignment) and str(alignment)[0] == ' ':
                    check = str(alignment)
                    sorter1.append(check)
                    if sorter1 != []:
                        for variant in sorter1:
                            if variant.count('|') > 1:
                                score1.append(variant.count('|'))
            alignments2 = aligner.align(seq0[-len(slices[i]):], slices[i])
            for alignment in alignments2:
                if "-" not in str(alignment) and str(alignment).count('|') > 1 and '.' not in str(alignment) and str(alignment)[0] != ' ':
                    check = str(alignment)
                    sorter2.append(check)
                    if sorter2 != []:
                        for variant in sorter2:
                            if variant.count('|') > 1:
                                score2.append(variant.count('|'))
            if score1 == []:
                if score2 != []:
                    seq0 = seq0[:-max(score2)]
                    seq0 = seq0 + slices[i]
                    sorter1 = []
                    sorter2 = []
                    score1 = []
                    score2 = []
                else:
                    restore.append(slices[i])
            elif score2 == []:
                if score1!= []:
                    slices[i] = slices[i][:-max(score1)]
                    seq0 = slices[i] + seq0
                    sorter1 = []
                    sorter2 = []
                    score1 = []
                    score2 = []
                else:
                    restore.append(slices[i])
            elif max(score1) > max(score2):
                slices[i] = slices[i][:-max(score1)]
                seq0 = slices[i] + seq0
                sorter1 = []
                sorter2 = []
                score1 = []
                score2 = []
            elif max(score2) > max(score1):
                seq0 = seq0[:-max(score2)]
                seq0 = seq0 + slices[i]
                sorter = []
                score1 = []
                score2 = []
            else:
                restore.append(slices[i])
            print(seq0)
    return seq0, restore
seq0, slices = recounting(slices, seq0, sorter1, sorter2, score1, score2, restore)
print(seq0)
print('Done!')
