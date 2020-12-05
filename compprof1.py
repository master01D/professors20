from Bio import Align
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
aligner = Align.PairwiseAligner()
aligner.mode = 'local'
epitopes = 'C:/Tester2.txt'
topro = open(epitopes, "r")
cont = topro.read()
part = ''
slices = []
sorter = []
score = []
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
def recounting(slices, seq0, sorter, score, restore):
    for i in range(len(slices)):
        print('Wait...')
        alignments = aligner.align(seq0, slices[i])
        for alignment in alignments:
            if "-" not in str(alignment) and str(alignment).count('|') > 1 and '.' not in str(alignment):
                check = str(alignment)
                sorter.append(check)
        if sorter != []:
            for variant in sorter:
                if variant.count('|') > 1:
                    score.append(variant.count('|'))
            if max(score) == len(seq0) or max(score) == len(slices[i]):
                if max(score) == len(seq0):
                    seq0 = slices[i]
                    sorter = []
                    score = []
                else:
                    sorter = []
                    score = []
            elif score != []:
                piece = sorter[score.index(max(score))]
                if piece[0] != ' ':
                    seq0 = seq0[:-max(score)]
                    seq0 = seq0 + slices[i]
                    sorter = []
                    score = []
                else:
                    slices[i] = slices[i][:-max(score)]
                    seq0 = slices[i] + seq0
                    sorter = []
                    score = []
        else:
             restore.append(slices[i])
    return seq0, restore
seq0, slices = recounting(slices, seq0, sorter, score,restore)
print(seq0)
print(slices)
print(recounting(slices, seq0, sorter, score, restore))
print('Done!')
    
    

