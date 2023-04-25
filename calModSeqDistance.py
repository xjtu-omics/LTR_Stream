from multiprocessing import Pool
import pickle
import numpy as np
import pandas as pd
import sys
danteD = sys.argv[1]
threads = int(sys.argv[2])
modSeqListDataFile = sys.argv[3]
# f'{danteD}/toCalRest.modSeqList.pkl'
modSeqList = ''
with open(modSeqListDataFile,'rb') as f:
    modSeqList = pickle.load(f)

def levenshtein(seq1, seq2):
    seq1 = seq1.split(',')
    seq2 = seq2.split(',')
    size_x = len(seq1) + 1
    size_y = len(seq2) + 1
    matrix = np.zeros((size_x, size_y))
    for x in range(size_x):
        matrix[x, 0] = x
    for y in range(size_y):
        matrix[0, y] = y

    for x in range(1,size_x):
        for y in range(1,size_y):
            if seq1[x - 1] == seq2[y - 1]:
                matrix[x, y] = min(
                    matrix[x - 1, y] + 1,
                    matrix[x - 1, y - 1],
                    matrix[x, y - 1] + 1
                )
            else:
                matrix[x, y] = min(
                    matrix[x - 1, y] + 1,
                    matrix[x - 1, y - 1] + 1,
                    matrix[x, y - 1] + 1
                )
    rel1 = (matrix[size_x - 1, size_y - 1])

    seq1.reverse()
    matrix = np.zeros((size_x, size_y))
    for x in range(size_x):
        matrix[x, 0] = x
    for y in range(size_y):
        matrix[0, y] = y

    for x in range(1,size_x):
        for y in range(1,size_y):
            if seq1[x - 1] == seq2[y - 1]:
                matrix[x, y] = min(
                    matrix[x - 1, y] + 1,
                    matrix[x - 1, y - 1],
                    matrix[x, y - 1] + 1
                )
            else:
                matrix[x, y] = min(
                    matrix[x - 1, y] + 1,
                    matrix[x - 1, y - 1] + 1,
                    matrix[x, y - 1] + 1
                )
    rel2 = (matrix[size_x - 1, size_y - 1])

    return min(rel1,rel2)

def calDisOfALine(lineNum):
    calRel = np.zeros(len(modSeqList))
    for j in range(len(modSeqList)):
        calRel[j] = levenshtein(modSeqList[lineNum],modSeqList[j])
    return calRel

if __name__=='__main__':
    calRel = []
    with Pool(int(threads)) as p:
        calRel = p.map(calDisOfALine,range(len(modSeqList)))
    disData = pd.DataFrame(calRel)
    disData.to_csv(f'{danteD}/toCalRest.distanceMat.csv',index=None,header=None)