import numpy as np
from copy import deepcopy
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


class varOperator:
    def __init__(self, cpP, tarSeq, yearStep):
        self.cpP = cpP
        self.tarSeq = tarSeq
        self.yearStep = yearStep
        self.minSvLen = 20

    def adjustCpNumByYear(self, te, nowYear, cpNum):
        return int(cpNum/((nowYear-te.bornTime)/self.yearStep))

    def geneSvedSeq(self, te, nowYear, cpNum, dirSvTypList, dirSvPosList):
        retList = []
        retList.append(te)
        cpNum = self.adjustCpNumByYear(te, nowYear, cpNum)
        for i in range(cpNum):
            nte = deepcopy(te)
            nte.id = f'{nte.id}_{te.cpedNum}'
            nte.seq = self.dirSv(nte.seq, dirSvTypList[i], dirSvPosList[i])
            nte.bornTime = nowYear
            retList.append(deepcopy(nte))
            te.addCpedNum(1)
        return retList

    def dirSv(self, oriSeq, dirSvTyp, dirSvPos):
        aln = pairwise2.align.globalxs(oriSeq, self.tarSeq, -1, -0.2)[0]
        alnedOriSeq, alnedTarSeq = aln.seqA, aln.seqB

        return oriSeq

