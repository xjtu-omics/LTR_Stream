import numpy as np
from copy import deepcopy
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


class varOperator:
    def __init__(self, cpP, tarSeq, yearStep, minSvLen):
        self.cpP = cpP
        self.tarSeq = tarSeq
        self.yearStep = yearStep
        self.minSvLen = minSvLen

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

    def getAlnScore(self, oriSeq):
        aln = pairwise2.align.globalxs(oriSeq, self.tarSeq, -1, -0.2)[0]
        return aln.score

    def dirSv(self, oriSeq, dirSvTyp, dirSvPos):
        aln = pairwise2.align.globalxs(oriSeq, self.tarSeq, -1, -0.2)[0]
        alnedOriSeq, alnedTarSeq = aln.seqA, aln.seqB

        ooPos = 0
        svOpFlag = False

        dirSvAlnPos = -1
        for i in range(len(alnedOriSeq)):
            if ooPos>=dirSvPos:
                dirSvOriPos = i
                break
            if alnedOriSeq[i] != '-':
                ooPos += 1

        oriWinFlag = False
        tarWinFlag = False
        oriWinSt = 0
        tarWinSt = 0
        opObj = None
        opAlnSt = None
        opAlnEd = None
        for i in range(dirSvAlnPos, len(alnedOriSeq)):
            if not oriWinFlag:
                if alnedOriSeq[i] == '-':
                    oriWinFlag = True
                    oriWinSt = i
            else:
                if alnedOriSeq[i] != '-':
                    if i-oriWinSt>=self.minSvLen:
                        opObj, opAlnSt, opAlnEd = 'ori', oriWinSt, i
                        break
                    else:
                        oriWinFlag = False
            if not tarWinFlag:
                if alnedTarSeq[i] == '-':
                    tarWinFlag = True
                    tarWinSt = i
            else:
                if alnedTarSeq[i] != '-':
                    if i-tarWinSt>=self.minSvLen:
                        opObj, opAlnSt, opAlnEd = 'tar', tarWinSt, i
                        break
                    else:
                        tarWinFlag = False
        retSeq = None
        if opObj is None:
            retSeq = oriSeq
        elif opObj == 'tar':
            retSeq = ''.join((alnedOriSeq[:opAlnSt]+alnedOriSeq[opAlnEd:]).split('-'))
        else:
            retSeq = ''.join((alnedOriSeq[:opAlnSt]+alnedTarSeq[opAlnSt:opAlnEd]+alnedOriSeq[opAlnEd:]).split('-'))
        return retSeq
