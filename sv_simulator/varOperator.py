import numpy as np
from copy import deepcopy
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


class varOperator:
    def __init__(self, cpP, tarSeq, yearStep, minSvLen, maxSvLen):
        self.cpP = cpP
        self.tarSeq = tarSeq
        self.yearStep = yearStep
        self.minSvLen = minSvLen
        self.maxSvLen = maxSvLen

    def adjustCpNumByYear(self, te, nowYear, cpNum):
        return int(cpNum/((nowYear-te.bornTime)/self.yearStep))

    def geneSvedSeq(self, te, nowYear, cpNum, dirSvPosList, dirSvLenList, svGenTyp):
        retList = []
        retList.append(te)
        cpNum = self.adjustCpNumByYear(te, nowYear, cpNum)
        for i in range(cpNum):
            nte = deepcopy(te)
            nte.id = f'{nte.id}_{te.cpedNum}'
            nte.cpedNum = 0
            nte.seq, nte.identity = self.dirSv(nte.seq, dirSvPosList[i], dirSvLenList[i], svGenTyp)
            nte.bornTime = nowYear
            retList.append(deepcopy(nte))
            te.addCpedNum(1)
        return retList

    def getAlnScore(self, oriSeq):
        # Dangerours to use for seq with lenth more than 10k
        aln = pairwise2.align.globalxs(oriSeq, self.tarSeq, -1, -0.2)[0]
        return aln.score

    def calIdentity(self, alnedSeq1, alnedSeq2, nomLen):
        identicalNum = 0
        for ca,cb in zip(alnedSeq1, alnedSeq2):
            if ca==cb:
                identicalNum+=1
        return identicalNum/nomLen

    def dirSv(self, oriSeq, dirSvPos, dirSvLen, svGenTyp):
        # 2023091815:19 Remove original global alignment that can be used for gap reconmendation
        # Only for substitutiom mode.
        if svGenTyp is None:
            raise 'Do not support gap reconmendation directed SV mode!'

        anchorLen = 1000
        toAlnEdPos = min(len(oriSeq), dirSvPos + dirSvLen + anchorLen)
        toAlnStPos = max(0, dirSvPos-anchorLen)
        preAnchorSize = min(anchorLen, dirSvPos)

        tarOriSeq = oriSeq[toAlnStPos:toAlnEdPos]
        aln = pairwise2.align.globalms(tarOriSeq, self.tarSeq, 2, -1, -100, -0.2)[0]
        alnedOriSeq, alnedTarSeq = aln.seqA, aln.seqB
        identity = self.calIdentity(alnedOriSeq, alnedTarSeq, toAlnEdPos-toAlnStPos)

        dirSvAlnPos = -1
        ooPos = 0
        for i in range(len(alnedOriSeq)):
            if ooPos>=preAnchorSize:
                dirSvAlnPos = i
                break
            if alnedOriSeq[i] != '-':
                ooPos += 1

        dirSvedSeq = oriSeq[:dirSvPos]
        for i in range(dirSvAlnPos, min(dirSvAlnPos+dirSvLen, len(alnedOriSeq))):
            if alnedTarSeq[i] != '-':
                dirSvedSeq += alnedTarSeq[i]
            if alnedOriSeq[i] != '-':
                ooPos += 1
        dirSvedSeq += oriSeq[(toAlnStPos+ooPos):]
        return dirSvedSeq, identity