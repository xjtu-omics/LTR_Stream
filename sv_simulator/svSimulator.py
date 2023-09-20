from multiprocessing import Pool
import random
import numpy as np
import sys

from meiPopulation import meiPopulation
from meier import meier
from varOperator import varOperator

class svSimulator:
    def __init__(self, maxPopSize, oriSeq, oriSeqName, tarSeq, ranSvP, yearStep, cpP, cpuNum, minSvLen, maxSvLen, dirSvTypNum=2):
        self.maxPopSize = maxPopSize
        self.oriSeq = oriSeq
        self.oriSeqName = oriSeqName
        self.tarSeq = tarSeq
        self.ranSvP = ranSvP
        self.yearStep = yearStep
        self.cpP = cpP
        self.cpuNum = cpuNum
        self.dirSvTypNum = dirSvTypNum
        self.minSvLen = minSvLen
        self.maxSvLen = maxSvLen
        self.meiPop = meiPopulation(maxPopSize=self.maxPopSize)
        self.oper = varOperator(cpP = self.cpP, tarSeq=self.tarSeq, yearStep=self.yearStep, minSvLen=self.minSvLen, maxSvLen=maxSvLen)
    def initCpNumList(self):
        return np.random.binomial(n=self.yearStep, p=self.cpP, size=self.meiPop.getPopSize())
    def initDirSvTypList(self, cpNumList):
        dirSvTypList = []
        for cpNum in cpNumList:
            dirSvTypList.append(np.random.randint(self.dirSvTypNum, size=cpNum))
        return dirSvTypList
    def initDirSvPosList(self, cpNumList):
        dirSvPosList = []
        for cpNum, te in zip(cpNumList, self.meiPop):
            dirSvPosList.append(np.random.randint(len(te.seq), size=cpNum))
        return dirSvPosList

    def initDirSvLenList(self, cpNumList, svGenTyp):
        dirSvLenList = []
        if svGenTyp=='sub':
            for cpNum in cpNumList:
                dirSvLenList.append(np.random.randint(self.minSvLen, self.maxSvLen, size=cpNum))
        else:
            for cpNum in cpNumList:
                dirSvLenList.append([ None for i in range(cpNum)])
        return dirSvLenList


    def simulate(self, svGenTyp):
        # np.random.seed(51)
        oriTe = meier(id=self.oriSeqName, seq=self.oriSeq, bornTime=0)
        self.meiPop.meiSet.add(oriTe)
        nowYear = self.yearStep
        while not self.meiPop.reachMaxPopSize():
            updatedTeList = []
            updatedMeiPop = meiPopulation(maxPopSize=self.maxPopSize)
            cpNumList = self.initCpNumList()
            dirSvTypList = self.initDirSvTypList(cpNumList)
            dirSvPosList = self.initDirSvPosList(cpNumList)
            dirSvLenList = self.initDirSvLenList(cpNumList, svGenTyp)

            pool = Pool(processes=self.cpuNum)
            for te, cpNum, dirSvPoss, dirSvLens in zip(self.meiPop, cpNumList, dirSvPosList, dirSvLenList):
                updatedTeList.append(pool.apply_async(self.oper.geneSvedSeq, (te, nowYear, cpNum, dirSvPoss, dirSvLens, svGenTyp)))
            pool.close()
            pool.join()

            for teList in updatedTeList:
                for te in teList.get():
                    updatedMeiPop.addMei(te)
            self.meiPop = updatedMeiPop
            print('PopSize:', self.meiPop.getPopSize())
            print('Mean Test Identity:', self.meiPop.getMeanIdentity(nowYear))
            nowYear += self.yearStep
