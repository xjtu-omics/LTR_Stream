from multiprocessing import Pool
import random
import numpy as np
import sys

from meiPopulation import meiPopulation
from meier import meier
from varOperator import varOperator

class svSimulator:
    def __init__(self, maxPopSize, oriSeq, tarSeq, ranSvP, yearStep, cpP, cpuNum, minSvLen, dirSvTypNum=2):
        self.maxPopSize = maxPopSize
        self.oriSeq = oriSeq
        self.tarSeq = tarSeq
        self.ranSvP = ranSvP
        self.yearStep = yearStep
        self.cpP = cpP
        self.cpuNum = cpuNum
        self.dirSvTypNum = dirSvTypNum
        self.minSvLen = minSvLen
        self.meiPop = meiPopulation(maxPopSize=self.maxPopSize)
        self.oper = varOperator(cpP = self.cpP, tarSeq=self.tarSeq, yearStep=self.yearStep, minSvLen=self.minSvLen)
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
    def simulate(self):
        np.random.seed(31)
        oriTe = meier(id='0', seq=self.oriSeq, bornTime=0)
        self.meiPop.meiSet.add(oriTe)
        nowYear = self.yearStep
        while not self.meiPop.reachMaxPopSize():
            updatedTeList = []
            updatedMeiPop = meiPopulation(maxPopSize=self.maxPopSize)
            cpNumList = self.initCpNumList()
            dirSvTypList = self.initDirSvTypList(cpNumList)
            dirSvPosList = self.initDirSvPosList(cpNumList)

            pool = Pool(processes=self.cpuNum)
            alnScoreRetList = []
            for te in self.meiPop:
                alnScoreRetList.append(pool.apply_async(self.oper.getAlnScore, (te.seq, )))
            pool.close()
            pool.join()
            alnScoreList = [ret.get() for ret in alnScoreRetList]
            print('Max Score:', np.max(alnScoreList))

            pool = Pool(processes=self.cpuNum)
            for te, cpNum, dirSvTyps, dirSvPoss in zip(self.meiPop, cpNumList, dirSvTypList, dirSvPosList):
                updatedTeList.append(pool.apply_async(self.oper.geneSvedSeq, (te, nowYear, cpNum, dirSvTyps, dirSvPoss)))
            pool.close()
            pool.join()

            for teList in updatedTeList:
                for te in teList.get():
                    updatedMeiPop.addMei(te)
            self.meiPop = updatedMeiPop
            nowYear += self.yearStep
            print('PopSize:', self.meiPop.getPopSize())
