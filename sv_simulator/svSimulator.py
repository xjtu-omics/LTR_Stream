from multiprocessing import Pool
import random
random.seed(41)
import numpy as np

from meiPopulation import meiPopulation
from meier import meier
from varOperator import varOperator

class svSimulator:
    def __init__(self, maxPopSize, oriSeq, tarSeq, ranSvP, yearStep, cpP, cpuNum=48):
        self.maxPopSize = maxPopSize
        self.oriSeq = oriSeq
        self.tarSeq = tarSeq
        self.ranSvP = ranSvP
        self.yearStep = yearStep
        self.cpP = cpP
        self.cpuNum = cpuNum
        self.meiPop = meiPopulation(maxPopSize=self.maxPopSize)
        self.oper = varOperator(cpP = self.cpP)

    def initCpNumList(self):
        # return np.random.binomial(n=self.yearStep, p=self.cpP, size=self.meiPop.getPopSize())
        return np.random.binomial(n=self.yearStep, p=self.cpP, size=100)

    def simulate(self):
        oriTe = meier(id='0', seq=self.oriSeq, bornTime=0)
        self.meiPop.meiSet.add(oriTe)
        nowYear = self.yearStep
        while not self.meiPop.reachMaxPopSize():
            pool = Pool(processes=self.cpuNum)
            updatedTeList = []
            updatedMeiPop = meiPopulation(maxPopSize=self.maxPopSize)
            print(self.initCpNumList())
            exit()
            for te in self.meiPop:
                # print(te)
                # self.oper.geneSvedSeq(te, nowYear)
                pool.apply_async(self.oper.geneSvedSeq, (te, nowYear,))
                # updatedTeList.append(pool.apply_async(self.oper.geneSvedSeq, (te, nowYear, )))
            pool.close()
            pool.join()
            break
            # for teList in updatedTeList:
            #     for te in teList:
            #         updatedMeiPop.addMei(te)
            # self.meiPop = updatedMeiPop
            # nowYear += self.yearStep
