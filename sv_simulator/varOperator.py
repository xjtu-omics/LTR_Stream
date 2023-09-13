import numpy as np
from copy import deepcopy

class varOperator:
    def __init__(self, cpP):
        self.cpP = cpP


    def geneSvedSeq(self, te, nowYear, cpNum):
        retList = []
        for i in range(cpNum):
            retList.append(deepcopy(te))
        # print(te)
        # a = np.random.binomial(n=1e5, p=self.cpP, size=1)
        # print(a[0])
