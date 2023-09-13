from meiPopulation import meiPopulation
class svSimulator:
    def __init__(self, maxPopSize, oriSeq, tarSeq, ranSvP, yearStep, cpP):
        self.meiPop = meiPopulation(maxPopSize=maxPopSize)
        self.maxPopSize = maxPopSize
        self.oriSeq = oriSeq
        self.tarSeq = tarSeq
        self.ranSvP = ranSvP
        self.yearStep = yearStep
        self.cpP = cpP
    def simulate(self):

