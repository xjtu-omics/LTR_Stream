class meiPopulation:
    def __init__(self, maxPopSize):
        self.maxPopSize = maxPopSize
        self.meiSet = set()

    def addMei(self, mei):
        self.meiSet.add(mei)

    def __iter__(self):
        teList = [te for te in self.meiSet]
        sortedTeList = sorted(teList, key=lambda x:x.id)
        for te in sortedTeList:
            yield te

    def reachMaxPopSize(self):
        return len(self.meiSet)>=self.maxPopSize

    def getPopSize(self):
        return len(self.meiSet)
