class meiPopulation:
    def __init__(self, maxPopSize):
        self.maxPopSize = maxPopSize
        self.meiSet = set()

    def addMei(self, mei):
        self.meiSet.add(mei)

    def __iter__(self):
        i = 0
        while i<100:
            yield i
            i += 1

    def reachMaxPopSize(self):
        return False

    def getPopSize(self):
        return len(self.meiSet)
