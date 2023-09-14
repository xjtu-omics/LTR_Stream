class meier:
    def __init__(self, id, seq, bornTime, updateTime=None):
        self.id = id
        self.seq = seq
        self.bornTime = bornTime
        self.updateTime = bornTime if updateTime is None \
            else updateTime
        self.cpedNum = 0

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        return self.id == other.id

    def addCpedNum(self, num):
        self.cpedNum += num


