class meier:
    def __init__(self, seq, bornTime, updateTime=None):
        self.seq = seq
        self.bornTime = bornTime
        self.updateTime = bornTime if self.updateTime is None \
            else updateTime


