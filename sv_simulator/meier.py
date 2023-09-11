class meier:
    def __init__(self, seq, born_time, update_time=None):
        self.seq = seq
        self.born_time = born_time
        self.update_time = born_time if self.update_time is None \
            else update_time


