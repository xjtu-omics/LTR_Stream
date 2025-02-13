from Bio import Align
from Bio.Seq import Seq
from itertools import combinations
from joblib import Parallel, delayed
class p2pValidator:
    def __init__(self, toValidateDf, cpuNum=4):
        # myDf should have two cols ["Class", "seq"]
        self.myDf = toValidateDf
        self.cpuNum = cpuNum
    @classmethod
    def alignSeqAndGetIdentity(cls, seqA, seqB):
        aligner = Align.PairwiseAligner()
        alnedA, alnedB = aligner.align(seqA, seqB)[0]
        id1 = cls.calculate_identity(alnedA, alnedB)
        seqC = str(Seq(seqB).reverse_complement())
        alnedA, alnedC = aligner.align(seqA, seqC)[0]
        id2 = cls.calculate_identity(alnedA, alnedC)
        return max(id1, id2)
    @classmethod
    def calculate_identity(cls, sequenceA, sequenceB):
        """
        Returns the percentage of identical characters between two sequences.
        Assumes the sequences are aligned.
        """
        sa, sb, sl = sequenceA, sequenceB, len(sequenceA)
        matches = [sa[i] == sb[i] for i in range(sl)]
        seq_id = (100 * sum(matches)) / sl
        # gapless_sl = sum([1 for i in range(sl) if (sa[i] != '-' and sb[i] != '-')])
        # gap_id = (100 * sum(matches)) / gapless_sl
        return seq_id
    def getAll2allIdentity(self, typ=None):

        def runBioPairAln(paraDict):
            aligner = Align.PairwiseAligner()
            alignment = aligner.align(paraDict['seqI'], paraDict['seqJ'])[0]
            # print(alignment)
            alnedA, alnedB = alignment
            id1 = self.calculate_identity(alnedA, alnedB)
            if ('typ' not in paraDict) or (paraDict['typ']=='nucl'):
                seqJ = str(Seq(paraDict['seqJ']).reverse_complement())
                alignment = aligner.align(paraDict['seqI'], seqJ)[0]
                alnedA, alnedB = alignment
                id2 = self.calculate_identity(alnedA, alnedB)
                return paraDict['i'], paraDict['j'], max(id1, id2)
            elif paraDict['typ'] == 'prot':
                return paraDict['i'], paraDict['j'], id1

        seqNum = self.myDf.shape[0]
        toAlnDictList = []
        for i,j in combinations(range(seqNum), 2):
            tmpParaDict = dict(i=i, j=j,
                               seqI=self.myDf.seq.iloc[i],
                               seqJ=self.myDf.seq.iloc[j])
            if typ is not None:
                tmpParaDict['typ'] = typ
            toAlnDictList.append(tmpParaDict)
        # relList = [runBioPairAln(alnDict) for alnDict in toAlnDictList]
        # return relList
        # print(seqNum)
        # print(len(toAlnDictList))
        # toAlnDictList = toAlnDictList[:1000]
        relList = Parallel(n_jobs=self.cpuNum)(
            delayed(runBioPairAln)(alnDict) for alnDict in toAlnDictList
        )
        return relList
    def ppline_plotIdentityBetweenGroups(self):
        self.getAll2allIdentity()