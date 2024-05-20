import pysam
import pandas as pd
import numpy as np
from itertools import combinations
from scipy.stats import wilcoxon, ranksums

from fileSysSupplier import fileSysSupplier

from validation.p2pValidation import p2pValidator
from validation.protDomainValidation import crossChangeValidator
class clusterValidator:
    def __init__(self, ii, ltrParaFile, filtCentPerc=20):
        self.ii = ii
        self.ii.filterCentPoint(perc=filtCentPerc)
        self.fss = fileSysSupplier(ltrParaFile)
    @classmethod
    def fq2df(cls, fastq, sortNames=True):
        names = []
        seqs = []
        with pysam.FastxFile(fastq) as inf:
            for rec in inf:
                names.append(rec.name)
                seqs.append(rec.sequence)
        relDf = pd.DataFrame(dict(name=names, seq=seqs))
        if sortNames:
            relDf = relDf.sort_values(by='name')
        return relDf
    def p2pAln(self, ranFa, cpuNum):
        df = self.fq2df(ranFa)
        df.columns = ['Class', 'seq']
        pv = p2pValidator(df, cpuNum=cpuNum)
        relList = pv.getAll2allIdentity()
        relList = [(df.Class.iloc[rel[0]], df.Class.iloc[rel[1]], rel[2]) for rel in relList]
        return crossChangeValidator.p2pAlnRelList2Mat(relList)
    def validateNucl(self, prefix='nucl', num=20, Classes=None, zoomInLevel=None, cpuNum=48):
        outTsv = f'{self.fss.figureD}/clusterNuclVali.tsv'
        ranFa = f'{self.fss.danteD}/validate.nucl.{prefix}.fa'
        self.ii.addNuclSeq(self.fss.ltrFasta)
        self.ii.ranNucl(num, ranFa, Classes, zoomInLevel)
        p2pIdenMat = self.p2pAln(ranFa, cpuNum)

        p2pIdenMat.to_pickle(f'{self.fss.figureD}/clusterNuclIdenMat.pkl')
        df = self.calCluDiffScore(p2pIdenMat, zoomInLevel=zoomInLevel, Classes=Classes)
        df.to_csv(outTsv, sep='\t', index=None)
    def calCluDiffScore(self, p2pIdenMat, Classes=None, zoomInLevel=None, onlyMinimum=True):
        def calEachClassDiffScore(tarClass, idenMat, toValClassSer, onlyMinimum):
            Class2idens = {Class:[] for ind,Class in toValClassSer.items()}
            for ca, cb in combinations(list(idenMat.index), 2):
                iden = idenMat.loc[ca,cb]
                ca = ca.split('__')[0]
                cb = cb.split('__')[0]
                if ca == tarClass or cb == tarClass:
                    if ca != tarClass:
                        ca, cb = cb, ca
                    Class2idens[cb].append(iden)
            Class2mean = {Class:np.mean(idens) for Class, idens in Class2idens.items()}
            selfMean = Class2mean[tarClass]
            Class2mean.pop(tarClass)
            otherMeans = [Class2mean[Class] for Class in Class2mean if Class != tarClass]
            if onlyMinimum:
                foldChange = min(otherMeans)/selfMean
                stat, pv = ranksums(Class2idens[min(Class2mean)], Class2idens[tarClass], alternative='greater')
                return pd.Series([foldChange, pv], index=['foldChange', 'p-value'])
            else:
                tarClasses = []
                foldChanges = otherMeans/selfMean
                stats, pvs, Classes = [], [], []
                for Class in Class2mean:
                    Classes.append(Class)
                    stat, pv = ranksums(Class2idens[Class], Class2idens[tarClass], alternative='greater')
                    stats.append(stat)
                    pvs.append(pv)
                    tarClasses.append(tarClass)
                return pd.DataFrame(dict(tarClass=tarClasses,
                                         otherClass=Classes,
                                         stat=stats,
                                         pv=pvs,
                                         foldChange=foldChanges))

        subDf = None
        if (not (Classes is None)) or (not (zoomInLevel is None)):
            subDf = self.ii.filterDf(Classes, zoomInLevel)
        else:
            subDf = self.ii.infoDf[self.ii.infoDf.isTerminalCluster==True]
        toValClassSer = subDf.Class.drop_duplicates()
        diffScoreDf = toValClassSer.apply(calEachClassDiffScore,
                                          idenMat=p2pIdenMat,
                                          toValClassSer=toValClassSer,
                                          onlyMinimum=onlyMinimum)
        diffScoreDf.index = list(toValClassSer)
        if onlyMinimum:
            return diffScoreDf
        else:
            relDf = pd.DataFrame(columns=list(diffScoreDf.iloc[0].columns))
            for ind, df in diffScoreDf.items():
                relDf = pd.concat([relDf, df], axis=0)
            return relDf