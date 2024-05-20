import pandas as pd

import ttConfiger


class fileSysSupplier:
    def __init__(self, ltrParaFile):
        self.ltrParaFile = ltrParaFile
        paraDict = ttConfiger.getParaDict(ltrParaFile)
        self.baseD = paraDict['workDir']
        self.ltrFasta = paraDict['ltrFasta']
        try:
            self.configFile = paraDict['refConfigTsv']
        except:
            self.configFile = None
        self.refD = f'{self.baseD}/ref'
        self.danteD = f'{self.baseD}/dante'
        self.figureD = f'{self.baseD}/figure'
        try:
            self.refParaFile = paraDict['refConfigTsv']
            self.refData = pd.read_table(self.refParaFile, comment='#', sep='\t')
        except:
            self.refParaFile = None
            self.refData = None

        self.disMatCsvFile = f'{self.danteD}/toCalRest.distanceMat.csv'
        self.oriId2modSeqTabFile = f'{self.danteD}/finalModSeq.tab'
        self.modSeq2modIdTabFile = f'{self.danteD}/toCalRest.modSeq2modId.tab'
        self.totalInfoCsvFile = f'{self.danteD}/toCalRest.totoalInfo.csv'
        self.ltrRetrieverD = f'{self.baseD}/ltrRetriever'
        self.disMatCsv = f'{self.danteD}/toCalRest.distanceMat.csv'

        self.classTsv = f'{self.figureD}/classInfo.tsv'
        self.d2InsPdf = f'{self.figureD}/2D_insTime.pdf'
    def getSpeId2faiFile(self):
        speId2faiFile = {}
        for ind, row in self.refData.iterrows():
            speId = row['speId']
            speId2faiFile[speId] = f'{self.refD}/{speId}.ref.fa.fai'
        return speId2faiFile
    def loadDisMat(self):
        return pd.read_csv(self.disMatCsv, header=None)
    def __eq__(self, other):
        return self.baseD == other.baseD