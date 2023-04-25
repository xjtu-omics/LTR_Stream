import sys
import yaml
import pandas as pd
with open('envConfig.yaml', 'r') as of:
    ttPara = yaml.safe_load(of)
    sys.path.append(ttPara['LTR_Stream'])

import ttConfiger


class fileSysSupplier:
    def __init__(self, ltrParaFile):
        paraDict = ttConfiger.getParaDict(ltrParaFile)
        self.baseD = paraDict['workDir']
        self.refD = f'{self.baseD}/ref'
        self.refParaFile = paraDict['refConfig']
        self.refData = pd.read_table(self.refParaFile, comment='#', sep='\t')

    def getSpeId2faiFile(self):
        speId2faiFile = {}
        for ind, row in self.refData.iterrows():
            speId = row['speId']
            speId2faiFile[speId] = f'{self.refD}/{speId}.ref.fa.fai'
        return speId2faiFile
