import pandas as pd

def getParaDict(configFile):
    paraDict = {}
    paraD = pd.read_table(configFile,header=None,sep='\t',comment='#')
    for ind,row in paraD.iterrows():
        paraDict[row[0]] = row[1]
    return paraDict

class ttConfiger:
    def __init__(self, ltrParaFile):
        tmpParaDict = getParaDict(ltrParaFile)
        self.baseD = tmpParaDict['workDir']
        self.ltrFa = tmpParaDict['ltrFasta']
        try:
            self.minOver = float(tmpParaDict['minOverLapForNovelModule'])
        except:
            self.minOver = 0.8
        try:
            self.topModNum = int(tmpParaDict['topModNum'])
        except:
            self.topModNum = 250
        try:
            self.epgLambda = float(tmpParaDict['epgLambda'])
        except:
            self.epgLambda = 0.002
        try:
            self.epgAlpha = float(tmpParaDict['epgAlpha'])
        except:
            self.epgAlpha = 0.2
        try:
            self.epgMu = float(tmpParaDict['epgMu'])
        except:
            self.epgMu = 0.02
        try:
            self.blastEvalue = tmpParaDict['blastEvalue']
        except:
            self.blastEvalue = '1e-10'
        try:
            self.tsneLearningRate = float(tmpParaDict['tsneLearningRate'])
        except:
            self.tsneLearningRate = 6
        try:
            self.tsnePerplexity = float(tmpParaDict['tsnePerplexity'])
        except:
            self.tsnePerplexity = 200
        try:
            self.tsneEarlyExaggeration = float(tmpParaDict['tsneEarlyExaggeration'])
        except:
            self.tsneEarlyExaggeration = 6
        try:
            self.outPrefix = tmpParaDict['outPrefix']
        except:
            self.outPrefix = 'ltrStream'
        try:
            self.refConfig = tmpParaDict['refConfigTsv']
        except:
            self.refConfig = None
        try:
            self.cluCentCut = tmpParaDict['cluCentCut']
        except:
            self.cluCentCut = 0.1
        try:
            self.maxZoomInLevel = float(tmpParaDict['maxZoomInLevel'])
        except:
            self.maxZoomInLevel = -1
        try:
            self.tsneLearningRate = int(tmpParaDict['tsneLearningRate'])
        except:
            self.tsneLearningRate = 7
        try:
            self.tsnePerplexity = float(tmpParaDict['tsnePerplexity'])
        except:
            self.tsnePerplexity = 100
        try:
            self.tsneEarlyExaggeration = float(tmpParaDict['tsneEarlyExaggeration'])
        except:
            self.tsneEarlyExaggeration = 6
        self.tarSpeStrSet = self.getTarSpeStrSet(tmpParaDict)
    def getTarSpeStrSet(self, tmpParaDict):
        tarSpeStrSet = set()
        if 'tarSpeIds' not in tmpParaDict:
            return None
        else:
            for group in tmpParaDict['tarSpeIds'].split(';'):
                tarSpeStrSet.add('&'.join(sorted(group.split(','))))
        return tarSpeStrSet
