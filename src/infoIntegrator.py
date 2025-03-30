import pandas as pd
import numpy as np
import ttUtils
import pysam
from copy import deepcopy
import ttDb
from fileSysSupplier import fileSysSupplier



# Do not change infoIntegrator.infoDf by operator =,
# only use setInfoDf() to change its value.
# Or, there will be problems with filter
class infoIntegrator:
    def __init__(self, infoDf, fss, randomState=42):
        self.infoDf = infoDf
        self.fss = fss
        self.randomState = randomState
        self.oriInfoDf = self.infoDf.copy(deep=True)
    def setInfoDf(self, newInfoDf=None):
        if newInfoDf is None:
            self.oriInfoDf = self.infoDf.copy(deep=True)
            return
        self.infoDf = newInfoDf.copy(deep=True)
        self.oriInfoDf = newInfoDf.copy(deep=True)
    @classmethod
    def init_Null(cls, ltrParaFile):
        fss = fileSysSupplier(ltrParaFile)
        infoDf = pd.DataFrame(dict(
            oriId = [],
            modId = [],
            Class = [],
            isTerminalCluster = [],
            zoomIn = [],
            x = [],
            y = [],
            z = [],
            xx = [],
            yy = []
        ))
        infoDf = cls.reformatInfoDf(infoDf)
        return infoIntegrator(infoDf, fss)
    @classmethod
    def init_from_adata(cls, adata, ltrParaFile):
        # addta should be subAdata of zoomIn
        fss = fileSysSupplier(ltrParaFile)
        zoomInPrefix = adata.uns['prefix']

        # First part, those modId that used to construct trajectories
        modIdList = list(adata.obs['modId'])
        cluIdList = list(adata.obs['cluId'])
        infoDf = pd.DataFrame(dict(modId=modIdList, cluId=cluIdList))
        posDf = pd.DataFrame(adata.obsm['X_dr'], columns=['x', 'y', 'z'])
        infoDf = pd.concat([infoDf, posDf], axis=1)
        infoDf = infoDf[ infoDf['cluId']!=0 ]
        classIdList = [ f'{zoomInPrefix}_{chr(ord("A")+int(cluId)-1)}' for cluId in list(infoDf['cluId'])]
        zoomInList = [zoomInPrefix for cluId in list(infoDf['cluId'])]
        isTerList = [True for cluId in list(infoDf['cluId'])]

        infoDf = pd.concat([infoDf,
                            pd.DataFrame(dict(Class=classIdList,
                                              zoomIn=zoomInList,
                                              isTerminalCluster=isTerList), index=infoDf.index)],
                           axis=1)

        # Second part, those modId that could be further zoomIn
        modIdList = list(adata.uns['oToPltPointPosDf'].index)
        cluIdList = list(adata.uns['oToPltPointPosDf']['cluId'])
        infoDf2 = pd.DataFrame(dict(modId=modIdList, cluId=cluIdList))
        infoDf2 = infoDf2[ infoDf2['cluId']!=0 ]
        classIdList = [ f'{zoomInPrefix}_{chr(ord("A")+int(cluId)-1)}' for cluId in list(infoDf2['cluId'])]
        zoomInList = [zoomInPrefix for cluId in list(infoDf2['cluId'])]
        isTerList = [False for cluId in list(infoDf2['cluId'])]
        infoDf2 = pd.concat([infoDf2,
                            pd.DataFrame(dict(Class=classIdList,
                                              zoomIn=zoomInList,
                                              isTerminalCluster=isTerList), index=infoDf2.index)],
                           axis=1)
        # Merge the two dataframes
        infoDf = pd.concat([infoDf, infoDf2], axis=0)
        infoDf = infoDf.reset_index(drop=True)
        infoDf = cls.expandInfoDfOriId(infoDf, fss)
        infoDf = cls.reformatInfoDf(infoDf)
        return infoIntegrator(infoDf, fss)
    @classmethod
    def init_from_posDf(cls, posDf, posDf_2d, zoomInPrefix, ltrParaFile, cluId2isTer):
        fss = fileSysSupplier(ltrParaFile)
        modIdList = list(posDf.index)
        cluIdList = list(posDf['cluId'])
        infoDf = pd.DataFrame(dict(modId=modIdList, cluId=cluIdList))
        infoDf = pd.concat([infoDf, posDf[['x', 'y', 'z']].reset_index(drop=True)], axis=1)
        tmpDf = posDf_2d.copy(deep=True)
        tmpDf.columns = ["xx", "yy"]
        infoDf = pd.concat([infoDf, tmpDf[['xx', 'yy']].reset_index(drop=True)], axis=1)
        infoDf = infoDf[ infoDf['cluId']!=0 ]
        classIdList = [ f'{zoomInPrefix}_{chr(ord("A")+int(cluId)-1)}' for cluId in list(infoDf['cluId']) ]
        zoomInList = [zoomInPrefix for cluId in list(infoDf['cluId'])]
        isTerList = [cluId2isTer[cluId] for cluId in list(infoDf['cluId'])]
        infoDf = pd.concat([infoDf,
                            pd.DataFrame(dict(Class=classIdList,
                                              zoomIn=zoomInList,
                                              isTerminalCluster=isTerList), index=infoDf.index)],
                           axis=1)
        infoDf = infoDf.reset_index(drop=True)
        infoDf = cls.expandInfoDfOriId(infoDf, fss)
        infoDf = cls.reformatInfoDf(infoDf)
        return infoIntegrator(infoDf, fss)
    @classmethod
    def init_fromTsv(cls, tsvFile, fss=None):
        infoDf = pd.read_csv(tsvFile, sep='\t')
        return infoIntegrator(infoDf, fss)
    @classmethod
    def expandInfoDfOriId(cls, infoDf, fss):
        modId2oriIdSet = ttUtils.getModId2OriIdSet(fss.oriId2modSeqTabFile, fss.modSeq2modIdTabFile)
        oriIdSetList = [modId2oriIdSet[modId] for modId in list(infoDf['modId'])]
        infoDf = pd.concat([infoDf,
                            pd.DataFrame(dict(oriId=oriIdSetList), index=infoDf.index)],
                           axis=1)
        infoDf = infoDf.explode('oriId')
        return infoDf
    @classmethod
    def reformatInfoDf(cls, infoDf):
        infoDf = infoDf[['oriId', 'modId', 'Class', 'isTerminalCluster', 'zoomIn', 'x', 'y', 'z', 'xx', 'yy']]
        infoDf = infoDf.astype({
            'modId':int,
            'isTerminalCluster':bool
        })
        infoDf = infoDf.reset_index(drop=True)
        return infoDf
    @classmethod
    def merge(cls, inte1, inte2):
        infoDf = pd.concat([inte1.infoDf, inte2.infoDf], axis=0)
        infoDf = infoDf.reset_index(drop=True)
        if (inte1.fss is None) or (inte2.fss is None) or (inte1.fss==inte2.fss):
            mergedFss = inte2.fss if inte1.fss is None else inte1.fss
            return infoIntegrator(infoDf, mergedFss)
        else:
            assert 'infoIntegrator with different fss!'
    def addValiScoreDf(self, scoreTsv):
        pass
    def isNull(self):
        return self.infoDf.shape[0] == 0
    def addZoomInLevel(self):
        print(self.infoDf.zoomIn.drop_duplicates())
        zoomInLevelList = [len(zoomInId.split('_')) for zoomInId in list(self.infoDf['zoomIn'])]
        zoomInLevelSer = pd.Series(zoomInLevelList, index=self.infoDf.index, dtype=int)
        self.infoDf['zoomInLevel'] = zoomInLevelSer
        self.setInfoDf()
    def addCreditLabel(self):
        def getCreditCutOffVal(infoDf):
            infoDf = deepcopy(infoDf)
            if not infoDf.isTerminalCluster.iloc[0]:
                return None
            else:
                disSer = infoDf.apply(
                    lambda row: np.sqrt(row.x*row.x + row.y*row.y + row.z*row.z),
                    axis=1
                )
                return 0.1*disSer.quantile(.95)
        def getCreditLabel(cutoff, x, y, z):
            if cutoff is None:
                return 'Not_final_zoomIn'
            dis = np.sqrt(x*x + y*y + z*z)
            if dis>=cutoff:
                return 'Confident'
            else:
                return 'Low_confidence'
        df = deepcopy(self.infoDf)
        cutOffSer = df.groupby('zoomIn').apply(
            getCreditCutOffVal
        )
        zoomIn2cutOff = cutOffSer.to_dict()
        labelSer = df.apply(
            lambda row: getCreditLabel(
                zoomIn2cutOff[row.zoomIn], row.x, row.y, row.z
            ),
            axis=1
        )
        self.infoDf['creditLabel'] = labelSer

    def writeRel(self, outTsv='classInfo.tsv'):
        self.addZoomInLevel()
        self.addCreditLabel()
        outPath = f'{self.fss.figureD}/{outTsv}'
        self.infoDf.to_csv(outPath, sep='\t', index=False)
    def getZoomInList(self):
        return list(self.infoDf['zoomIn'].drop_duplicates())
    '''
    def getKNeibourPerc(self, num=.1):
        # Do not use this function!
        # Not Complete!
        pandarallel.initialize(nb_workers=48, )
        def calOulaDis(ser, tarSer):
            return np.linalg.norm(ser[['x', 'y', 'z']] - tarSer)

        def calSameClusterPerc(ser, infoDf):
            infoDf['oulaDis'] = infoDf.apply(calOulaDis, tarSer=ser[['x', 'y', 'z']], axis=1)

        infoDf = self.infoDf.copy(deep=True)
        self.infoDf.parallel_apply(calSameClusterPerc, infoDf=infoDf, axis=1)
    '''
    def filterCentPoint(self, perc=20):
        def filt(infoDf, perc):
            cordDis = infoDf[['x', 'y', 'z']].apply(np.linalg.norm, axis=1)
            cordDisPercValue = cordDis.quantile(q=.99)*(perc/100)
            return infoDf[cordDis>cordDisPercValue]
        # print(self.oriInfoDf.shape[0])
        self.infoDf = self.oriInfoDf.groupby('Class').apply(filt, perc=perc).reset_index(drop=True)
        # print(self.infoDf.shape[0])
    def getSeqFromFasta(self, oriId, fasta):
        return fasta.fetch(oriId)
    def addNuclSeq(self, faFile):
        with pysam.FastaFile(faFile) as fasta:
            self.infoDf['seq'] = self.infoDf.oriId.apply(self.getSeqFromFasta, fasta=fasta)
            self.setInfoDf()
    def ranNucl(self, num, outFa, Classes=None, zoomInLevel=None):
        classDf = None
        if (not (Classes is None)) or (not (zoomInLevel is None)):
            classDf = self.filterDf(Classes, zoomInLevel)
        else:
            classDf = self.infoDf[self.infoDf.isTerminalCluster==True]
        classDf['Class'] = classDf['finalClass']
        nuclDf = classDf.groupby('Class').apply(
            lambda x: x.sample(n=min(num, x.shape[0]), replace=False, random_state=self.randomState))
        nuclDf = nuclDf.reset_index(drop=True)
        with open(outFa, 'w') as of:
            for ind, row in nuclDf.iterrows():
                print(f'>{row.Class}__{ind}', file=of)
                print(row.seq, file=of)
        return outFa
    def geneFinalClassLabel(self):
        al = 'A'
        classSer = self.infoDf[self.infoDf.isTerminalCluster==True]['Class'].drop_duplicates().reset_index(drop=True)
        transDict = {Class:chr(ord(al)+ind) for ind, Class in  classSer.items()}
        self.infoDf['finalClass'] = self.infoDf.Class.map(transDict)
        self.infoDf.finalClass = self.infoDf.finalClass.fillna('unClustered')
        self.setInfoDf()
    def filterDf(self, Classes=None, zoomInLevel=None):
        classDf = self.infoDf.copy(deep=True)
        if Classes is not None:
            classDf = classDf[classDf.Class.isin(Classes)]
        if zoomInLevel is not None:
            classDf = classDf[classDf.zoomInLevel == zoomInLevel]
        return classDf
    def addUserLabel(self, prefix):
        finalClassSer = self.infoDf.Class[self.infoDf.isTerminalCluster==True].drop_duplicates()
        finalClassSer = finalClassSer.reset_index(drop=True)
        al = 'A'
        Class2outLabel = {Class:f'{prefix}_{chr(ord(al)+ind)}' for ind, Class in finalClassSer.items()}
        self.infoDf['Cluster'] = self.infoDf.Class.map(Class2outLabel)
        self.setInfoDf()
    def writeUserRel(self, prefix, outTsv='clusterRel.tsv'):
        self.addUserLabel(prefix)
        userInfoDf = self.infoDf.copy(deep=True)
        userInfoDf = userInfoDf[userInfoDf.isTerminalCluster==True]
        userInfoDf = userInfoDf[['oriId', 'Cluster']]
        outPath = f'{self.fss.figureD}/{outTsv}'
        userInfoDf.to_csv(outPath, sep='\t', index=False)
    def reformatZoomInOrder(self, prefix, zoomInAdataDict):
        # def updateTransDict(tmpDf, transDict):
        #     tmpDf = tmpDf.sort_values(['isTerminalCluster', 'Class'])
        #     aOrder = tmpDf.Class.drop_duplicates().tolist()
        #     bOrder = sorted(aOrder)
        #     for a,b in zip(aOrder, bOrder):
        #         transDict[b] = a
        # self.infoDf.groupby('zoomIn').apply(updateTransDict, transDict)
        # transDict[prefix] = prefix
        # self.infoDf.Class = self.infoDf.Class.map(transDict)
        # self.infoDf.zoomIn = self.infoDf.zoomIn.map(transDict)
        # newDict = {transDict[zoomIn]:adata for zoomIn,adata in zoomInAdataDict.items()}
        # return newDict
        def updateTransDict(tmpDf, transDict):
            tmpDf = tmpDf.sort_values(by=['isTerminalCluster', 'Class'])
            aOrder = tmpDf.Class.drop_duplicates().tolist()
            bOrder = sorted(aOrder)
            for a, b in zip(aOrder, bOrder):
                transDict[a] = b
        def updatePrefix(ser, transDict, zoomInLevel):
            if ser.zoomInLevel<zoomInLevel:
                return ser
            for a,b in transDict.items():
                if ser.Class.startswith(a):
                    ser.Class = b+ser.Class.split(a)[1]
                    if ser.zoomIn.startswith(a):
                        ser.zoomIn = b+ser.zoomIn.split(a)[1]
                    break
            return ser
        def check(ser):
            return ser.Class.startswith(ser.zoomIn)
        self.addZoomInLevel()
        oldInfoDf = self.infoDf.copy(deep=True)
        zls = self.infoDf.zoomInLevel.drop_duplicates().sort_values().tolist()
        for zl in zls:
            transDict = {}
            self.infoDf[self.infoDf.zoomInLevel==zl].groupby('zoomIn').apply(updateTransDict, transDict=transDict)
            self.infoDf = self.infoDf.apply(updatePrefix, transDict=transDict, zoomInLevel=zl, axis=1)
        self.setInfoDf()
        transDict = {}
        for (ind1, row1), (ind2, row2) in zip(oldInfoDf.iterrows(), self.infoDf.iterrows()):
            transDict[row1.zoomIn] = row2.zoomIn
        print(transDict)
        newDict = {transDict[zoomIn]:adata for zoomIn,adata in zoomInAdataDict.items()}
        print('check:', np.sum(self.infoDf.apply(check, axis=1)))
        return newDict
    def addInsTime(self, fss, inId):
        oriId2insTime = ttDb.getOriId2insertTime(fss.configFile, fss.ltrRetrieverD, None)
        if inId:
            inRangeId2oriId = ttDb.getInRangeId2oriId(fss.configFile, fss.ltrRetrieverD)
            inRangeId2insTime = {inRangeId:oriId2insTime[oriId] for inRangeId, oriId in inRangeId2oriId.items()}
            oriId2insTime = inRangeId2insTime
        self.infoDf['insTime'] = self.infoDf.oriId.map(oriId2insTime)
        self.setInfoDf()
    def getClass2finalClass(self):
        tmpDf = self.infoDf[self.infoDf.isTerminalCluster==True]
        relDict = {}
        for ind, row in tmpDf.iterrows():
            relDict[row.Class] = row.finalClass
        return relDict
    def addPosInfo(self):
        from ttUtils import seqName2bainfo
        self.infoDf['chr'] = self.infoDf.oriId.apply(
            lambda x:seqName2bainfo(x).chr
        )
        self.infoDf['st'] = self.infoDf.oriId.apply(
            lambda x:seqName2bainfo(x).st
        )
        self.infoDf['ed'] = self.infoDf.oriId.apply(
            lambda x:seqName2bainfo(x).ed
        )
        self.setInfoDf()
    def getClass2finalClass(self):
        class2finalClass = {}
        for ind, row in self.infoDf.iterrows():
            class2finalClass[row.Class] = row.finalClass
        return class2finalClass





























































































