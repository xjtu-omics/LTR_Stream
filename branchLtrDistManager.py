import sys
import pandas as pd
from sklearn.cluster import KMeans

import yaml

with open('envConfig.yaml', 'r') as of:
    ttPara = yaml.safe_load(of)
    sys.path.append(ttPara['LTR_Stream'])
from collections import defaultdict


from ltrDistManager import ltrDistManager

class branchLtrDistManager(ltrDistManager):
    def __init__(self, ltrParaFile, outDir):
        super().__init__(ltrParaFile=ltrParaFile, outDir=outDir)
        # Dist detection is based on a DFS search.
        # So class variations for recoding detection results are needed.
        self.detectRelDictList = []
        self.divideIdList = []

    def halveOriIdList(self, originOriIdList, originModIdSeries):
        toHalveData = self._3dPos.loc[originModIdSeries,:]
        kmeans = KMeans(n_clusters=2, init='k-means++', random_state=42).fit(toHalveData)
        halveRel = pd.Series(kmeans.labels_, index=toHalveData.index)
        oriIdSubList1 = []
        oriIdSubList2 = []
        for oriId in originOriIdList:
            modId = self.oriId2modId[oriId]
            if halveRel[modId]==0:
                oriIdSubList1.append(oriId)
            elif halveRel[modId]==1:
                oriIdSubList2.append(oriId)
            else:
                raise RuntimeError(f'Wrong index in halveRel!')
        return oriIdSubList1,oriIdSubList2

    def detectASet(self, toDetectOriIdList, dfsDepth, divideIdSuffix):
        # dfsDepth should not be bigger than 2.
        # (Means at most four parts for a final branch could be divided.)
        if dfsDepth>2:
            return None
        nowSetDetectRel = self.distDetector.detectOriIdList(toDetectOriIdList)
        toDetectModIdSeries = pd.Series([self.oriId2modId[oriId] for oriId in toDetectOriIdList])
        toDetectModIdSeries.drop_duplicates(inplace=True)
        oriIdSubList1, oriIdSubList2 = self.halveOriIdList(originOriIdList=toDetectOriIdList,
                                                           originModIdSeries=toDetectModIdSeries)
        subSet1DetectRel = self.distDetector.detectOriIdList(oriIdSubList1)
        subSet2DetectRel = self.distDetector.detectOriIdList(oriIdSubList2)

        if dfsDepth<2:
            if nowSetDetectRel['validOriIdList'] is None:
                if not (subSet1DetectRel['validOriIdList'] is None):
                    self.detectASet(oriIdSubList1, dfsDepth+1, f'{divideIdSuffix}A')
                if not (subSet2DetectRel['validOriIdList'] is None):
                    self.detectASet(oriIdSubList2, dfsDepth+1, f'{divideIdSuffix}B')
            else:
                furtherDfsFlag = False
                if not (subSet1DetectRel['validOriIdList'] is None):
                    if subSet1DetectRel['enrichScore'] > 2*nowSetDetectRel['enrichScore']:
                        furtherDfsFlag = True
                        self.detectASet(oriIdSubList1, dfsDepth+1, f'{divideIdSuffix}A')
                if not (subSet2DetectRel['validOriIdList'] is None):
                    if subSet2DetectRel['enrichScore'] > 2*nowSetDetectRel['enrichScore']:
                        furtherDfsFlag = True
                        self.detectASet(oriIdSubList2, dfsDepth+1, f'{divideIdSuffix}B')
                if not furtherDfsFlag:
                    myBranchId = self.modId2branchId[int(toDetectModIdSeries.to_list()[0])]
                    divideId = f'{myBranchId}__{divideIdSuffix}'
                    self.divideIdList.append(divideId)
                    self.detectRelDictList.append(nowSetDetectRel)
        else:
            if not (nowSetDetectRel['validOriIdList'] is None):
                myBranchId = self.modId2branchId[int(toDetectModIdSeries.to_list()[0])]
                divideId = f'{myBranchId}__{divideIdSuffix}'
                self.divideIdList.append(divideId)
                self.detectRelDictList.append(nowSetDetectRel)

    def detectSubLtrDist(self):
        # First clear all the detection results. (Though it is not permitted logically.)
        self.divideIdList = []
        self.detectRelDictList = []

        branchId2oriIdList = self.getBranchId2oriIdList()
        for branchId in self.finalBranchIdSet:
            self.detectASet(toDetectOriIdList = branchId2oriIdList[branchId],
                            dfsDepth = 0,
                            divideIdSuffix = '')

        self.detectRelWritter(self.divideIdList, self.detectRelDictList, f'branchSearch')
        self.detectRelPlotter(self.divideIdList, self.detectRelDictList, f'branchSearch')

# myDistManager = branchLtrDistManager(baseD='/data/home/xutun/papaverEvolution/mei/pro2',
#                                      outDir='/data/home/xutun/papaverEvolution/mei/pro2/testLtrDistDetector/branchSearch')
# myDistManager.detectSubLtrDist()
