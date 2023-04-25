import sys
sys.path.append('/data/home/xutun/mySrc/mei')

from ltrDistManager import ltrDistManager


class kmeansBranchLtrDistManager(ltrDistManager):
    # The detectSubLtrDist rule in this class was performed at the branch level.
    # Each branch would be divided into k parts, and only the part with enrichScore
    # higher than the enrichScore of the total branch would be selected as validOriId.

    def __init__(self, ltrParaFile, outDir, k, rootNode=None):
        super().__init__(ltrParaFile=ltrParaFile, outDir=outDir, rootNode=rootNode)
        self.branchId2oriIdList = self.getBranchId2oriIdList()
        self.k = k

    def detectSubLtrDist(self, foldChange=1, onlyDetectFinalBranch=True):
        # For each branch, divide the branch into k components, and to check whether any of the
        # component has enrich score more foldChange (usually 1 in the modified strategy) than
        # the total branch.
        k = self.k
        divideIdList = []
        detectRelDictList = []
        branchSet = {}
        if onlyDetectFinalBranch:
            branchSet = self.finalBranchIdSet
        else:
            branchSet = self.branchId2oriIdList
        for finalBranchId in branchSet:
            oriIdList = self.branchId2oriIdList[finalBranchId]
            totDetectRelDict = self.distDetector.detectOriIdList(oriIdList)
            dividedOriIdLists = self.divideOriIdListByKmeans(oriIdList,k)
            dividedDetectRelDictList = [self.distDetector.detectOriIdList(dividedOriIdList) \
                                        for dividedOriIdList in dividedOriIdLists]
            validDividedDetectRelDictList = []
            kthIdxList = []
            for kthIdx,dividedDetectRelDict in enumerate(dividedDetectRelDictList):
                if not (dividedDetectRelDict['validOriIdList'] is None):
                    if totDetectRelDict['validOriIdList'] is None or \
                       dividedDetectRelDict['enrichScore'] >= foldChange*totDetectRelDict['enrichScore']:
                        validDividedDetectRelDictList.append(dividedDetectRelDict)
                        kthIdxList.append(kthIdx)

            if 0 < len(validDividedDetectRelDictList) < k:
                mergedOriIdList = []
                for kthId in kthIdxList:
                    mergedOriIdList += dividedOriIdLists[kthId]
                mergedDetectRelDict = self.distDetector.detectOriIdList(mergedOriIdList)
                if not (mergedDetectRelDict['validOriIdList'] is None):
                    detectRelDictList.append(mergedDetectRelDict)
                    divideIdList.append(f'{finalBranchId}')
            elif not (totDetectRelDict['validOriIdList'] is None):
                detectRelDictList.append(totDetectRelDict)
                divideIdList.append(f'{finalBranchId}')

        return divideIdList, detectRelDictList, f'kmeans_{k}'

    def plotResults(self, divideIdList, detectRelDictList, outPre):
        self.detectRelWritter(divideIdList, detectRelDictList, outPre)
        self.detectRelPlotter(divideIdList, detectRelDictList, outPre)


# myDistManager = kmeansBranchLtrDistManager(baseD='/data/home/xutun/papaverEvolution/mei/pro2',
#                                            outDir='/data/home/xutun/papaverEvolution/mei/pro2/testLtrDistDetector/kmeansBranch')
# myDistManager.detectSubLtrDist(k=2, foldChange=2)
