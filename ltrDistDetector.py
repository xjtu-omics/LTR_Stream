from collections import defaultdict
import pandas as pd
import numpy as np
from scipy import stats
from bisect import bisect_left
import sys
import yaml

with open('envConfig.yaml', 'r') as of:
    ttPara = yaml.safe_load(of)
    sys.path.append(ttPara['LTR_Stream'])

from fileSysSupplier import fileSysSupplier
import ttUtils
from ttlib import bainfo
import ttDb


class ltrDistDetector:

    def __init__(self, ltrParaFile):
        self.sysSupplier = fileSysSupplier(ltrParaFile)
        self.init_ltrBaseInfo()
        self.species2chrNum = self.init_getSpecies2chrNum()
        self.species2genomeInfo = self.init_getSpecies2genomeInfo()
        self.init_genomeInfo()

    def init_ltrBaseInfo(self):
        self.seqName2bainfo = ttUtils.seqName2bainfo

    def init_getSpecies2genomeInfo(self):
        """
        # An example for infoDict
        info = {
            'ymr': '/data/home/xutun/papaverEvolution/mei/pro2/ref/ymr.ref.fa.fai',
            'hn1': '/data/home/xutun/papaverEvolution/mei/pro2/ref/hn1.ref.fa.fai',
            'dcw': '/data/home/xutun/papaverEvolution/mei/pro2/ref/dcw.ref.fa.fai'
        }
        """
        info = self.sysSupplier.getSpeId2faiFile()
        return info

    def init_getSpecies2chrNum(self):
        """
        # An example for species2chrNum (Actually speId2chrNum)
        species2chrNum = {
            'ymr': 7,
            'hn1': 11,
            'dcw': 22
        }
        """
        species2chrNum = ttDb.getSpeId2chrNum(self.sysSupplier.refParaFile)
        return species2chrNum

    def init_genomeInfo(self):
        species2genomeInfo = self.species2genomeInfo
        species2chrNum = self.species2chrNum
        totGenomeLen = 0
        chrmList = []
        self.chrm2len = defaultdict(int)
        self.sortedChrmList = []
        for spe in species2genomeInfo:
            genomeInfoData = pd.read_table(species2genomeInfo[spe], sep='\t', header=None)
            i = 0
            for ind, row in genomeInfoData.iterrows():
                i += 1
                if i > species2chrNum[spe]:
                    break
                totGenomeLen += row[1]
                chrmList.append(row[0])
                self.chrm2len[row[0]] = row[1]
                self.sortedChrmList.append(row[0])
        self.genomeLen = totGenomeLen
        self.firstChrm = chrmList[0]
        self.lastChrm = chrmList[-1]

    def getLinePos2offNum(self, edgeList, removePrefixNum):
        # Step I: get all the edges to be removed
        # Step II: sort them by the ed position, because only the
        #     the ed positions of removed edges matter
        relDict = {0: 0}
        edgeLenSortedIdx = sorted(range(len(edgeList)),
                                  key=lambda idx: (edgeList[idx]['edPos'] - edgeList[idx]['stPos']),
                                  reverse=True)
        toRemoveEdgeList = [edgeList[edgeLenSortedIdx[i]] for i in range(removePrefixNum)]
        edSortedToRemoveEdgeList = sorted(toRemoveEdgeList,
                                          key=lambda edg: edg['edPos'])
        preOffSet = 0
        for edge in edSortedToRemoveEdgeList:
            preOffSet += (edge['edPos'] - edge['stPos'])
            relDict[edge['edPos']] = preOffSet
        return relDict

    def getUnifiedPosDis(self, edgeList, removePrefixNum, linePos2offNum):
        # Remove the segments that to be removed, and recalculate the linePos (defined as unifiedPosDis)
        edgeLenSortedIdx = sorted(range(len(edgeList)),
                                  key=lambda idx: (edgeList[idx]['edPos'] - edgeList[idx]['stPos']),
                                  reverse=True)
        toRemoveEdgeList = [edgeList[edgeLenSortedIdx[i]] for i in range(removePrefixNum)]
        linePos2removeNum = defaultdict(int)
        for edge in toRemoveEdgeList:
            linePos2removeNum[edge['stPos']] += 1
            linePos2removeNum[edge['edPos']] += 1

        unUnifiedPosDis = []
        for edge in edgeList:
            if edge['stOriIdx'] == -1:
                if linePos2removeNum[edge['stPos']] < 1:
                    unUnifiedPosDis.append(edge['stPos'])
            else:
                if linePos2removeNum[edge['stPos']] < 2:
                    unUnifiedPosDis.append(edge['stPos'])
            if edge['edOriIdx'] == -1:
                if linePos2removeNum[edge['edPos']] < 1:
                    unUnifiedPosDis.append(edge['edPos'])

        linePos = sorted(list(linePos2offNum.keys()))
        unifiedPosDis = []
        for pos in unUnifiedPosDis:
            offNum = linePos2offNum[pos] if pos in linePos2offNum \
                else linePos2offNum[linePos[bisect_left(linePos, pos) - 1]]
            unifiedPosDis.append(pos - offNum)
        return unifiedPosDis

    def getVirtualDis(self, linePos2offNum, genomeLen, toSimulateNum):
        # Return a random distribution on the linePos2offNum range
        maxOffNum = max(list(linePos2offNum.values()))
        simulateRange = genomeLen - maxOffNum
        relDis = list(np.random.randint(0, simulateRange, toSimulateNum))
        return relDis

    def getRandomDisList(self, linePos2offNum, genomeLen, toSimulateNumInEachDis, toSimulateDisNum):
        # Return a list of random distribution
        relList = []
        for i in range(toSimulateDisNum):
            relList.append(self.getVirtualDis(linePos2offNum, genomeLen, toSimulateNumInEachDis))
        return relList

    def getPercentileOfNowSet(self, toTestDis, bakDis, randomDisList):
        # Get the percentile of the distance of testDist to bakDist among randomDist to bakDist
        # This usually seriously affected by a bad bakDist.
        testWasDistance = stats.wasserstein_distance(toTestDis, bakDis)
        randomWasDistance = sorted([stats.wasserstein_distance(randomDis, bakDis)
                                    for randomDis in randomDisList])
        relPercentile = float(bisect_left(randomWasDistance, testWasDistance)) / len(randomWasDistance)
        return relPercentile

    def getPercentileOfNowSet2(self, toTestDis, randomDisList1, randomDisList2):
        # Warning: not modified for suiting class here, only adding 'self' in parameter list
        testWasDistanceList = [stats.wasserstein_distance(toTestDis, bakDis)
                               for bakDis in randomDisList2]
        randomWasDistanceList = [stats.wasserstein_distance(randomDis, bakDis)
                                 for randomDis, bakDis in zip(randomDisList1, randomDisList2)]
        rel = 0
        for i, j in zip(testWasDistanceList, randomWasDistanceList):
            if i > j:
                rel += 1
        return rel

    def getPercentageOfKsTest(self, toTestDis, genomeLen, linePos2offNum):
        maxOffNum = max(list(linePos2offNum.values()))
        simulateRange = genomeLen - maxOffNum
        pv = stats.kstest(toTestDis, stats.randint.cdf, args=(0, simulateRange)).pvalue
        return pv

    def getRemovePercentage(self, pvalueList):
        # Return the percentage of removed top lens so that the rest
        # of lens is similar with a [稳定的] random dist.
        pBordary = 0.05
        consecutiveNum = 5
        findOutlier = False
        for i in range(consecutiveNum):
            if pvalueList[i] < pBordary:
                findOutlier = True
                break
        if not findOutlier:
            return 0
        for i in range(consecutiveNum, len(pvalueList)):
            findOutlier = False
            for j in range(i - consecutiveNum + 1, i):
                if pvalueList[j] <= pBordary:
                    findOutlier = True
            if not findOutlier:
                return i
        return 100

    def getRemovedOriIdxDict(self, edgeList, toRemovePercentage):
        # OriIdx is not oriId.
        # It's the original index of linePos
        if toRemovePercentage == 0 or toRemovePercentage == 100:
            return None
        removePrefixNum = int(toRemovePercentage / 100 * len(edgeList))
        edgeLenSortedIdx = sorted(range(len(edgeList)),
                                  key=lambda idx: (edgeList[idx]['edPos'] - edgeList[idx]['stPos']),
                                  reverse=True)

        toRemoveEdgeList = [edgeList[edgeLenSortedIdx[i]] for i in range(removePrefixNum)]
        oriIdx2removeNum = defaultdict(int)
        for edge in toRemoveEdgeList:
            oriIdx2removeNum[edge['stOriIdx']] += 1
            oriIdx2removeNum[edge['edOriIdx']] += 1

        removedOriIdxDict = defaultdict(int)
        for oriIdx in oriIdx2removeNum:
            removeNum = oriIdx2removeNum[oriIdx]
            if removeNum >= 2:
                removedOriIdxDict[oriIdx] = 1
        return removedOriIdxDict

    def getRemovedEdgeIdxSet(self, edgeList, toRemovePercentage):
        if toRemovePercentage == 0 or toRemovePercentage == 100:
            return None
        removePrefixNum = int(toRemovePercentage / 100 * len(edgeList))
        edgeLenSortedIdx = sorted(range(len(edgeList)),
                                  key=lambda idx: (edgeList[idx]['edPos'] - edgeList[idx]['stPos']),
                                  reverse=True)
        removedEdgeIdxSet = defaultdict(int)
        for i in range(removePrefixNum):
            removedEdgeIdxSet[edgeLenSortedIdx[i]] = 1
        return removedEdgeIdxSet

    def getEnrichScore(self, edgeList, toRemovePercentage, removedOriIdxDict):
        genomeLen = self.genomeLen
        if toRemovePercentage == 0 or toRemovePercentage == 100:
            return 1
        totLtrNum = len(edgeList) - 1
        removedLtrNum = len(removedOriIdxDict)
        ltrPerc = 1 - (removedLtrNum / totLtrNum)
        removePrefixNum = int(toRemovePercentage / 100 * len(edgeList))
        edgeLenSortedIdx = sorted(range(len(edgeList)),
                                  key=lambda idx: (edgeList[idx]['edPos'] - edgeList[idx]['stPos']),
                                  reverse=True)
        toRemoveEdgeList = [edgeList[edgeLenSortedIdx[i]] for i in range(removePrefixNum)]
        toRemoveEdgeLen = 0
        for edg in toRemoveEdgeList:
            toRemoveEdgeLen += (edg['edPos'] - edg['stPos'] + 1)
        rangePerc = 1 - (toRemoveEdgeLen / genomeLen)
        return ltrPerc / rangePerc

    def transOriId2linePos(self, oriIdList):
        filteredOriIdList = []
        linePosList = []
        chr2prefix = {}
        species2genomeInfo = self.species2genomeInfo
        species2chrNum = self.species2chrNum
        tmpPrefix = 0
        for spe in species2genomeInfo:
            genomeInfoData = pd.read_table(species2genomeInfo[spe], sep='\t', header=None)
            i = 0
            for ind, row in genomeInfoData.iterrows():
                i += 1
                if i > species2chrNum[spe]:
                    break
                chr2prefix[row[0]] = tmpPrefix
                tmpPrefix += row[1]
        for oriId in oriIdList:
            myBain = self.seqName2bainfo(oriId)
            if not (myBain.chr in chr2prefix):
                continue
            filteredOriIdList.append(oriId)
            linePosList.append(chr2prefix[myBain.chr] + myBain.st)
        return filteredOriIdList, linePosList

    def transLinePos2edge(self, toDetectLinePosList, linePosSortedIdx):
        edgeList = []
        edgeI = 0
        # Add the 0 to the first LTR-RT
        edgeList.append({
            'edgeId': edgeI,
            'stPos': 0,
            'edPos': toDetectLinePosList[linePosSortedIdx[0]],
            'stOriIdx': -1,
            'edOriIdx': linePosSortedIdx[0]
        })
        edgeI += 1
        for i in range(len(linePosSortedIdx) - 1):
            tmpEdge = {
                'edgeId': edgeI,
                'stPos': toDetectLinePosList[linePosSortedIdx[i]],
                'edPos': toDetectLinePosList[linePosSortedIdx[i + 1]],
                'stOriIdx': linePosSortedIdx[i],
                'edOriIdx': linePosSortedIdx[i + 1]
            }
            edgeI += 1
            edgeList.append(tmpEdge)

        # Add the last LTR-RT to the end of the genome
        edgeList.append({
            'edgeId': edgeI,
            'stPos': toDetectLinePosList[linePosSortedIdx[-1]],
            'edPos': self.genomeLen,
            'stOriIdx': -1,
            'edOriIdx': linePosSortedIdx[-1]
        })
        edgeI += 1
        return edgeList

    def getEnrichRangeBainfoList(self, edgeList, removedEdgeIdxSet, oriIdList):
        # Get genomic ranges of enriched regions.
        def addGenomicRangeBainfo(relList, rangeStOriIdx, rangeEdOriIdx):
            rangeStOriId = f'{self.firstChrm}:0-0(+)' if rangeStOriIdx < 0 \
                      else oriIdList[rangeStOriIdx]
            lastChrmLen = self.chrm2len[self.lastChrm]
            rangeEdOriId = f'{self.lastChrm}:{lastChrmLen}-{lastChrmLen}(+)' if rangeEdOriIdx < 0 \
                      else oriIdList[rangeEdOriIdx]

            rangeStBainfo, rangeEdBainfo = self.seqName2bainfo(rangeStOriId), \
                                           self.seqName2bainfo(rangeEdOriId)
            if rangeStBainfo.chr != rangeEdBainfo.chr:
                relList.append(bainfo(_chr=rangeStBainfo.chr,
                                      _st=rangeStBainfo.st,
                                      _ed=self.chrm2len[rangeStBainfo.chr]))
                relList.append(bainfo(_chr=rangeEdBainfo.chr,
                                      _st=0,
                                      _ed=rangeEdBainfo.st))
            else:
                relList.append(bainfo(_chr=rangeStBainfo.chr,
                                      _st=rangeStBainfo.st,
                                      _ed=rangeEdBainfo.st))

        enrichRangeList = []
        if removedEdgeIdxSet is None:
            return None
        rangeStOriIdx = None
        rangeEdOriIdx = None
        for i in range(len(edgeList)):
            if i in removedEdgeIdxSet:
                if not (rangeStOriIdx is None):
                    addGenomicRangeBainfo(relList=enrichRangeList,
                                          rangeStOriIdx=rangeStOriIdx,
                                          rangeEdOriIdx=rangeEdOriIdx)
                rangeStOriIdx = None
                rangeEdOriIdx = None
            else:
                edge = edgeList[i]
                if rangeStOriIdx is None:
                    rangeStOriIdx = edge['stOriIdx']
                rangeEdOriIdx = edge['edOriIdx']
        if not (rangeStOriIdx is None):
            addGenomicRangeBainfo(relList=enrichRangeList,
                                  rangeStOriIdx=rangeStOriIdx,
                                  rangeEdOriIdx=rangeEdOriIdx)
        return enrichRangeList

    def getValidOriIdList(self, removedOriIdxDict, oriIdList, linePosSortedIdx):
        if removedOriIdxDict is None:
            return None
        validOriIdList = []
        for i in range(len(oriIdList)):
            ind = linePosSortedIdx[i]
            if not (ind in removedOriIdxDict):
                validOriIdList.append(oriIdList[ind])
        return validOriIdList

    def detectOriIdList(self, oriIdList):
        genomeLen = self.genomeLen
        # oriIdList would be filtered to preclude LTR-RTs on scaffolds.
        # oriIdList and toDetectLinePosList could be indexed by linePosSortedIdx
        if len(oriIdList) < 10:
            # To prevent all (at the most 80%) the LTR-RTs were removed.
            return {'percentile': None,
                    'toRemovePercentage': None,
                    'enrichScore': None,
                    'enrichRangeBainfoList': None,
                    'validOriIdList': None}
        oriIdList, toDetectLinePosList = self.transOriId2linePos(oriIdList)
        linePosSortedIdx = sorted(range(len(toDetectLinePosList)),
                                  key=lambda idx: toDetectLinePosList[idx])
        edgeList = self.transLinePos2edge(toDetectLinePosList=toDetectLinePosList,
                                          linePosSortedIdx=linePosSortedIdx)
        removePrefixNumList = [0]
        relPercentileList = []
        for i in range(80):
            removePrefixNumList.append(int(i / 100 * len(edgeList)))

        for removePrefixNum in removePrefixNumList:
            linePos2offNum = self.getLinePos2offNum(edgeList, removePrefixNum)
            unifiedPosDis = self.getUnifiedPosDis(edgeList, removePrefixNum, linePos2offNum)

            testPercentile = self.getPercentageOfKsTest(toTestDis=unifiedPosDis,
                                                        genomeLen=genomeLen,
                                                        linePos2offNum=linePos2offNum)
            relPercentileList.append(testPercentile)

        toRemovePercentage = self.getRemovePercentage(pvalueList=relPercentileList)
        # removedOriIdxDict is also the flag of whether the distibution this subtype
        # is partially enriched in any genomic regions.
        # If removedOriIdxDict is None, this usually means that the subtype is not
        # partially distributed. (Especially for some subtypes with little LTR-RTs)
        removedOriIdxDict = self.getRemovedOriIdxDict(edgeList=edgeList,
                                                      toRemovePercentage=toRemovePercentage)
        removedEdgeIdxSet = self.getRemovedEdgeIdxSet(edgeList=edgeList,
                                                      toRemovePercentage=toRemovePercentage)
        enrichRangeBainfoList = self.getEnrichRangeBainfoList(edgeList=edgeList,
                                                              removedEdgeIdxSet=removedEdgeIdxSet,
                                                              oriIdList=oriIdList)
        validOriIdList = self.getValidOriIdList(removedOriIdxDict=removedOriIdxDict,
                                                oriIdList=oriIdList,
                                                linePosSortedIdx=linePosSortedIdx)
        enrichScore = self.getEnrichScore(edgeList=edgeList,
                                          toRemovePercentage=toRemovePercentage,
                                          removedOriIdxDict=removedOriIdxDict)

        return {'percentile': relPercentileList,
                'toRemovePercentage': toRemovePercentage,
                'enrichScore': enrichScore,
                'enrichRangeBainfoList': enrichRangeBainfoList,
                'validOriIdList': validOriIdList}
