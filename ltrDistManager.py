from collections import defaultdict
import sys
from sklearn.cluster import KMeans
import pandas as pd
import yaml

with open('envConfig.yaml', 'r') as of:
    ttPara = yaml.safe_load(of)
    sys.path.append(ttPara['LTR_Stream'])

from ltrDistDetector import ltrDistDetector
from ttStream import loadWorkData
from meiPlottter import meiPlotter
import ttUtils
from fileSysSupplier import fileSysSupplier


class ltrDistManager:
    def __init__(self, ltrParaFile, outDir, rootNode=None):
        self.rootNode = rootNode
        self.sysSupplier = fileSysSupplier(ltrParaFile)
        self.init_workDir(outDir)
        self.distDetector = ltrDistDetector(ltrParaFile)
        self.modId2branchId = ttUtils.getModId2branchId(f'{self.danteD}/tot.pgl.hdf5')
        self.modId2oriIdSet = ttUtils.getModId2OriIdSet(f'{self.danteD}/finalModSeq.tab',
                                                        f'{self.danteD}/toCalRest.modSeq2modId.tab')
        self.oriId2modId = ttUtils.getOriId2modId(f'{self.danteD}/finalModSeq.tab',
                                                  f'{self.danteD}/toCalRest.modSeq2modId.tab')
        self._3dPos = ttUtils.getPcoaDisMat(f'{self.danteD}/toCalRest.pcoaDis.csv')
        self.finalBranchIdSet = self.init_getFinalBranchIdSet()
        self.plotter = meiPlotter(self, initTyp='ltrDistManager', ltrParaFile=ltrParaFile)
        self.maxPlotTypNum = self.plotter.maxTypNum
        self.adata = loadWorkData(f'{self.danteD}/tot.pgl.hdf5')

    def init_workDir(self, outDir):
        self.baseD = self.sysSupplier.baseD
        self.outDir = outDir
        self.danteD = f'{self.baseD}/dante'

    def init_getFinalBranchIdSet(self):
        # Select branchId related to leaf nodes.
        branchIdSet = defaultdict(int)
        for modId in self.modId2branchId:
            branchId = self.modId2branchId[modId]
            branchIdSet[branchId] = 1
        nodeId2edgeNum = defaultdict(int)
        for branchId in branchIdSet:
            nodeIdList = branchId.split('__')
            for nodeId in nodeIdList:
                nodeId2edgeNum[nodeId] += 1
        finalBranchIdSet = defaultdict(int)
        for branchId in branchIdSet:
            nodeIdList = branchId.split('__')
            for nodeId in nodeIdList:
                if (not (self.rootNode is None)) and nodeId == self.rootNode:
                    continue
                if nodeId2edgeNum[nodeId] == 1:
                    finalBranchIdSet[branchId] = 1
        return finalBranchIdSet

    def getBranchId2oriIdList(self):
        finalBranchId2oriIdSet = defaultdict(list)
        for oriId in self.oriId2modId:
            modId = self.oriId2modId[oriId]
            branchId = self.modId2branchId[modId]
            finalBranchId2oriIdSet[branchId].append(oriId)
        return finalBranchId2oriIdSet

    def divideOriIdListByKmeans(self, oriIdList, k):
        originModIdList = [self.oriId2modId[oriId] for oriId in oriIdList]
        originModIdSeries = pd.Series(originModIdList)
        originModIdSeries.drop_duplicates(inplace=True)
        toDivideData = self._3dPos.loc[originModIdSeries, :]
        kmeans = KMeans(n_clusters=k, init='k-means++', random_state=42).fit(toDivideData)
        divideRel = pd.Series(kmeans.labels_, index=toDivideData.index)
        oriIdSubLists = []
        for i in range(k):
            oriIdSubLists.append([])
        for oriId in oriIdList:
            modId = self.oriId2modId[oriId]
            oriIdSubLists[divideRel[modId]].append(oriId)
        return oriIdSubLists

    def getDivideId2typ(self, divideIdList):
        uniqDivideIdList = list(pd.Series(divideIdList).unique())
        divideId2typ = defaultdict(int)
        for i in range(len(uniqDivideIdList)):
            divideId2typ[uniqDivideIdList[i]] = i + 1
        return divideId2typ

    def plotCircosDetectRel(self, outPrefix, divideIdList, validOriIdList, divideId2typ):
        # divideIdList: list of a str
        # validOriIdList: oriIdList should be same len with divideIdList
        outPdf = f'{self.outDir}/{outPrefix}.circos.pdf'
        typList = [divideId2typ[divideId] for divideId in divideIdList]
        typ2divideId = {divideId2typ[divideId]: divideId for divideId in divideId2typ}
        self.plotter.circosPlot(typList=typList,
                                oriIdList=validOriIdList,
                                typId2annot=typ2divideId,
                                outFigFile=outPdf)

    def plot3DSpaceDetectRel(self, outPrefix, divideIdList, validOriIdList, divideId2typ):
        # divideIdList: list of a str
        # validOriIdList: oriIdList should be same len with divideIdList
        outGif = f'{self.outDir}/{outPrefix}.3D.gif'
        modId2plotTyp = defaultdict(int)
        for divideId, oriId in zip(divideIdList, validOriIdList):
            modId2plotTyp[self.oriId2modId[oriId]] = divideId2typ[divideId]
        toPlotTypList = []
        for modId in self.modId2branchId:
            if not (modId in modId2plotTyp):
                # 0 for grey background scatters.
                toPlotTypList.append(0)
            else:
                # Setting the same color with that in circosPlot.
                toPlotTypList.append(modId2plotTyp[modId])
        typ2divideId = {divideId2typ[divideId]: divideId for divideId in divideId2typ}
        self.plotter.p3DPlot(adata=self.adata,
                             typList=toPlotTypList,
                             typId2annot=typ2divideId,
                             outGif=outGif)

    def plotBatchRel(self, tmpOutPrefix, toPlotDivideIdList, toPlotValidOriIdList, divideId2typ):
        self.plotCircosDetectRel(tmpOutPrefix, toPlotDivideIdList, toPlotValidOriIdList, divideId2typ)
        self.plot3DSpaceDetectRel(tmpOutPrefix, toPlotDivideIdList, toPlotValidOriIdList, divideId2typ)

    def detectRelPlotter(self, divideIdList, detectRelDictList, outPrefix, divideIdPlotOrderDict=None):
        plotBatchId = 0
        tmpInd = 0
        toPlotDivideIdList = []
        toPlotValidOriIdList = []
        for divId, relDict in zip(divideIdList, detectRelDictList):
            if (len(toPlotDivideIdList) > 0) and (tmpInd % self.maxPlotTypNum == 0):
                tmpOutPrefix = f'{outPrefix}_{plotBatchId}'
                divideId2typ = divideIdPlotOrderDict
                if divideId2typ is None:
                    divideId2typ = self.getDivideId2typ(toPlotDivideIdList)
                self.plotBatchRel(tmpOutPrefix, toPlotDivideIdList, toPlotValidOriIdList, divideId2typ)
                plotBatchId += 1
                toPlotDivideIdList = []
                toPlotValidOriIdList = []
            toPlotDivideIdList += [divId for i in range(len(relDict['validOriIdList']))]
            toPlotValidOriIdList += relDict['validOriIdList']
            tmpInd += 1
        if len(toPlotDivideIdList) > 0:
            tmpOutPrefix = f'{outPrefix}_{plotBatchId}'
            divideId2typ = divideIdPlotOrderDict
            if divideId2typ is None:
                divideId2typ = self.getDivideId2typ(toPlotDivideIdList)
            self.plotBatchRel(tmpOutPrefix, toPlotDivideIdList, toPlotValidOriIdList, divideId2typ)

    def detectRelWritter(self, divideIdList, detectRelDictList, outPrefix):
        # A summary tabFile. # divideId, enrichScore, ltrNumber
        # A oriId tabFile. # oriId, divideId
        # A bedFile. # Range of enrich regions of each divideId
        summaryTabFile = f'{self.outDir}/{outPrefix}_summary.tab'
        oriIdTabFile = f'{self.outDir}/{outPrefix}_oriId.tab'
        bedFile = f'{self.outDir}/{outPrefix}_enrichRegion.bed'
        with open(summaryTabFile, 'w') as of:
            print('# Cluster', 'EnrichScore', 'No. of LTR-RT', sep='\t', file=of)
            for divId, relDict in zip(divideIdList, detectRelDictList):
                print(divId, relDict['enrichScore'], len(relDict['validOriIdList']), sep='\t', file=of)
        with open(oriIdTabFile, 'w') as of:
            print('# LTR-RT', 'Cluster', sep='\t', file=of)
            for divId, relDict in zip(divideIdList, detectRelDictList):
                oriIdList = relDict['validOriIdList']
                for oriId in oriIdList:
                    print(oriId, divId, sep='\t', file=of)
        with open(bedFile, 'w') as of:
            for divId, relDict in zip(divideIdList, detectRelDictList):
                bainList = relDict['enrichRangeBainfoList']
                for bain in bainList:
                    print(bain.chr, bain.st, bain.ed, divId, '.', '+', sep='\t', file=of)

    def detectSubLtrDist(self):
        pass
