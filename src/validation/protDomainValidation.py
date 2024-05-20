import pandas as pd
import numpy as np
import re
import os
import pysam

import ttUtils
from infoIntegrator import infoIntegrator as classInfoIntegrator
from ttlib import constructEnv
from validation.p2pValidation import p2pValidator

def fq2df(fastq, sortNames=True):
    names = []
    seqs  =[]
    with pysam.FastxFile(fastq) as inf:
        for rec in inf:
            names.append(rec.name)
            seqs.append(rec.sequence)
    relDf = pd.DataFrame(dict(name=names, seq=seqs))
    if sortNames:
        relDf = relDf.sort_values(by='name')
    return relDf
def df2fq(infoDf, outFa):
    with open(outFa, 'w') as of:
        for ind, row in infoDf.iterrows():
            print(f'>{row["name"]}', file=of)
            print(row.seq, file=of)
class danteRelIntegrator:
    infoDfColumns = ['oriId', 'st', 'ed', 'len', 'strand', 'domName', 'domSeq']
    def __init__(self, infoDf):
        self.infoDf = infoDf
    @classmethod
    def init_readDanteRelGff(cls, danteRelGff):
        gffDf = pd.read_table(danteRelGff, header=None, comment='#', sep='\t')
        gffDf.columns = ['oriId', 'dante', 'proDom', 'st', 'ed', 'len', 'strand', 'dot', 'annot']
        domNameList = []
        domSeqList = []
        for ind, row in gffDf.iterrows():
            domLinName = re.findall('Final_Classification=(.+?);', row['annot'])[0] + '|' + \
                         re.findall('^Name=(.+?);', row['annot'])[0]
            domNameList.append(domLinName.split('|')[-1])

            domSeq = re.findall('Query_Seq=(.+?);', row['annot'])
            if len(domSeq)>0:
                domSeq = re.sub(r'[^a-zA-Z]', '', domSeq[0])
            else:
                domSeq = ''
            domSeqList.append(domSeq)
        infoDf = pd.concat([gffDf[['oriId', 'st', 'ed', 'len', 'strand']],
                            pd.DataFrame(dict(domName=domNameList, domSeq=domSeqList), index=gffDf.index)],
                           axis=1)
        infoDf = infoDf[infoDf['domSeq']!='']
        return infoDf
    @classmethod
    def init_fromDanteGff(cls, danteRelGff):
        return danteRelIntegrator(cls.init_readDanteRelGff(danteRelGff))
    @classmethod
    def init_fromMultiDanteGff(cls, danteRelGffList):
        infoDf = pd.DataFrame(columns=cls.infoDfColumns)
        for danteRelGff in danteRelGffList:
            infoDf = pd.concat([infoDf, cls.init_readDanteRelGff(danteRelGff)],
                               axis=0)
        infoDf = infoDf.reset_index(drop=True)
        return danteRelIntegrator(infoDf)
    def oriId2inRangeId(self, ltrRetrievListFiles):
        def getCaliValue(ser, oriId2inRangeId):
            oriId = ser['oriId']
            inRangeId = oriId2inRangeId[oriId]
            oriBain = ttUtils.seqName2bainfo(oriId)
            inBain = ttUtils.seqName2bainfo(inRangeId)
            return inBain.st-oriBain.st
        oriId2inRangeId = {}
        for ltrRetrievListF in ltrRetrievListFiles:
            oriId2inRangeId = ttUtils.getOriId2inRangeId(ltrRetrievListF, relDict=oriId2inRangeId)
        caliValueSer = self.infoDf.apply(getCaliValue, oriId2inRangeId=oriId2inRangeId, axis=1)
        self.infoDf['st'] = self.infoDf['st']-caliValueSer
        self.infoDf['ed'] = self.infoDf['ed']-caliValueSer
        inRangeIdSer = pd.Series([oriId2inRangeId[oriId] for oriId in list(self.infoDf['oriId'])], index=self.infoDf.index)
        self.infoDf['oriId'] = inRangeIdSer
    def getDomainPercent(self, classInfoDf):
        subDf = self.infoDf[self.infoDf['oriId'].isin(classInfoDf['oriId'])]
        subDf = subDf.drop_duplicates(subset=['oriId', 'domName'])
        countSer = subDf['domName'].value_counts()
        return countSer/classInfoDf.shape[0]
    def selectDomSeq(self, classInfoDf, domName, retDf=False):
        subDf = self.infoDf[self.infoDf['oriId'].isin(classInfoDf['oriId'])]
        subDf = subDf[ subDf['domName']==domName ]
        subDf = subDf.sort_values(by='oriId')
        if retDf:
            return subDf
        return subDf['domSeq']
    def getValidOriIdSer(self, classInfoDf, domNameSer):
        def checkHaveAllDom(df, domNameSer):
            df = df.drop_duplicates(subset=['domName'])
            return domNameSer.isin(df.domName).all()
        subDf = self.infoDf[self.infoDf['oriId'].isin(classInfoDf['oriId'])]
        indexSer = subDf.groupby('oriId').apply(checkHaveAllDom, domNameSer)
        indexDf = indexSer.reset_index()
        indexDf.columns = ['oriId', 'bol']
        return indexDf[indexDf.bol==True].oriId
    def dedupOriId(self):
        self.infoDf = self.infoDf.drop_duplicates(['oriId', 'domName'])
class domainValidator:
    def __init__(self, classTsvFile, domainInfo, workDir, removeCentPerc, filter=True):
        self.classInfo = classInfoIntegrator.init_fromTsv(classTsvFile)
        if filter:
            self.classInfo.filterCentPoint(perc=removeCentPerc)
        self.domainInfo = domainInfo
        self.randomState = 42
        self.init_dir(workDir)
    def init_dir(self, workDir):
        self.baseD = workDir
        self.protFaD = f'{self.baseD}/protFasta'
        self.treeD = f'{self.baseD}/tree'
        self.alnD = f'{self.baseD}/aln'
        constructEnv([self.baseD, self.protFaD, self.treeD, self.alnD])
    def checkNonAutoPerc(self):
        def isNonAutoTe(infoDf):
            if infoDf.shape[0]<6:
                return True
            myDoms = infoDf.domName.tolist()
            for dom in ['GAG', 'INT', 'PROT', 'RH', 'RT', 'aRH']:
                if dom not in myDoms:
                    return True
            return False
        myDomainInfoDf = self.domainInfo.infoDf[self.domainInfo.infoDf.oriId.isin(self.classInfo.infoDf.oriId)]
        flagIsNonAutoTe = myDomainInfoDf.groupby('oriId').apply(isNonAutoTe)
        return np.sum(flagIsNonAutoTe)/self.classInfo.infoDf.drop_duplicates('oriId').shape[0]
        # print(self.classInfo.infoDf.drop_duplicates('oriId').shape)
        # print(myDomainInfoDf.drop_duplicates('oriId').shape)


    def ranDomSeq(self, domName=None):
        ranDomSer = None
        if domName is None:
            ranDomSer = self.classInfo.infoDf.groupby('zoomIn').apply(self.selectMostCommonDomain)
        else:
            zoomInList = self.classInfo.getZoomInList()
            ranDomSer = pd.Series([domName for zoomIn in zoomInList], index = zoomInList)
        ranDomSeqDf = self.classInfo.infoDf.groupby('zoomIn', as_index=False).apply(self.ranZoomInDomSeq, ranDomSer=ranDomSer)
        # In case that this zoomIn level only contains one Class, there is no necessary to construct the phylo-tree.
        # This function has been discarsded since the multi-level cluster is proved to be incorrect.
        onlyOneClassZoomInSer = ranDomSeqDf.groupby('zoomIn').apply(lambda x:len(list(x.Class.drop_duplicates())))
        onlyOneClassZoomInSer = onlyOneClassZoomInSer[onlyOneClassZoomInSer==1]
        ranDomSeqDf = ranDomSeqDf[~(ranDomSeqDf.zoomIn.isin(onlyOneClassZoomInSer.index))]
        ranDomSeqDf.groupby('zoomIn').apply(self.writeRanDomSeqFa)
        return ranDomSeqDf
    def selectMostCommonDomain(self, classInfoDf):
        dmPercDf = classInfoDf.groupby('Class').apply(self.domainInfo.getDomainPercent)
        dmPercDf = dmPercDf.reset_index()
        dmPercDf = dmPercDf.rename(columns={'level_1':'domName', 'domName':'perc'})
        dmPercDf = dmPercDf.pivot(index='Class', columns='domName', values='perc')
        dmPercDf = dmPercDf.fillna(0)
        commPercSer = dmPercDf.apply(min, axis=0)
        return commPercSer.idxmax()
    def getDomSeqPerc(self, myClassDf, domSeqDf):
        tarClass = myClassDf.Class.iloc[0]
        return pd.DataFrame(dict(Class=[tarClass],
                                 perc=[np.sum(domSeqDf.Class==tarClass)/myClassDf.shape[0]]))
    def ranZoomInDomSeq(self, classInfoDf, ranDomSer, selectSize=10):
        # domName is zoomIn
        domName = ranDomSer[classInfoDf['zoomIn'].iloc[0]]
        # Get each Class in this zoomIn, ret is dataframe containing protein domain sequences
        domSeqSer = classInfoDf.groupby('Class').apply(self.domainInfo.selectDomSeq, domName=domName)
        domSeqDf = domSeqSer.reset_index(level='Class')
        # Filter those classes with rare protein domains
        classDomSeqPercDf = classInfoDf.groupby('Class', as_index=False).apply(self.getDomSeqPerc, domSeqDf=domSeqDf).reset_index(drop=True)
        classDomSeqPercDf = classDomSeqPercDf[classDomSeqPercDf.perc>0.2]
        domSeqDf = domSeqDf[domSeqDf.Class.isin(classDomSeqPercDf.Class)]
        # Random protein sequences.
        randomSeqDf = domSeqDf.groupby('Class').apply(lambda x:x.sample(n=min(selectSize, x.shape[0]),
                                                                        random_state=self.randomState, replace=False))
        randomSeqDf = randomSeqDf.reset_index(drop=True)
        randomSeqDf['zoomIn'] = pd.Series([classInfoDf['zoomIn'].iloc[0] for i in range(randomSeqDf.shape[0])],
                                          index=randomSeqDf.index)

        return randomSeqDf
    def writeRanDomSeqFa(self, seqInfoDf):
        zoomIn = seqInfoDf['zoomIn'].iloc[0]
        outFa = f'{self.protFaD}/{zoomIn}.randomProtSeq.fa'
        with open(outFa, 'w') as of:
            for i, (ind, row) in enumerate(seqInfoDf.iterrows()):
                print(f'>{row["Class"]}__{i}', file=of)
                print(row['domSeq'], file=of)
    @classmethod
    def runMuscle(cls, inputFa, outFile, outFormat='clw', calDis=True, alnFlag='-align'):
        softPath = '/data/home/xutun/miniconda3/envs/alntree/bin'
        flag = ''
        if outFormat=='clw':
            flag += '-clwstrict'
        os.system(f'''
            PATH={softPath}:$PATH
            muscle {alnFlag} {inputFa} {flag} -output {outFile} >/dev/null 2>&1
        ''')
        if calDis:
            os.system(f'''
                PATH={softPath}:$PATH
                clustalo -i {outFile} --percent-id --distmat-out={outFile}.dis --full --force >/dev/null 2>&1
            ''')
    def runIqTree(self, inputClw, outPrefix):
        softPath = '/data/home/xutun/miniconda3/envs/alntree/bin'
        softPath = '/data/home/testXT/miniconda3/envs/phyloTree/bin/'
        os.system(f'''
            PATH={softPath}:$PATH
            iqtree2 -s {inputClw} --prefix {outPrefix} --mem 2000G -T 48 --quiet -redo --keep-ident --fast >/dev/null 2>&1
        ''')
    def constructPhyloTree(self, ranDomSeqDf):
        zoomInList = list(ranDomSeqDf['zoomIn'].drop_duplicates())
        for zoomIn in zoomInList:
            print(f'Construct {zoomIn}...')
            inputFa = f'{self.protFaD}/{zoomIn}.randomProtSeq.fa'
            outClw = f'{self.alnD}/{zoomIn}.muscle.clw'
            outPrefix = f'{self.treeD}/{zoomIn}'
            self.runMuscle(inputFa, outClw)
            self.runIqTree(outClw, outPrefix)
    def plotPhyloTree(self, ranDomSeqDf, outPrefix=''):
        zoomInList = list(ranDomSeqDf['zoomIn'].drop_duplicates())
        for zoomIn in zoomInList:
            treeFile = f'{self.treeD}/{zoomIn}.treefile'
            softPath = '/data/home/testXT/miniconda3/envs/ggtree/bin/'
            outPdf = f'{self.treeD}/{outPrefix}{zoomIn}.pdf'
            os.system(f'''
                PATH={softPath}:$PATH
                Rscript /data/home/testXT/LTR_Stream/validation/test_ggtree.R {treeFile} {zoomIn} {outPdf}
            ''')
    def ppline_geneProtDomPhy(self, domList, onlyStaPerc=False):
        for dom in domList:
            outPdf = f'{self.treeD}/{dom}_*.pdf'
            os.system(f'''
                            rm -rf {outPdf}
                      ''')
            ranDomSeqDf = self.ranDomSeq(dom)
            if ranDomSeqDf.shape[0] == 0:
                print(f'No {dom} should be constructed!')
                continue
            if onlyStaPerc:
                continue
            self.constructPhyloTree(ranDomSeqDf)
            self.plotPhyloTree(ranDomSeqDf,f'{dom}_')
    def plotDomainPercHeat(self, zoomInLevel=None, ax=None, figSize=None, Class2x=None, ClassColors=None):
        infoDf = None
        if zoomInLevel is None:
            infoDf = self.classInfo.infoDf[self.classInfo.infoDf.isTerminalCluster==True]
        else:
            infoDf = self.classInfo.infoDf.query(f'zoomInLevel=={zoomInLevel}')
        percDf = infoDf.groupby('Class').apply(self.domainInfo.getDomainPercent).reset_index()
        percDf.columns = ['Class', 'domName', 'perc']
        percDf = percDf.pivot(columns='Class', index='domName', values='perc')
        percDf = percDf.fillna(0)
        percDf = percDf[percDf.apply(lambda x:np.max(x)>=0.5, axis=1)]
        percDf = percDf * 100
        Class2finalClass = self.classInfo.getClass2finalClass()
        percDf.columns = [Class2finalClass[x] for x in percDf.columns]
        percDf = percDf[sorted(Class2x, key=Class2x.get)]

        import matplotlib.pyplot as plt
        import seaborn as sns
        from palettable.cartocolors.diverging import Geyser_3
        if figSize is None:
            figsize=(9, 2)

        # ax = sns.heatmap(percDf,
        #                  cmap=Geyser_3.mpl_colormap,
        #                  linewidths=1,
        #                  linecolor='#cfcfcf',
        #                  ax=ax,
        #                  cbar_kws=dict(use_gridspec=False, location="right", shrink=1, pad=0.01, label='Percentage'))
        # ax.set_xlabel('Sub-lineage Cluster')
        # ax.set_ylabel('Protein Domain')
        # ax.tick_params(pad=-2)
        # if not (ax is None):
        #     return ax
        # fig, ax = plt.subplots()
        # cax = inset_axes(ax,
        #                  width="10%",  # width: 40% of parent_bbox width
        #                  height="100%",  # height: 10% of parent_bbox height
        #                  loc='center right',
        #                  bbox_to_anchor=(1, 0.5, 1, 1),
        #                  bbox_transform=ax.transAxes,
        #                  borderpad=0,
        #                  )
        g = sns.clustermap(data=percDf,
                           figsize=figSize,
                           cmap=Geyser_3.mpl_colormap,
                           row_cluster=False,
                           col_cluster=False,
                           linecolor='#cfcfcf',
                           linewidths=1,
                           col_colors=ClassColors,
                           colors_ratio=0.08,
                           dendrogram_ratio=0.01,
                           cbar_pos=(0.9, 0.175, 0.02, 0.7),
                           cbar_kws=dict(label='Percentage'))
        g.ax_heatmap.tick_params(axis='y', labelright=False, labelleft=True,
                                 right=False, left=True)
        g.ax_heatmap.yaxis.set_label_position("left")
        g.ax_heatmap.set_ylabel('Protein Domain')
        g.ax_heatmap.set_xlabel('Sub-lineage Cluster')
        if (ax is None) and (figSize is None):
            plt.show()
            plt.close()
        # plt.savefig('/data/home/testXT/Elements/LTR_Stream/Heat_ClusterVsSubgGenome.pdf', bbox_inches='tight')
class crossChangeValidator:
    def __init__(self, classTsvFile, domainInfo, workDir, removeCentPerc):
        self.classInfo = classInfoIntegrator.init_fromTsv(classTsvFile)
        self.classInfo.filterCentPoint(perc=removeCentPerc)
        self.domainInfo = domainInfo
        # Final PROT Tree is 41
        self.randomState = 41
        self.init_dir(workDir)
    def init_dir(self, workDir):
        self.baseD = workDir
        self.protFaD = f'{self.baseD}/protFasta'
        self.nuclFaD = f'{self.baseD}/nuclFasta'
        self.alnD = f'{self.baseD}/aln'
        constructEnv([self.baseD, self.protFaD, self.nuclFaD, self.alnD])
    def ranDomSeq(self, num, domNameSer, Classes=None, outPrefix='', zoomInLevel=None, excludeLabels=None):

        classDf = self.classInfo.infoDf.copy(deep=True)
        if Classes is not None:
            classDf = classDf[classDf.Class.isin(Classes)]
        if zoomInLevel is not None:
            classDf = classDf[classDf.zoomInLevel==zoomInLevel]
        if (Classes is None) and (zoomInLevel is None):
            classDf = classDf[classDf.isTerminalCluster==True]
        if excludeLabels is not None:
            classDf = classDf[~classDf.finalClass.isin(excludeLabels)]
        validOriIdSer = self.domainInfo.getValidOriIdSer(classDf, domNameSer)
        classDf = classDf[classDf.oriId.isin(validOriIdSer)]
        oriId2Class = {row.oriId:row.Class for ind, row in classDf.iterrows()}
        tmpOriIdDf = pd.DataFrame(dict(oriId=[oriId for oriId in oriId2Class],
                                       Class=[Class for oriId,Class in oriId2Class.items()]))
        tmpOriIdDf = tmpOriIdDf.groupby('Class').apply(lambda x: x.sample(n=min(num, x.shape[0]), replace=False, random_state=self.randomState))
        outFaList = []
        for ii, dom in domNameSer.items():
            domSeqDf = self.domainInfo.selectDomSeq(tmpOriIdDf, dom, retDf=True)
            domSeqDf['Class'] = domSeqDf.oriId.map(oriId2Class)
            domSeqDf = domSeqDf.reset_index(drop=True)
            Class2finalClass = self.classInfo.getClass2finalClass()
            domSeqDf['Class'] = domSeqDf.Class.map(Class2finalClass)
            outFa = f'{self.protFaD}/{outPrefix}.{dom}.fastq'
            outFaList.append(outFa)
            with open(outFa, 'w') as of:
                for ind, row in domSeqDf.iterrows():
                    print(f'>{row.Class}__{ind}', file=of)
                    print(row.domSeq, file=of)
        return outFaList
    def ranNucl(self, num, Classes=None, outPrefix='', zoomInLevel=None):
        outFa = f'{self.nuclFaD}/{outPrefix}.fastq'
        self.classInfo.ranNucl(num, outFa, Classes, zoomInLevel)
        return outFa
    def runPairwiseAln(self, domNameSer, Classes=None, outPrefix='', cpuNum=48, typ='dom'):
        def runAlnEachDom(dom, outPrefix, cpuNum):
            fastq = f'{self.protFaD}/{outPrefix}.{dom}.fastq'
            df = fq2df(fastq)
            df.columns = ['Class', 'seq']
            pv = p2pValidator(df, cpuNum=cpuNum)
            relList = pv.getAll2allIdentity(typ='prot')
            relList = [(df.Class.iloc[rel[0]], df.Class.iloc[rel[1]], rel[2]) for rel in relList]
            return relList
        def runNuclAln(outPrefix, cpuNum):
            fastq = f'{self.nuclFaD}/{outPrefix}.fastq'
            df = fq2df(fastq)
            df.columns = ['Class', 'seq']
            pv = p2pValidator(df, cpuNum=cpuNum)
            relList = pv.getAll2allIdentity()
            relList = [(df.Class.iloc[rel[0]], df.Class.iloc[rel[1]], rel[2]) for rel in relList]
            return relList

        if typ=='dom':
            return self.p2pSer2p2pIdenMatSer(domNameSer.apply(runAlnEachDom, outPrefix=outPrefix, cpuNum=cpuNum))
        elif typ=='nucl':
            return self.p2pAlnRelList2Mat(runNuclAln(outPrefix=outPrefix, cpuNum=48))
    @classmethod
    def p2pAlnRelList2Mat(cls, relList):
        id1s, id2s, idens = [], [], []
        visSet = set()
        for id1, id2, iden in relList:
            if id1 != id2:
                id1s.append(id1)
                id2s.append(id2)
                id1s.append(id2)
                id2s.append(id1)
                idens.append(100 - iden)
                idens.append(100 - iden)
        df = pd.DataFrame(dict(id1=id1s, id2=id2s, iden=idens))
        return df.pivot(index='id1', columns='id2', values='iden').fillna(0)
    @classmethod
    def p2pSer2p2pIdenMatSer(cls, p2pSer):
        return p2pSer.apply(cls.p2pAlnRelList2Mat)
    @classmethod
    def runMds(cls, df, n_components=1):
        from sklearn.manifold import MDS
        embedding = MDS(n_components=n_components, dissimilarity='precomputed')
        relDf = pd.DataFrame(embedding.fit_transform(df))
        relDf.columns = ['MDS_1', 'MDS_2']
        relDf.index = list(df.index)
        return relDf
    @classmethod
    def plotMdsScatter(cls, df):
        import matplotlib.pyplot as plt
        import seaborn as sns
        df = df.reset_index()
        df.columns = ['Class'] + list(df.columns[1:])
        df.Class = df.Class.apply(lambda x:'.'.join(x.split('__')[0].split('_')[1:]))
        sns.scatterplot(data=df, x='MDS_1', y='MDS_2', hue='Class', linewidth=0)
        plt.show()
        plt.close()
    @classmethod
    def plotDiffHeat(cls, df, colorDict):
        import matplotlib.pyplot as plt
        import seaborn as sns
        def toClass(x):
            return '.'.join(x.split('__')[0].split('_')[1:])
        if colorDict is not None:
            colors = list(map(lambda x:colorDict[x], list([toClass(x) for x in list(df.index)])))
            sns.clustermap(data=df,
                           row_colors=colors,
                           col_colors=colors)
            plt.show()
            plt.close()
        else:
            sns.clustermap(data=df,
                           row_cluster=False,
                           col_cluster=False)
    def pp_protDomPlot(self, Classes, domSer, colorDict, zoomInLevel=None):
        self.ranDomSeq(num=50, domNameSer=pd.Series(domSer), Classes=Classes, zoomInLevel=zoomInLevel)
        p2pIdenMatSer = self.runPairwiseAln(domNameSer=pd.Series(domSer), Classes=Classes)
        p2pMdsRelSer = p2pIdenMatSer.apply(self.runMds, n_components=2)
        p2pMdsRelSer.apply(self.plotMdsScatter)
        p2pIdenMatSer.apply(self.plotDiffHeat, colorDict=colorDict)
    @classmethod
    def runFasttree(cls, inputFa, outTree):
        softPath = '/data/home/testXT/miniconda3/envs/jup/bin/fasttree'
        os.system(f'''
            PATH={softPath}:$PATH
            fasttree {inputFa} >{outTree}
        ''')
    def plotTree(self, treeFile, outPdf):
        softPath = '/data/home/testXT/miniconda3/envs/ggtree/bin/'
        os.system(f'''
            PATH={softPath}:$PATH
            Rscript /data/home/testXT/LTR_Stream/validation/test_ggtree2.R {treeFile} {outPdf}
        ''')
