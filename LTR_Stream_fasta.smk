import re
import os
import sys

configfile: 'envConfig.yaml'

sys.path.append(config['LTR_Stream'])

import ttConfiger
import ttDb
import ttlib

paraDict = ttConfiger.getParaDict(config['ltrParaFile'])
inputFasta = paraDict['ltrFasta']

baseD = paraDict['workDir']

minOverLapForNovelModule = 0.8
if 'minOverLapForNovelModule' in paraDict:
    minOverLapForNovelModule = float(paraDict['minOverLapForNovelModule'])

minNumInClusterForInitTree = 10
if 'minNumInClusterForInitTree' in paraDict:
    minNumInClusterForInitTree = int(paraDict['minNumInClusterForInitTree'])

nClustersForInitTree = 50
if 'nClustersForInitTree' in paraDict:
    nClustersForInitTree = int(paraDict['nClustersForInitTree'])

kForGeneticMarker = 2
if 'kForGeneticMarker' in paraDict:
    kForGeneticMarker = int(paraDict['kForGeneticMarker'])

onlyDetectFinalBranch = 'True'
if 'onlyDetectFinalBranch' in paraDict:
    onlyDetectFinalBranch = paraDict['onlyDetectFinalBranch']

maxCNSPercentage = 1
if 'maxCNSPercentage' in paraDict:
    maxCNSPercentage = paraDict['maxCNSPercentage']

epgLambda = 0.002
if 'epgLambda' in paraDict:
    epgLambda = float(paraDict['epgLambda'])

epgAlpha = 0.2
if 'epgAlpha' in paraDict:
    epgAlpha = float(paraDict['epgAlpha'])

epgMu = 0.02
if 'epgMu' in paraDict:
    epgMu = float(paraDict['epgMu'])

def prepareWorkDir():
    global danteEnv
    global statusD,testD,blastD,configD,rawRefD,refD,ltrFinderD,ltrHarvestD,ltrRetrieverD,danteD,plotD,figureD,tesorterD,subphaserD
    danteEnv = config['danteEnv']
    envDirs = []

    envDirs.append(baseD)

    statusD = f'{baseD}/status'
    envDirs.append(statusD)

    testD = f'{baseD}/test'
    envDirs.append(testD)

    blastD = f'{baseD}/blast'
    envDirs.append(blastD)

    configD = f'{baseD}/config'
    envDirs.append(configD)

    rawRefD = f'{baseD}/rawRef'
    envDirs.append(rawRefD)

    refD = f'{baseD}/ref'
    envDirs.append(refD)

    ltrFinderD = f'{baseD}/ltrFinder'
    envDirs.append(ltrFinderD)

    ltrHarvestD = f'{baseD}/ltrHarvest'
    envDirs.append(ltrHarvestD)

    ltrRetrieverD = f'{baseD}/ltrRetriever'
    envDirs.append(ltrRetrieverD)

    danteD = f'{baseD}/dante'
    envDirs.append(danteD)

    plotD = f'{baseD}/plot'
    envDirs.append(plotD)

    figureD = f'{baseD}/figure'
    envDirs.append(figureD)

    tesorterD = f'{baseD}/tesorter'
    envDirs.append(tesorterD)

    subphaserD = f'{baseD}/subphaser'
    envDirs.append(subphaserD)
    ttlib.constructEnv(envDirs)
prepareWorkDir()

ltrFasta = f'{baseD}/ltr.fa'

rule tar:
    input: f'{statusD}/tesorter.ok',
           f'{statusD}/dante.ok',
           f'{statusD}/getNovelFasta.ok',
           f'{statusD}/calModSeqDistance.ok',
           f'{statusD}/stream.ok'

rule geneLtrFasta:
    input: inputFasta
    output: f'{statusD}/geneLtrFasta.ok'
    run:
        import pysam
        with pysam.FastxFile(inputFasta) as iFa, open(ltrFasta,'w') as of:
            for te in iFa:
                print(f'>{te.name}:0-{len(te.sequence)}(+)', file=of)
                print(te.sequence, file=of)
        shell('''
            touch {output}
        ''')

rule tesorter:
    input:  f'{statusD}/geneLtrFasta.ok'
    output: f'{statusD}/tesorter.ok'
    threads:    44
    run:
        outPre = f'{tesorterD}/tesorter.out'
        workTmp = f'{tesorterD}/workTmp'
        shell('''
            rm -rf {outPre}*
            rm -rf {workTmp}
        ''')
        shell(f'''
            TEsorter -pre {outPre} -p {threads} -tmp {workTmp} -nocln {ltrFasta}
            touch {output}
        ''')
rule dante:
    input:  f'{statusD}/geneLtrFasta.ok'
    output: f'{statusD}/dante.ok'
    threads:    4
    run:
        pdbFile = f'{danteEnv}/tool-data/protein_domains/Viridiplantae_v3.0_pdb'
        classTblFile = f'{danteEnv}/tool-data/protein_domains/Viridiplantae_v3.0_class'
        outGff = f'{danteD}/ltrAnnot.gff'
        shell('''
            PATH={danteEnv}:$PATH
            dante.py -q {ltrFasta} -pdb {pdbFile} -cs {classTblFile} --domain_gff {outGff}
            touch {output}
        ''')
rule getNovelFasta:
    input:  f'{statusD}/dante.ok'
    output: f'{statusD}/getNovelFasta.ok'
    run:
        import pandas as pd
        from ttlib import bainfo,getRestRange,bainfoList2BedFile
        from ttUtils import seqName2bainfo,bainfo2seqName
        from ttDb import getOriId2inRangeId
        oriId2insideBainStr = getOriId2inRangeId(None,ltrRetrieverD, ltrFasta=ltrFasta)
        minRestLen = 200
        toCalRestBedFile = f'{danteD}/toCalRest.bed'
        toCalRestFaFile = f'{danteD}/toCalRest.fa'
        shell(f'''
            rm -rf {toCalRestFaFile}
            rm -rf {toCalRestBedFile}
            rm -rf {danteD}/ori2restId.tab
        ''')

        annotGff = f'{danteD}/ltrAnnot.gff'
        annotData = pd.read_table(annotGff,sep='\t',comment='#',header=None)
        annotData.columns = ['oriId','soft','prot','st','ed','len','strand','dot','class']
        prevOriId = ''
        tmpBainList = []
        #用于输出从oriId到toCalRestId的映射
        oriId = []
        toCalRestId = []
        toCalRestBain = []
        for ind,row in annotData.iterrows():
            if row['oriId']!=prevOriId:
                if len(tmpBainList)>0:
                    totalRange = seqName2bainfo(oriId2insideBainStr[prevOriId])
                    restRange = getRestRange(totalRange,tmpBainList)
                    toCalList = []
                    for bain in restRange:
                        if bain.ed - bain.st>minRestLen:
                            toCalRestBain.append(bain)
                            oriId.append(prevOriId)
                            toCalRestId.append(bainfo2seqName(bain))
                    tmpBainList = []
                prevOriId = row['oriId']
            oriBain = seqName2bainfo(row['oriId'])
            if oriBain.strand == '+':
                tmpBainList.append(bainfo(_chr=oriBain.chr,_st=oriBain.st+row['st'],_ed=oriBain.st+row['ed'],_strand=oriBain.strand))
            else:
                tmpBainList.append(bainfo(_chr=oriBain.chr,_st=oriBain.ed-row['ed'],_ed=oriBain.ed-row['st'],_strand=oriBain.strand))

        if len(tmpBainList) > 0:
            totalRange = seqName2bainfo(prevOriId)
            restRange = getRestRange(totalRange,tmpBainList)
            toCalList = []
            for bain in restRange:
                if bain.ed - bain.st > minRestLen:
                    toCalRestBain.append(bain)
                    oriId.append(prevOriId)
                    toCalRestId.append(bainfo2seqName(bain))
        with open(f'{danteD}/ori2restId.tab','a') as of:
            for i in range(len(oriId)):
                print(f'{oriId[i]}\t{toCalRestId[i]}',file=of)
        tmpBedFile = f'{danteD}/tmp.bed'
        bainfoList2BedFile(toCalRestBain,tmpBedFile)
        shell('''
            bedtools getfasta -fi {inputFasta} -bed {tmpBedFile} -s >>{toCalRestFaFile}
            cat {tmpBedFile} >>{toCalRestBedFile}
            rm -rf {tmpBedFile}
            touch {output}
        ''')
rule novelBlast:
    input: f'{statusD}/getNovelFasta.ok'
    output: f'{statusD}/novelBlast.ok'
    threads: 48
    run:
        fastaF = f'{danteD}/toCalRest.fa'
        blastWorkD = f'{danteD}/blast'
        ttlib.constructEnv([blastWorkD])
        blastDb = f'{blastWorkD}/toCalRest'
        outF = f'{danteD}/toCalRest.blastSelf.tab'
        shell('''
                makeblastdb -in {fastaF} -dbtype nucl -title tot.ref.fa.pass.list -parse_seqids -out {blastDb}
                blastn -query {fastaF} -db {blastDb} -evalue 1e-10 -out {outF} -num_threads {threads} -outfmt 6 -max_target_seqs 100
                touch {output}
        ''')
rule selectNovelUnits:
    input:  f'{statusD}/novelBlast.ok'
    output: f'{statusD}/selectNovelUnits.ok'
    run:
        import pandas as pd
        from collections import defaultdict
        from ttlib import ttUFS
        from ttlib import bainfo,overlapMin,overlapA
        from ttUtils import getAlignedId, seqName2bainfo, seqName2len, bainfo2seqName
        minUnitLen = 100
        minMinOverlap = minOverLapForNovelModule
        blastRelF = f'{danteD}/toCalRest.blastSelf.tab'
        alignedIdSet = defaultdict(int)
        alignedId2originalId = defaultdict(str)
        for line in open(blastRelF,'r'):
            line = line.strip().split('\t')
            id1 = line[0]
            id2 = line[1]
            if id1==id2:
                continue
            st1 = int(line[6])
            ed1 = int(line[7])
            if st1>ed1:
                ed1,st1 = st1,ed1
            st2 = int(line[8])
            ed2 = int(line[9])
            if st2>ed2:
                ed2,st2 = st2,ed2
            alignedId1 = getAlignedId(id1,st1,ed1,rev = (seqName2bainfo(id1).strand=='-'))
            alignedId2 = getAlignedId(id2,st2,ed2,rev = (seqName2bainfo(id2).strand=='-'))
            alignedIdSet[alignedId1] = 1
            alignedIdSet[alignedId2] = 1
            alignedId2originalId[alignedId1] = id1
            alignedId2originalId[alignedId2] = id2
        alignedBainList = []
        for id in alignedIdSet:
            if seqName2len(id)<minUnitLen:
                continue
            alignedBainList.append(seqName2bainfo(id))
        alignedBainList.sort()
        originalIdList = []
        unitNum = len(alignedBainList)
        ufsId2alignedId = defaultdict(str)
        ufsId2originalId = defaultdict(str)
        for i in range(unitNum):
            alnId = bainfo2seqName(alignedBainList[i])
            if alnId not in alignedId2originalId:
                print(f'{alnId} not in alignedId2originalId')
                exit(-1)
            originalIdList.append(alignedId2originalId[alnId])
            if alnId not in alignedIdSet:
                print(f'{alnId} not in !')
                exit(-1)
            alignedIdSet[alnId] = i
            ufsId2alignedId[i] = alnId
            ufsId2originalId[i] = alignedId2originalId[alnId]
        myUfs = ttUFS(unitNum)
        for line in open(blastRelF,'r'):
            line = line.strip().split('\t')
            id1 = line[0]
            id2 = line[1]
            if id1 == id2:
                continue
            st1 = int(line[6])
            ed1 = int(line[7])
            if st1 > ed1:
                ed1, st1 = st1, ed1
            st2 = int(line[8])
            ed2 = int(line[9])
            if st2 > ed2:
                ed2, st2 = st2, ed2
            alignedId1 = getAlignedId(id1,st1,ed1,rev=(seqName2bainfo(id1).strand=='-'))
            alignedId2 = getAlignedId(id2,st2,ed2,rev=(seqName2bainfo(id2).strand=='-'))
            myUfs.union(alignedIdSet[alignedId1],alignedIdSet[alignedId2])
        for i in range(unitNum):
            for j in range(i+1,unitNum):
                if alignedBainList[j].st>alignedBainList[i].ed:
                    break
                if overlapMin(alignedBainList[i],alignedBainList[j])>minMinOverlap:
                    myUfs.union(i,j)
        ufsId2size = myUfs.getRotIdSet()

        with open(f'{danteD}/toCalRest.alignedInfo2ufsId.tab','w') as of:
            for i in range(len(alignedBainList)):
                print(bainfo2seqName(alignedBainList[i]),myUfs.find_set(i),file=of,sep='\t')

        with open(f'{danteD}/toCalRest.conserveUnits.tab','w') as of:
            for ui in ufsId2size:
                print(bainfo2seqName(alignedBainList[ui]),ufsId2size[ui],originalIdList[ui],alignedBainList[ui].ed-alignedBainList[ui].st,sep='\t',file=of)

        with open(f'{danteD}/toCalRest.ufsId2num.tab','w') as of:
            for ui in ufsId2size:
                print(ui,ufsId2size[ui],sep='\t',file=of)

        originalId2rotuiSet = defaultdict(lambda :defaultdict(int))
        alignedId2rotui = defaultdict(int)
        for i in range(unitNum):
            rotui = myUfs.find_set(i)
            originalId2rotuiSet[ufsId2originalId[i]][rotui] = 1
            alignedId2rotui[bainfo2seqName(alignedBainList[i])] = rotui
        with open(f'{danteD}/toCalRest.originalId2rotIdList.tab','w') as of:
            for ogid in originalId2rotuiSet:
                print(ogid,end='\t',file=of)
                toPrintList = list(originalId2rotuiSet[ogid].keys())
                toPrintList.sort()
                toPrintNum = len(toPrintList)
                for i in range(toPrintNum):
                    endFlag=','
                    if i==toPrintNum-1:
                        endFlag='\n'
                    print(toPrintList[i],end=endFlag,file=of)

        originalId2alignedIdInfoList = defaultdict(list)

        tmpOriIdList = []
        tmpRotuiList = []
        tmpStList = []
        tmpEdList = []
        for alnId in alignedId2rotui:
            rotui = alignedId2rotui[alnId]
            tmpBainfo = seqName2bainfo(alnId)
            st = tmpBainfo.st
            ed = tmpBainfo.ed
            oriId = alignedId2originalId[alnId]
            tmpOriIdList.append(oriId)
            tmpRotuiList.append(rotui)
            tmpStList.append(st)
            tmpEdList.append(ed)
        toWriteCsvData = pd.concat([pd.Series(tmpOriIdList),
                                    pd.Series(tmpRotuiList),
                                    pd.Series(tmpStList),
                                    pd.Series(tmpEdList)],
                                    axis = 1)
        toWriteCsvData.columns = ['oriId','rotui','st','ed']
        toWriteCsvData.drop_duplicates(['oriId','rotui'],'first',inplace=True)
        toWriteCsvData.to_csv(f'{danteD}/toCalRest.totoalInfo.csv')
        shell('''
            touch {output}
        ''')
rule getNovelModelSeq:
    input:  f'{statusD}/selectNovelUnits.ok'
    output: f'{statusD}/getNovelModelSeq.ok'
    run:
        import pandas as pd
        from collections import defaultdict
        topNum = 300
        ufsId2numFile = f'{danteD}/toCalRest.ufsId2num.tab'
        numData = pd.read_table(ufsId2numFile,sep='\t',header=None)
        numData.columns = ['rotUi','num']
        numData.sort_values(by='num',ascending=False,inplace=True)
        rotIdSet = defaultdict(int)
        i = 0
        for ind,row in numData.iterrows():
            if i>topNum:
                break
            rotIdSet[row['rotUi']] = 1
            i = i+1
        totalInfoFile = f'{danteD}/toCalRest.totoalInfo.csv'
        infoData = pd.read_csv(totalInfoFile)
        infoData.sort_values(by='st',inplace=True)
        oriId2oriIdList = defaultdict(list)
        for ind,row in infoData.iterrows():
            if row['rotui'] in rotIdSet:
                oriId2oriIdList[row['oriId']].append(row['rotui'])
        with open(f'{danteD}/toCalRest.oriId2modelSeq.tab','w') as of:
            for oriId in oriId2oriIdList:
                print(oriId,end='\t',file=of)
                tmpList = oriId2oriIdList[oriId]
                num = len(tmpList)
                for i in range(num):
                    endFlag = ','
                    if i == num-1:
                        endFlag='\n'
                    print(tmpList[i],end=endFlag,file=of)
        shell('''
            touch {output}
        ''')
rule staNovelSelectedNumVsCovered:
    input: f'{statusD}/getNovelModelSeq.ok'
    output: f'{statusD}/staNovelSelectedNumVsCovered.ok'
    run:
        import matplotlib as mpl

        mpl.use('agg')
        mpl.rcParams['pdf.fonttype'] = 42
        from matplotlib import pyplot as plt
        import pandas as pd
        import seaborn as sns
        from collections import defaultdict

        numData = pd.read_table(f'{danteD}/toCalRest.ufsId2num.tab',sep='\t',header=None)
        numData.columns = ['ui', 'num']
        numData.sort_values(by='num',axis=0,ascending=False,inplace=True)
        conserveUnitData = pd.read_table(f'{danteD}/toCalRest.originalId2rotIdList.tab',sep='\t',header=None)
        conserveUnitData.columns = ['orId', 'ui']
        ui2orIdList = defaultdict(lambda: defaultdict(int))
        totOrIdNum = conserveUnitData.shape[0]
        for ind, row in conserveUnitData.iterrows():
            uiList = row['ui'].split(',')
            myOrId = row['orId']
            for ui in uiList:
                intUi = int(ui)
                ui2orIdList[intUi][myOrId] = 1
        preSum2num = defaultdict(int)
        orId2preSum = defaultdict(int)
        toPlotDict = defaultdict(list)
        for ind, row in numData.iterrows():
            ui = row['ui']
            orIdSet = ui2orIdList[ui]
            for orId in orIdSet:
                myPreSum = orId2preSum[orId] + 1
                if myPreSum <= 10:
                    orId2preSum[orId] = myPreSum
                    preSum2num[myPreSum] = preSum2num[myPreSum] + 1
            for i in range(1,11):
                toPlotDict[i].append(preSum2num[i] / totOrIdNum)

        xList = []
        colorList = []
        yList = []
        totNum = len(toPlotDict[1])
        for i in range(1,11):
            xList = xList + list(range(1,totNum + 1))
            colorList = colorList + ([i] * totNum)
            yList = yList + toPlotDict[i]
        plotData = pd.concat([pd.Series(xList),
                              pd.Series(yList),
                              pd.Series(colorList)],
            axis=1)
        plotData.columns = ['x', 'y', 'lineColor']
        plotData = plotData[plotData['x'] < 500]
        sns.lineplot(x='x',y='y',hue='lineColor',data=plotData)
        plt.savefig(f'{figureD}/toCalRest.coverLine.pdf',bbox_inches='tight')
        plt.close()
        tmpData = plotData.loc[plotData['lineColor'] == 1, :]
        tmpData = tmpData.sort_values(by='x')
        if maxCNSPercentage < 1:
            with open(f'{danteD}/novelMoludesUsedNum.tab','w') as of:
                tmpFlag = False
                for ind, row in tmpData.iterrows():
                    if (row['y'] > 0.7):
                        tmpFlag = True
                        print(row['x'],file=of)
                        break
                if not tmpFlag:
                    print(tmpData.shape[0],file=of)
        else:
            with open(f'{danteD}/novelMoludesUsedNum.tab','w') as of:
                print('300',file=of)

        shell(f"""touch {output}""")
rule mergeModules:
    input:  f'{statusD}/staNovelSelectedNumVsCovered.ok'
    output: f'{statusD}/mergeModules.ok'
    run:
        from ttUtils import getRest2oriId,getAnnotModid2Modnum,getTopNovelModSet,seqName2bainfo,bainfo2seqName
        from collections import defaultdict
        import pandas as pd
        from ttlib import bainfo

        topNovelModUsed = -1
        for line in open(f'{danteD}/novelMoludesUsedNum.tab'):
            line = line.strip()
            topNovelModUsed = int(float(line))
        novelModSet = getTopNovelModSet(topNovelModUsed,f'{danteD}/toCalRest.ufsId2num.tab')
        danteGffList = [f'{danteD}/ltrAnnot.gff']
        ori2annot, annot2modnum, modnum2descript = getAnnotModid2Modnum(danteGffList)

        with open(f'{danteD}/toCalRest.annotId2descript.tab','w') as of:
            for modnum in modnum2descript:
                print(modnum,modnum2descript[modnum],sep='\t',file=of)

        rest2ori = getRest2oriId(f'{danteD}/ori2restId.tab')

        ori2all = ori2annot
        all2modnum = annot2modnum

        novelData = pd.read_csv(f'{danteD}/toCalRest.totoalInfo.csv')
        for ind,row in novelData.iterrows():
            if row['rotui'] not in novelModSet:
                continue
            oriId = rest2ori[row['oriId']]
            tmpBain = seqName2bainfo(oriId)
            novelId = bainfo2seqName(bainfo(_chr=tmpBain.chr,_st=row['st'],_ed=row['ed'],_strand=tmpBain.strand))
            ori2all[oriId].append(novelId)
            all2modnum[novelId] = f'b{row["rotui"]}'
        oriId2mergeMod = defaultdict(str)
        for oriId in ori2all:
            bainList = []
            modList = []
            seqNameList = ori2all[oriId]
            for seqName in seqNameList:
                bainList.append(seqName2bainfo(seqName))
            bainList.sort()
            modSeq = ''
            i = 0
            for bain in bainList:
                seqName = bainfo2seqName(bain)
                if seqName not in all2modnum:
                    print(f'{seqName} not in all2modnum!')
                    exit(-1)
                sepC = ','
                if i==0:
                    sepC = ''
                modSeq = f'{modSeq}{sepC}{all2modnum[seqName]}'
                i = i+1
            oriId2mergeMod[oriId] = modSeq
        with open(f'{danteD}/finalModSeq.tab','w') as of:
            for oriId in oriId2mergeMod:
                print(oriId,oriId2mergeMod[oriId],sep='\t',file=of)
        shell(f'''
            touch {output}
        ''')
rule getTeclassEachModId:
    input: f'{statusD}/tesorter.ok',
           f'{statusD}/mergeModules.ok'
    output: f'{statusD}/getTeclassEachModId.ok'
    run:
        import ttUtils
        from collections import defaultdict
        import pandas as pd
        tesorterOutFile = f'{tesorterD}/tesorter.out.cls.tsv'
        oriId2teClass = ttDb.getOriId2tesorterClass(None, None, tesorterOutFile)
        oriId2modId = ttUtils.getOriId2modId(f'{danteD}/finalModSeq.tab',f'{danteD}/toCalRest.modSeq2modId.tab')
        modIdSet = defaultdict(int)
        for oriId in oriId2modId:
            modIdSet[oriId2modId[oriId]] = 1
        modId2teClass = ['_' for _ in range(len(modIdSet))]
        tmpData = pd.read_table(tesorterOutFile, sep='\t')
        ltrTypSer = tmpData['Clade'].value_counts().index.tolist()
        plotTeSet = {}
        for i in range(min(3, len(ltrTypSer))):
            plotTeSet[ltrTypSer[i]] = 1
        for oriId in oriId2modId:
            modId = oriId2modId[oriId]
            teClass = oriId2teClass[oriId]
            if teClass not in plotTeSet:
                teClass = 'Other'
            modId2teClass[modId] = teClass
        with open(f'{danteD}/plotTeClass.csv','w') as of:
            for teClass in modId2teClass:
                print(teClass,file=of)
        shell(f"""
            touch {output}
        """)
rule calModSeqDistancePre:
    input: f'{statusD}/mergeModules.ok'
    output: f'{statusD}/calModSeqDistancePre.ok'
    run:
        import pandas as pd
        from collections import defaultdict
        import pickle

        modSeqFile = f'{danteD}/finalModSeq.tab'
        modSeqData = pd.read_table(modSeqFile,sep='\t',header=None)
        modSeqData.columns = ['oriId', 'modSeq']
        modSeq2oriId = defaultdict(lambda: defaultdict(str))
        modSeq2num = defaultdict(int)
        modSeq2modId = defaultdict(int)
        modSeqSet = defaultdict(int)
        modSeqList = []
        idUsed = -1
        for ind, row in modSeqData.iterrows():
            modSeq2oriId[row['modSeq']][row['oriId']] = 1
            modSeq2num[row['modSeq']] = modSeq2num[row['modSeq']] + 1
            if row['modSeq'] not in modSeqSet:
                modSeqSet[row['modSeq']] = 1
                idUsed = idUsed + 1
                modSeq2modId[row['modSeq']] = idUsed
                modSeqList.append(row['modSeq'])
        with open(f'{danteD}/toCalRest.modSeq2num.tab','w') as of:
            for modSeq in modSeq2num:
                print(modSeq,modSeq2num[modSeq],sep='\t',file=of)
        with open(f'{danteD}/toCalRest.modSeq2modId.tab','w') as of:
            for modSeq in modSeq2modId:
                print(modSeq,modSeq2modId[modSeq],sep='\t',file=of)

        modSeqListDataFile = f'{danteD}/toCalRest.modSeqList.pkl'
        with open(modSeqListDataFile,'wb') as of:
            pickle.dump(modSeqList,of)

        shell('''
             touch {output}
         ''')
rule calModSeqDistance:
    input: f'{statusD}/calModSeqDistancePre.ok'
    output: f'{statusD}/calModSeqDistance.ok'
    threads: 48
    run:
        shell(f"""
             python calModSeqDistance.py {danteD} {threads} {f'{danteD}/toCalRest.modSeqList.pkl'}
             touch {output}
         """)
rule pcoaPlot:
    input:  f'{statusD}/calModSeqDistance.ok'
    output: f'{statusD}/pcoaPlot.ok'
    run:
        from skbio.stats.ordination import pcoa
        import pandas as pd
        import matplotlib as mpl
        from ttlib import plotRotate3DScatter
        mpl.use('agg')
        mpl.rcParams['pdf.fonttype'] = 42
        from matplotlib import pyplot as plt
        import seaborn as sns
        from ttUtils import getModId2modNum

        twoDimensionData = None
        twoDataFile = f'{danteD}/toCalRest.twoDimensionData.pcoa.csv'
        disData = pd.read_csv(f'{danteD}/toCalRest.distanceMat.csv',header=None)
        pcoaRel = pcoa(disData)
        pcoaData = pcoaRel.samples
        pcoaData.to_csv(f'{danteD}/toCalRest.pcoaDis.csv',header=None,index=None)
        plotRotate3DScatter(pcoaData,f'{figureD}/pcoaDistance.3d.gif')
        pcoaData = pcoaData.iloc[:,0:3]
        pcoaData.columns = ['x','y','z']
        print(type(pcoaData.iloc[0,0]))
        print(pcoaData.head())
        #twoDimensionData.plot()
        modId2modNum = getModId2modNum(f'{danteD}/toCalRest.modSeq2num.tab',f'{danteD}/toCalRest.modSeq2modId.tab')
        modIdList = []
        modNumList = []
        for modId in modId2modNum:
            modIdList.append(modId)
            modNumList.append(modId2modNum[modId])
        sizeData = pd.concat([pd.Series(modIdList),
                              pd.Series(modNumList)],
                             axis=1)
        sizeData.columns = ['id','num']
        sizeData.sort_values(by='id',inplace=True)
        print(pcoaData['x'].head())
        pcoaData.reset_index(drop=True,inplace=True)
        sizeData.reset_index(drop=True,inplace=True)
        twoDimensionData = pd.concat([pcoaData['x'],
                                      pcoaData['y'],
                                      sizeData['num']],
                                      axis = 1)
        twoDimensionData.columns = ['x','y','num']
        plotData = twoDimensionData[twoDimensionData['num']>1]
        sns.scatterplot(data=plotData,x='x',y='y',hue='num',size='num')
        plt.savefig(f'{figureD}/twoDimensionData.pcoa.pdf',height=10,width=10)
        shell('''
            touch {output}
        ''')
rule stream:
    input:  f'{statusD}/pcoaPlot.ok',f'{statusD}/getTeclassEachModId.ok'
    output: f'{statusD}/stream.ok'
    threads: 48
    run:
        tesorterOutFile = f'{tesorterD}/tesorter.out.cls.tsv'
        with open(f'{danteD}/streamLabels.csv','w') as of:
            for _ in open(f'{danteD}/toCalRest.pcoaDis.csv','r'):
                print('LTR_RT',file=of)
        shell(f'''
            python streamCollDetect.py {baseD} {danteD}/plotTeClass.csv None \
                                       {minNumInClusterForInitTree} {epgLambda} {epgAlpha} {epgMu} {nClustersForInitTree} \
                                       {tesorterOutFile} {ltrFasta}
            touch {output}
        ''')
