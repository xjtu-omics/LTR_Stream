import re
from collections import defaultdict
import pandas as pd
import sys
import yaml
try:
    with open('../config/envConfig.yaml', 'r') as of:
        ttPara = yaml.safe_load(of)
        if ttPara['LTR_Stream'] not in sys.path:
            sys.path.append(ttPara['LTR_Stream'])
except:
    pass
import ttStream
import networkx as nx
import pysam

def getChrm2speId(chrm):
    return chrm.split('_')[0]
def seqName2bainfo(seqName):
    from ttlib import bainfo
    st = int(seqName.split(':')[1].split('-')[0])
    ed = int(seqName.split(':')[1].split('-')[1].split('(')[0])
    chr = seqName.split(':')[0]
    strand = re.findall('\((.)\)',seqName)[0]
    return bainfo(_chr=chr,_st=st,_ed=ed,_strand=strand)
def bainfo2seqName(bain):
    return f'{bain.chr}:{bain.st}-{bain.ed}({bain.strand})'
def getAlignedId(rawId, stOff, edOff,rev=False):
    tmpBain = seqName2bainfo(rawId)
    tmpSt = tmpBain.st+stOff-1
    tmpEd = tmpBain.st+edOff-1
    if rev:
        tmpSt = tmpBain.ed-edOff
        tmpEd = tmpBain.ed-stOff
    tmpBain.st = tmpSt
    tmpBain.ed = tmpEd
    return bainfo2seqName(tmpBain)
def seqName2len(seqName):
    st = int(seqName.split(':')[1].split('-')[0])
    ed = int(seqName.split(':')[1].split('-')[1].split('(')[0])
    return ed-st+1
def getRest2oriId(ori2restIdTabFile):
    rest2oriId = defaultdict(str)
    for line in open(ori2restIdTabFile,'r'):
        line = line.strip()
        ll = line.split('\t')
        rest2oriId[ll[1]] = ll[0]
    return rest2oriId
def getAnnotModid2Modnum(danteAnnotGffFileList):
    descript2num = defaultdict(int)
    num2descript = defaultdict(str)
    annotModid2modnum = defaultdict(str)
    oriId2annotModid = defaultdict(list)
    usedNum = 0
    classData = None
    for danteAnnotGffFile in danteAnnotGffFileList:
        if classData is None:
            classData = pd.read_table(danteAnnotGffFile, sep='\t', header=None, comment='#')
            classData.columns = ['oriId', 'soft', 'prot', 'st', 'ed', 'len', 'strand', 'dot', 'class']
        else:
            tmpData = pd.read_table(danteAnnotGffFile, sep='\t', header=None, comment='#')
            tmpData.columns = ['oriId', 'soft', 'prot', 'st', 'ed', 'len', 'strand', 'dot', 'class']
            classData = pd.concat([classData,tmpData],axis=0)

    for ind, row in classData.iterrows():
        treeInfoStr = re.findall('(Name=.+?;Final_Classification=.+?);', row['class'])[0]
        if treeInfoStr not in descript2num:
            usedNum = usedNum +1
            descript2num[treeInfoStr] = usedNum
            num2descript[usedNum] = treeInfoStr
        rev = False
        if (seqName2bainfo(row['oriId'])).strand == '-':
            rev=True
        annotModid = getAlignedId(row['oriId'],min(row['st'],row['ed']),max(row['st'],row['ed']),rev=rev)
        annotModid2modnum[annotModid] = f'a{descript2num[treeInfoStr]}'
        oriId2annotModid[row['oriId']].append(annotModid)
    return oriId2annotModid,annotModid2modnum,num2descript
def getTopNovelModSet(num,ufsId2numTabFile):
    relSet = defaultdict(int)
    numData = pd.read_table(ufsId2numTabFile,sep='\t',header=None)
    numData.columns = ['ui', 'num']
    numData.sort_values(by='num', ascending=False, inplace=True)
    i = 0
    for ind,row in numData.iterrows():
        relSet[row['ui']] = 1
        i = i+1
        if i>num:
            break
    return relSet
def getOriId2inRangeId(ltrRetRelFile,relDict = None):
    def getFormatId(myLtrRecord):
        strand = '+'
        if myLtrRecord['Strand'] == '-':
            strand = '-'
        myId = myLtrRecord['LTR_loc'].replace('..', '-')
        return f'{myId}({strand})'
    rel = defaultdict(str)
    if not (relDict is None):
        rel = relDict
    totData = pd.read_table(ltrRetRelFile, sep='\t', comment='#', header=None)
    totData.columns = ['LTR_loc', 'Category', 'Motif', 'TSD', '5_TSD', '3_TSD', 'Internal', 'Identity', 'Strand',
                       'SuperFamily', 'TE_type', 'Insertion_Time']
    for ind, row in totData.iterrows():
        oriId = getFormatId(row)
        insideBainfo = seqName2bainfo(oriId)
        insideBainfo.st = int(re.findall('^IN\:([0-9]+)?\.\..+', row['Internal'])[0])
        insideBainfo.ed = int(re.findall('^IN\:[0-9]+?\.\.(.+)?$', row['Internal'])[0])
        rel[oriId] = bainfo2seqName(insideBainfo)
    return rel
def getModSeq2modId(modSeq2modIdTabFile):
    tabData = pd.read_table(modSeq2modIdTabFile,sep='\t',header=None)
    tabData.columns = ['seq','id']
    rel = defaultdict(str)
    for ind,row in tabData.iterrows():
        rel[row['seq']] = row['id']
    return rel
def getModSeq2modNum(modSeq2numTabFile):
    tabData = pd.read_table(modSeq2numTabFile,sep='\t',header=None)
    tabData.columns = ['seq','num']
    rel = defaultdict(str)
    for ind,row in tabData.iterrows():
        rel[row['seq']] = row['num']
    return rel
def getModId2modNum(modSeq2modIdTabFile):
    rel = defaultdict(int)
    modSeq2modId = getModSeq2modId(modSeq2modIdTabFile)
    for modSeq in modSeq2modId:
        rel[modSeq2modId[modSeq]] = len(modSeq.split(','))
    return rel
def getModSeq2modId(modSeq2modIdTabFile):
    tabData = pd.read_table(modSeq2modIdTabFile,sep='\t',header=None)
    tabData.columns = ['seq','id']
    rel = defaultdict(str)
    for ind,row in tabData.iterrows():
        rel[row['seq']] = row['id']
    return rel
def getModId2modSeq(modSeq2modIdTabFile):
    modSeq2modId = getModSeq2modId(modSeq2modIdTabFile)
    rel = defaultdict(str)
    for modSeq in modSeq2modId:
        rel[modSeq2modId[modSeq]] = modSeq
    return rel
def getModId2modLen(modSeq2modIdTabFile):
    modId2modSeq = getModId2modSeq(modSeq2modIdTabFile)
    rel = defaultdict(int)
    for modId in modId2modSeq:
        rel[modId] = len(modId2modSeq[modId].split(','))
    return rel
def getModId2branchId(hdf5RelFile, adata=None):
    if adata is None:
        adata = ttStream.loadWorkData(hdf5RelFile)
    branchIdList = adata.obs['branch_id'].to_list()
    ft_node_label = nx.get_node_attributes(adata.uns['flat_tree'], 'label')
    modId2branchId = defaultdict(str)
    for i in range(len(branchIdList)):
        tmpList = [ft_node_label[branchIdList[i][0]], ft_node_label[branchIdList[i][1]]]
        tmpList.sort()
        modId2branchId[i] = f'{tmpList[0]}__{tmpList[1]}'
    return modId2branchId
def getBranchNodePair2branchId(hdf5RelFile, adata=None):
    if adata is None:
        adata = ttStream.loadWorkData(hdf5RelFile)
    branchNodePairList = list(adata.uns['flat_tree'].edges())
    ft_node_label = nx.get_node_attributes(adata.uns['flat_tree'], 'label')
    branchNodePair2branchId = defaultdict(str)
    for branchNodePair in branchNodePairList:
        tmpList = [ft_node_label[branchNodePair[0]], ft_node_label[branchNodePair[1]]]
        tmpList.sort()
        branchNodePair2branchId[branchNodePair] = f'{tmpList[0]}__{tmpList[1]}'
    return branchNodePair2branchId
def getModId2subBranchId(hdf5RelFile):
    modId2subBranchNodePair = defaultdict(str)
    adata = ttStream.loadWorkData(hdf5RelFile)
    subBranchPairList = adata.obs['epg_branch_node'].to_list()
    for i in range(len(subBranchPairList)):
        tmpList = [subBranchPairList[i][0], subBranchPairList[i][1]]
        tmpList.sort()
        modId2subBranchNodePair[i] = f'{tmpList[0]}__{tmpList[1]}'
    return modId2subBranchNodePair
def getModId2OriIdSet(oriId2modSeqTabFile, modSeq2modIdTabFile):
    modId2OriIdSet = defaultdict(lambda :defaultdict(int))
    oriId2modId = getOriId2modId(oriId2modSeqTabFile,modSeq2modIdTabFile)
    for oriId in oriId2modId:
        modId = oriId2modId[oriId]
        modId2OriIdSet[modId][oriId] = 1
    return modId2OriIdSet
def getModId2speIdSet(oriId2modSeqTabFile, modSeq2modIdTabFile):
    modId2oriIdSet = getModId2OriIdSet(oriId2modSeqTabFile, modSeq2modIdTabFile)
    modId2speIdSet = {}
    for modId in modId2oriIdSet:
        oriIdSet = modId2oriIdSet[modId]
        modId2speIdSet[modId] = set()
        for oriId in oriIdSet:
            speId = getChrm2speId(seqName2bainfo(oriId).chr)
            modId2speIdSet[modId].add(speId)
    return modId2speIdSet
def getOriId2modSeq(oriId2modSeqTabFile):
    tabData = pd.read_table(oriId2modSeqTabFile,sep='\t',header=None)
    tabData.columns = ['oriId','modSeq']
    rel = defaultdict(str)
    for ind,row in tabData.iterrows():
        rel[row['oriId']] = row['modSeq']
    return rel
def getOriId2modId(oriId2modSeqTabFile,modSeq2modIdTabFile):
    oriId2modSeq = getOriId2modSeq(oriId2modSeqTabFile)
    modSeq2modId = getModSeq2modId(modSeq2modIdTabFile)
    rel = defaultdict(int)
    for oriId in oriId2modSeq:
        rel[oriId] = modSeq2modId[oriId2modSeq[oriId]]
    return rel
def getModId2oriId(oriId2modSeqTabFile,modSeq2modIdTabFile):
    oriId2modSeq = getOriId2modSeq(oriId2modSeqTabFile)
    modSeq2modId = getModSeq2modId(modSeq2modIdTabFile)
    rel = defaultdict(int)
    for oriId in oriId2modSeq:
        rel[modSeq2modId[oriId2modSeq[oriId]]] = oriId
    return rel
def getOriId2tesorterClass(tesorterOutFile,preDict=None):
    import pandas as pd
    from collections import defaultdict
    aData = pd.read_table(tesorterOutFile,sep='\t')
    rel = defaultdict(str)
    if not (preDict is None):
        rel = preDict
    for ind,row in aData.iterrows():
        rel[row['#TE']] = row['Clade']
    return rel
def getPcoaDisMat(pcoaDisCsvFile):
    rel = pd.read_csv(pcoaDisCsvFile,header=None)
    rel = rel.iloc[:,0:3]
    rel.columns = ['x','y','z']
    return rel
def getOriId2insertTime(ltrRetFile,preDict = None):
    aData = pd.read_table(ltrRetFile,sep='\t',comment='#',header=None)
    aData.columns = ['LTR_loc','Category','Motif','TSD','5_TSD','3_TSD','Internal','Identity','Strand','SuperFamily','TE_type','Insertion_Time']
    rel = defaultdict(float)
    if not(preDict is None):
        rel = preDict
    def getFormatId(myLtrRecord):
        strand = '+'
        if myLtrRecord['Strand'] == '-':
            strand = '-'
        myId = myLtrRecord['LTR_loc'].replace('..','-')
        return f'{myId}({strand})'
    for ind,row in aData.iterrows():
        oriId = getFormatId(row)
        insertTime = row['Insertion_Time']
        rel[oriId] = insertTime
    return rel
def getNovelRotui2oriIdSet(oriId2allRotIdTabFile):
    novelRotui2oriIdSet = defaultdict(lambda: defaultdict(int))
    myData = pd.read_table(oriId2allRotIdTabFile, sep='\t', header=None)
    myData.columns = ['oriId', 'rotIdList']
    for ind,row in myData.iterrows():
        rotIdList = row['rotIdList'].split(',')
        oriId = row['oriId']
        for rotId in rotIdList:
            novelRotui2oriIdSet[int(rotId)][oriId] = 1
    return novelRotui2oriIdSet
def getBranchId2epgNodeFinalPath(hdf5RelFile, adata=None):
    if adata is None:
        adata = ttStream.loadWorkData(hdf5RelFile)
    epgTree = adata.uns['epg']
    flatTree = adata.uns['flat_tree']
    epgLeafSet = getLeafSet(epgTree)
    # This gives a terminal epg (node1, node2) and a path from leaf to root
    totBranch2epgNodePath = nx.get_edge_attributes(flatTree, 'nodes')
    finalBranch2epgNodePath = {}
    for branchNodePair in totBranch2epgNodePath:
        if any([(node in epgLeafSet) for node in branchNodePair]):
            finalBranch2epgNodePath[branchNodePair] = list(reversed(totBranch2epgNodePath[branchNodePair]))
    branchNodePair2branchId = getBranchNodePair2branchId(hdf5RelFile, adata)
    branchId2epgNodeFinalPath = {}
    for branchNodePair in finalBranch2epgNodePath:
        branchId = branchNodePair2branchId[branchNodePair]
        branchId2epgNodeFinalPath[branchId] = finalBranch2epgNodePath[branchNodePair]
    return branchId2epgNodeFinalPath
def getLeafSet(tree):
    leafSet = set()
    for n in tree.nodes:
        if tree.degree(n) == 1:
            leafSet.add(n)
    return leafSet
def getDisMat(disMatCsvFile):
    disDf = pd.read_csv(disMatCsvFile, header=None)
    return disDf
def getFaSeqNameDict(refFa, prefix):
    import pysam
    relDict = {}
    with pysam.FastxFile(refFa) as fa:
        usedNum = 0
        for entry in fa:
            usedNum += 1
            relDict[entry.name] = f'{prefix}_{usedNum}'
    return relDict
def transChrName(oriFaFile, prefix, b2a=False, chrNum=None):

    oriChr2chr = {}
    with pysam.FastxFile(oriFaFile) as oriFa:
        for ind, entry in enumerate(oriFa):
            if (not chrNum is None) and (ind>=chrNum):
                break
            oriChr2chr[entry.name] = f'{prefix}_{int(ind+1)}'
    if b2a:
        chr2oriChr = {b:a for a,b in oriChr2chr.items()}
        return chr2oriChr
    return oriChr2chr
def transJuicerStyleChr(oriFaFile, b2a=False):
    oriChr2chr = {}
    with pysam.FastxFile(oriFaFile) as oriFa:
        for ind, entry in enumerate(oriFa):
            reRel = re.findall('chr([0-9]+)$', entry.name)
            if len(reRel)>0:
                chrNum = reRel[0]
                oriChr2chr[entry.name] = chrNum
            else:
                oriChr2chr[entry.name] = str(entry.name).upper()
    if b2a:
        chr2oriChr = {b:a for a,b in oriChr2chr.items()}
        return chr2oriChr
    return oriChr2chr
def getChr2juicerStyle(oriFaFile, prefix):
    chr2oriChr = transChrName(oriFaFile, prefix, b2a=True)
    oriChr2juicerStyle = transJuicerStyleChr(oriFaFile)
    chr2juicerStyle = {chr: oriChr2juicerStyle[chr2oriChr[chr]] for chr in chr2oriChr}
    return chr2juicerStyle