import argparse
from collections import defaultdict,Counter
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import pycircos
import numpy as np
import networkx as nx
from ttUtils import getOriId2modId,seqName2bainfo,getPcoaDisMat,getModId2modLen,getModId2oriId
from ttDb import getOriId2insertTime
import seaborn as sns
import sys
from sklearn.cluster import KMeans
import argparse
import yaml
with open('envConfig.yaml','r') as of:
    ttPara = yaml.safe_load(of)
    sys.path.append(ttPara['LTR_Stream'])
import ttStream
from ttlib import plotRotate3DScatterWithColor
from ttDb import getSpeId2chrNum,getSpeIdList,getSpeId2speName,getOriId2tesorterClass
from ttDb import getSpeId2speName
import ttConfiger

plotColors = ['#BF7067','#66CC96','#8AA0CF','#D1B36D','#BB65E0','#F570A3','#000000']

parser = argparse.ArgumentParser(description='Select the branch to plot.')
parser.add_argument('-w','--workDir',type=str)
parser.add_argument('-b','--branch',nargs='+',action='append',type=str)
parser.add_argument('-c','--contain',nargs='+',action='append',type=str)
parser.add_argument('-l','--lineage',nargs='+',action='append',type=str)
parser.add_argument('-o','--outPre',type=str)
parser.add_argument('-r','--refConfig',type=str)
parser.add_argument('-f','--configFile',type=str)
parser.add_argument('-t','--maxTime', type=float,default=2e6)
args = parser.parse_args()
if args.branch is None:
    args.branch = []
if args.contain is None:
    args.contain = []
if args.lineage is None:
    args.lineage = []
classNum = 0
for b in args.branch:
    if len(b)==2 or ('k' not in b[2]):
        classNum += 1
    elif 'k' in b[2]:
        classNum += int(b[2][1:])
for c in args.contain:
    if len(c) == 1 or ('k' not in c[1]):
        classNum += 1
    else:
        classNum += int(c[1][1:])
toPlotLineageList = []
for l in args.lineage:
    for ll in l:
        toPlotLineageList.append(ll)
    args.lineage = toPlotLineageList
if classNum>7:
    print('At most seven classes were supported!')
    print('Please plot the additional on other plots!')
    exit(-1)
baseD = ''
if args.workDir is not None:
    baseD = args.workDir
else:
    baseD = ttConfiger.getParaDict(args.configFile)['workDir']
refConfig = ''
if args.refConfig is not None:
    refConfig = args.refConfig
else:
    refConfig = ttConfiger.getParaDict(args.configFile)['refConfig']

danteD = f'{baseD}/dante'
refD = f'{baseD}/ref'
figureD = f'{baseD}/figure'
tesorterD = f'{baseD}/tesorter'
speId2speName = getSpeId2speName(refConfig)
class2minSeqLen = defaultdict(int)
branchList = []
legendPatchList = []
kDict = defaultdict(int)
lDict = defaultdict(int)
for tl in args.branch:
    if tl[0] > tl[1]:
        branchList.append((tl[1],tl[0]))
    else:
        branchList.append((tl[0],tl[1]))
    branchName = f'{branchList[-1][0]}__{branchList[-1][1]}'
    kDict[branchName] = 1
    lDict[branchName] = 0
    if len(tl) == 3:
        if 'k' in tl[2]:
            kDict[branchName] = int(tl[2][1:])
        elif 'l' in tl[2]:
            lDict[branchName] = int(tl[2][1:])

containList = []
for cc in args.contain:
    containList.append(cc[0])
    kDict[cc[0]] = 1
    lDict[cc[0]] = 0
    if len(cc)==2:
        if 'k' in cc[1]:
            kDict[cc[0]] = int(cc[1][1:])
        elif 'l' in cc[1]:
            lDict[cc[0]] = int(cc[1][1:])

modId2modLen = getModId2modLen(f'{danteD}/toCalRest.modSeq2modId.tab')

garc = pycircos.Garc
gcircle = pycircos.Gcircle

adata = ttStream.loadWorkData(f'{danteD}/tot.pgl.hdf5')
branchIdList = adata.obs['branch_id'].to_list()
ft_node_label = nx.get_node_attributes(adata.uns['flat_tree'],'label')
modId2branchId = defaultdict(str)
for i in range(len(branchIdList)):
    tmpList = [ft_node_label[branchIdList[i][0]],ft_node_label[branchIdList[i][1]]]
    tmpList.sort()
    if tmpList[0] > tmpList[1]:
        modId2branchId[i] = (tmpList[1],tmpList[0])
    else:
        modId2branchId[i] = (tmpList[0],tmpList[1])

toAnnClass = []
for cc in containList:
    if kDict[cc]==1:
        toAnnClass.append(cc)
    else:
        for i in range(kDict[cc]):
            toAnnClass.append(f'{cc}__{i}')
for br in branchList:
    branchName = f'{br[0]}__{br[1]}'
    if kDict[branchName] == 1:
        toAnnClass.append(branchName)
    else:
        for i in range(kDict[branchName]):
            toAnnClass.append(f'{branchName}__{i}')
plotColorDict = {toAnnClass[i]:plotColors[i] for i in range(len(toAnnClass))}
plotColorDict['other'] = '#cdcdcd'

oriId2teClass = getOriId2tesorterClass(refConfig,tesorterD)
oriId2modId = getOriId2modId(f'{danteD}/finalModSeq.tab', f'{danteD}/toCalRest.modSeq2modId.tab')
modId2oriId = getModId2oriId(f'{danteD}/finalModSeq.tab', f'{danteD}/toCalRest.modSeq2modId.tab')
br2modIdList = defaultdict(list)
_3dClassList = [0 for _ in range(len(modId2branchId))]
for i in range(len(modId2branchId)):
    branchId = modId2branchId[i]
    teclass = oriId2teClass[modId2oriId[i]]
    if len(args.lineage)==0 or (teclass in args.lineage):
        for br in branchList:
            if br == branchId and modId2modLen[i]>=lDict[f'{br[0]}__{br[1]}']:
                _3dClassList[i] = f'{br[0]}__{br[1]}'
                br2modIdList[f'{br[0]}__{br[1]}'].append(i)
        for cc in containList:
            if cc in branchId and modId2modLen[i]>=lDict[cc]:
                _3dClassList[i] = cc
                br2modIdList[cc].append(i)
    if _3dClassList[i] == 0:
        _3dClassList[i] = 'other'

pcoaDisMat = getPcoaDisMat(f'{danteD}/toCalRest.pcoaDis.csv')
for br in kDict:
    if kDict[br] > 1:
        tmpDisMat = pcoaDisMat.iloc[br2modIdList[br],0:3]
        km = KMeans(n_clusters=kDict[br], init='k-means++', random_state=42).fit(np.array(tmpDisMat))
        cluster_labels = km.labels_
        for i in range(len(cluster_labels)):
            _3dClassList[br2modIdList[br][i]] = f'{br}__{int(cluster_labels[i])}'
pcoaDisMat.reset_index(drop=True,inplace=True)
pcoaDisMat = pd.concat([pcoaDisMat,
                        pd.Series(_3dClassList)],
                       axis=1)
pcoaDisMat.columns = ['x','y','z','class']
plotRotate3DScatterWithColor(pcoaDisMat,f'{baseD}/figure/{args.outPre}.3dClass.gif',colorTyp='str',colorDict=plotColorDict)

oriId2insertTime = getOriId2insertTime(refConfig,f'{baseD}/ltrRetriever')
speIdList = getSpeIdList(refConfig)
speId2speName = getSpeId2speName(refConfig)
for spe in speIdList:
    timeList = []
    conClassList = []
    for oriId in oriId2insertTime:
        if len(args.lineage)>0 and (oriId2teClass[oriId] not in args.lineage):
            continue
        if spe in oriId and _3dClassList[oriId2modId[oriId]] != 'other':
            timeList.append(oriId2insertTime[oriId])
            conClassList.append(_3dClassList[oriId2modId[oriId]])
    toPlotData = pd.concat([pd.Series(timeList),
                            pd.Series(conClassList)],
                           axis=1)
    toPlotData.columns = ['insertTime', 'typ']
    totTyp = toPlotData['typ'].drop_duplicates()
    for typ in totTyp:
        tmpData = toPlotData.loc[toPlotData['typ']==typ,:]
        if tmpData.shape[0] < 5:
            toPlotData = toPlotData.loc[toPlotData['typ']!=typ,:]
    g = sns.displot(data=toPlotData, x='insertTime', kind="kde", hue='typ', palette=plotColorDict, legend=False,
                    linewidth=2)
    g.fig.set_size_inches((3, 1))
    plt.xlim(0, args.maxTime)
    g.ax.text(g.ax.get_xlim()[1]*0.5,g.ax.get_ylim()[1]*0.85,speId2speName[spe])
    plt.savefig(f'{figureD}/insertTime.{spe}.pdf', bbox_inches='tight')
    plt.close()
