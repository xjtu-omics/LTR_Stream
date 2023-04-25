import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['pdf.fonttype'] = 42
from matplotlib import pyplot as plt
from collections import defaultdict
import seaborn as sns
import sys
import pandas as pd
import yaml
with open('envConfig.yaml','r') as of:
    ttPara = yaml.safe_load(of)
    sys.path.append(ttPara['LTR_Stream'])
import ttDb
import ttUtils


baseD = sys.argv[1]
refConfigFile = sys.argv[2]
tesorterD = f'{baseD}/tesorter'
figureD = f'{baseD}/figure'
danteD = f'{baseD}/dante'
speId2speName = ttDb.getSpeId2speName(refConfigFile)

oriId2teClass = ttDb.getOriId2tesorterClass(refConfigFile,tesorterD)
oriId2modId = ttUtils.getOriId2modId(f'{danteD}/finalModSeq.tab', f'{danteD}/toCalRest.modSeq2modId.tab')
modIdSet = defaultdict(int)
for oriId in oriId2modId:
    modIdSet[oriId2modId[oriId]] = 1

oriIdList = []
ltrTypList = []
speList = []
for speId in speId2speName:
    tesOut = f'{tesorterD}/{speId}.tesorter.out.cls.tsv'
    tmpData = pd.read_table(tesOut,sep='\t')
    oriIdList = oriIdList + tmpData['#TE'].tolist()
    ltrTypList = ltrTypList + tmpData['Clade'].tolist()
    speList = speList + [speId2speName[speId] for i in range(tmpData.shape[0])]
totData = pd.concat([pd.Series(oriIdList),
                     pd.Series(ltrTypList),
                     pd.Series(speList)],
                    axis=1)
totData.columns = ['oriId','ltrTyp','spe']
countData = totData.groupby(['ltrTyp','spe']).count()
countData.reset_index(inplace=True)
countData = countData.pivot(index='ltrTyp',columns='spe',values='oriId')
countData.fillna(0,inplace=True)
countData['Total'] = countData.apply(lambda x:x.sum(),axis=1)
countData.sort_values(by='Total',ascending=False,inplace=True)
countData.to_csv(f'{figureD}/tesorter.sta.csv')
plotTopNum = 4
plotTeSet = {}
i = 0
for ind,row in countData.iterrows():
    if i>=plotTopNum:
        break
    i = i + 1
    plotTeSet[ind] = 1

modId2teClass = ['_' for _ in range(len(modIdSet))]
for oriId in oriId2modId:
    modId = oriId2modId[oriId]
    teClass = oriId2teClass[oriId]
    if teClass not in plotTeSet:
        teClass = 'Other'
    modId2teClass[modId] = teClass
with open(f'{danteD}/plotTeClass.csv','w') as of:
    for teClass in modId2teClass:
        print(teClass,file=of)
if countData.shape[0]>plotTopNum:
    restData = countData.loc[countData.index[4:],:]
    restData.loc['Other',:] = restData.sum(axis=0)
    countData = pd.concat([countData.loc[countData.index[0:4],:],
                           restData.loc[['Other'],:]],
                          axis=0)
countData.drop(columns=['Total'],inplace=True)
countData.reset_index(inplace=True)
countData = countData.melt(id_vars=['ltrTyp'],var_name='spe',value_name='count')
p = sns.barplot(data=countData,x='spe',hue='ltrTyp',y='count')
p.set_xlabel('Species')
plt.legend(title='')
plt.savefig(f'{figureD}/tesorter.sta.bar.pdf',bbox_inches='tight')
