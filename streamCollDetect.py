import matplotlib as mpl
mpl.use('agg')
mpl.rcParams['pdf.fonttype'] = 42
import sys
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')
import yaml
with open('envConfig.yaml','r') as of:
    ttPara = yaml.safe_load(of)
    sys.path.append(ttPara['LTR_Stream'])
import stream as st
import ttStream
import ttDb
import ttUtils

envDirs = []
baseD = sys.argv[1]
labelFile = sys.argv[2]
refConfig = sys.argv[3]
if refConfig=='None':
    refConfig = None
minNumInClusterForInitTree = int(sys.argv[4])
epgLambda = float(sys.argv[5])
epgAlpha = float(sys.argv[6])
epgMu = float(sys.argv[7])
nClusterForInitTree = int(sys.argv[8])

tesorterOutFile = None
try:
    tesorterOutFile = sys.argv[9]
except:
    pass

ltrFastaFile = None
try:
    ltrFastaFile = sys.argv[10]
except:
    pass



envDirs.append(baseD)
danteD = f'{baseD}/dante'
envDirs.append(danteD)
streamD = f'{danteD}/stream'
envDirs.append(streamD)
tesorterD = f'{baseD}/tesorter'
envDirs.append(tesorterD)
figureD = f'{baseD}/figure'
envDirs.append(figureD)

for myD in envDirs:
    if not os.path.exists(myD):
        os.system(f'''
            mkdir -p {myD}
        ''')

st.set_figure_params(dpi=80,style='white',figsize=[5.4,4.8],
                     rc={'image.cmap': 'viridis'})
adata=st.read(file_name=f'toCalRest.distanceMat.csv',file_path=danteD,workdir=streamD,delimiter=',')
st.add_cell_labels(adata,file_name=labelFile)

# st.dimension_reduction(adata,method='pcoa',feature='all',n_components=10,nb_pct=0.025)
# adata_low = st.switch_to_low_dimension(adata,n_components=3)
# st.seed_elastic_principal_graph(adata_low,outGif=f'{danteD}/kmeans_class.gif',n_clusters=50)

modSeq2modIdData = pd.read_table(f'{danteD}/toCalRest.modSeq2modId.tab',sep='\t',header=None)
modSeq2modIdData.columns = ['modSeq','modNum']
modLenList = [len(row[0].split(',')) for ind,row in modSeq2modIdData.iterrows()]

# tar1
st.ltr_dimension_reduction(adata)
adata_low = st.switch_to_low_dimension(adata,n_components=3)
# tar2
# default is 50
# for cotton 30
st.ltr_seed_elastic_principal_graph(adata_low,n_clusters=nClusterForInitTree,modLenList=modLenList,minNumInClusterForInitTree=minNumInClusterForInitTree)
st.plot_dimension_reduction(adata_low,color=['label'],n_components=3,show_graph=True,show_text=True,gif_name=f'{figureD}/kmeans_initTree.gif')
st.elastic_principal_graph(adata_low,epg_lambda=epgLambda,epg_alpha=epgAlpha,epg_mu=epgMu)
# tar3
st.plot_dimension_reduction(adata_low,color=['label'],n_components=3,show_graph=True,show_text=True,gif_name=f'{figureD}/finalResult.gif')
ttStream.saveWorkData(adata_low,f'{danteD}/tot.pgl.hdf5')
# Output finalInfo.csv
branchIdList = adata_low.obs['branch_id'].to_list()
oriId2insertTime = ttDb.getOriId2insertTime(refConfig,f'{baseD}/ltrRetriever', ltrFastaFile)
oriId2teClass = ttDb.getOriId2tesorterClass(refConfig, tesorterD, tesorterOutFile)
oriId2modId = ttUtils.getOriId2modId(f'{danteD}/finalModSeq.tab', f'{danteD}/toCalRest.modSeq2modId.tab')
pcoaData = pd.read_table(f'{danteD}/toCalRest.distanceMat.csv',header=None,sep=',')
pcoaData = pcoaData.iloc[:,0:3]

with open(f'{figureD}/finalInfo.tab','w') as of:
    print('LTR-RT_Name','Pcoa1','Pcoa2','Pcoa3','Branch','Clade','InsertTime',file=of,sep='\t')
    for oriId in oriId2insertTime:
        insertTime = oriId2insertTime[oriId]
        teClass = oriId2teClass[oriId]
        modId = oriId2modId[oriId]
        print(oriId,pcoaData.iloc[modId,0],pcoaData.iloc[modId,1],pcoaData.iloc[modId,2],branchIdList[modId],teClass,insertTime,file=of,sep='\t')





