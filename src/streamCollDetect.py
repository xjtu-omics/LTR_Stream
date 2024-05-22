import matplotlib as mpl
mpl.use('agg')
mpl.rcParams['pdf.fonttype'] = 42
import sys
import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')
import yaml
with open('../config/envConfig.yaml', 'r') as of:
    ttPara = yaml.safe_load(of)
    if ttPara['LTR_Stream'] not in sys.path:
        sys.path.append(ttPara['LTR_Stream'])
import stream as st
import ttStream
import ttDb
import ttUtils
from collections import defaultdict
from clusterValidator import clusterValidator

envDirs = []
baseD = sys.argv[1]
labelFile = sys.argv[2]
refConfig = sys.argv[3]
if refConfig=='None':
    refConfig = None
epgLambda = float(sys.argv[4])
epgAlpha = float(sys.argv[5])
epgMu = float(sys.argv[6])
maxZoomInLevel = float(sys.argv[7])

tsneLearningRate = float(sys.argv[8])
tsnePerplexity = float(sys.argv[9])
tsneEarlyExaggeration = float(sys.argv[10])
cluCentCut = float(sys.argv[11])
prefix = sys.argv[12]
tesorterOutFile = None
try:
    tesorterOutFile = sys.argv[13]
except:
    pass

ltrFastaFile = None
try:
    ltrFastaFile = sys.argv[14]
except:
    pass
ltrParaFile = sys.argv[15]

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
adata=st.read(file_name=f'toCalRest.distanceMat.csv', file_path=danteD, workdir=streamD, delimiter=',')
adata.uns['ltrParaFile'] = ltrParaFile
oriDisMat = np.array(pd.read_csv(f'{danteD}/toCalRest.distanceMat.csv', header=None))

st.add_cell_labels(adata,file_name=labelFile)



# st.dimension_reduction(adata,method='pcoa',feature='all',n_components=10,nb_pct=0.025)
# adata_low = st.switch_to_low_dimension(adata,n_components=3)
# st.seed_elastic_principal_graph(adata_low,outGif=f'{danteD}/kmeans_class.gif',n_clusters=50)

modSeq2modIdData = pd.read_table(f'{danteD}/toCalRest.modSeq2modId.tab',sep='\t',header=None)
modSeq2modIdData.columns = ['modSeq','modNum']
modLenList = [len(row[0].split(',')) for ind,row in modSeq2modIdData.iterrows()]

# tar1
st.ltr_dimension_reduction(adata, tsneLearningRate=tsneLearningRate, tsnePerplexity=tsnePerplexity, tsneEarlyExaggeration=tsneEarlyExaggeration)

adata_low = st.switch_to_low_dimension(adata,n_components=3)

modId2oriIdSet = ttUtils.getModId2OriIdSet(f'{danteD}/finalModSeq.tab', f'{danteD}/toCalRest.modSeq2modId.tab')

modId2insTimeSet = None
try:
    inRangeId2oriId = ttDb.getInRangeId2oriId(refConfig, f'{baseD}/ltrRetriever')
    # print(modId2oriIdSet[1])
    oriId2insTime = ttDb.getOriId2insertTime(refConfig,f'{baseD}/ltrRetriever')
    # print(len(oriId2insTime))
    # print(list(oriId2insTime.keys())[0])
    # print(oriId2insTime[list(oriId2insTime.keys())[0]])
    modId2insTimeSet = {}
    for modId in modId2oriIdSet:
        modId2insTimeSet[modId] = set()
        for oriId in modId2oriIdSet[modId]:
            oriId = inRangeId2oriId[oriId]
            modId2insTimeSet[modId].add(oriId2insTime[oriId])
except:
    pass

zoomInAdataDict = {}
inte = st.zoomIn(adata_low, oriDisMat, list(range(oriDisMat.shape[0])), tsneLearningRate, tsnePerplexity, tsneEarlyExaggeration, epgLambda, epgAlpha, epgMu, prefix, modId2insTimeSet, ltrFastaFile, zoomInAdataDict=zoomInAdataDict, zoomInLevel=1, maxZoomInLevel=maxZoomInLevel, cluCentCut=cluCentCut)
zoomInAdataDict = inte.reformatZoomInOrder(prefix, zoomInAdataDict)

inte.geneFinalClassLabel()

if inte.isNull():
    print('The parameters are not supported for clustering the LTR-RTs. Please check the coverLine and adjust parameters.')
    exit()
inte.writeRel()
inte.writeUserRel(prefix=prefix)
st.plotAllZoomIn(inte, zoomInAdataDict, prefix)
cv = clusterValidator(ii=inte, ltrParaFile=ltrParaFile)
cv.validateNucl()
