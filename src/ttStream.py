import anndata
import networkx as nx
import re
import ast
import pandas as pd
def saveWorkData(workData,outFile):
    strBranchId = []
    for ind,row in workData.obs.iterrows():
        strBranchId.append(str(row['branch_id']))

    workData.obs['branch_id'] = workData.obs['branch_id'].astype(str)
    workData.obs['branch_id_alias'] = workData.obs['branch_id_alias'].astype(str)
    try:
        workData.obs['epg_branch_node'] = workData.obs['epg_branch_node'].astype(str)
    except:
        pass

    workData.uns['epg_obj'] = None
    workData.uns['ori_epg_obj'] = None

    sepStorageArr = ['epg','seed_epg','flat_tree','seed_flat_tree','ori_epg']
    extTyp = outFile.split('.')[-1]
    outPre = re.findall(f'(.+)\.{extTyp}$',outFile)[0]
    for arr in sepStorageArr:
        nx.write_gpickle(workData.uns[arr],f'{outPre}.{arr}.pkl')
        workData.uns[arr] = None
    workData.write(outFile)
def loadWorkData(toLoadFile):
    workData = anndata.read_h5ad(toLoadFile,backed='r+')
    extTyp = toLoadFile.split('.')[-1]
    toLoadPre = re.findall(f'(.+)\.{extTyp}$', toLoadFile)[0]
    sepStorageArr = ['epg','seed_epg','flat_tree','seed_flat_tree','ori_epg']
    branch_id = []
    branch_id_alias = []
    epg_branch_node = []
    for ind,row in workData.obs.iterrows():
        branch_id.append(ast.literal_eval(row['branch_id']))
        branch_id_alias.append(ast.literal_eval((row['branch_id_alias'])))
        try:
            epg_branch_node.append(ast.literal_eval((row['epg_branch_node'])))
        except:
            pass
    workData.obs['branch_id'] = branch_id
    workData.obs['branch_id_alias'] = branch_id_alias
    try:
        workData.obs['epg_branch_node'] = epg_branch_node
    except:
        pass
    for arr in sepStorageArr:
        workData.uns[arr] = nx.read_gpickle(f'{toLoadPre}.{arr}.pkl')
    workData.obs['label'] = workData.obs['label'].astype(str)
    # workData.obs['kmeans'] = workData.obs['kmeans'].astype(str)
    return workData