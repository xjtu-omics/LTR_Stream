import numpy as np
import pandas as pd
from pandas.api.types import (
    is_string_dtype,
    is_numeric_dtype
)
import anndata as ad
import networkx as nx
import re
from statsmodels.nonparametric.smoothers_lowess import lowess
from statsmodels.sandbox.stats.multicomp import multipletests
import seaborn as sns
import pylab as plt
import plotly.graph_objects as go
import plotly.express as px
import multiprocessing
import os
from sklearn.decomposition import PCA as sklearnPCA
# from sklearn import preprocessing
from sklearn.manifold import (
    LocallyLinearEmbedding,
    TSNE,
    SpectralEmbedding,
    MDS
)
from sklearn.cluster import (
    SpectralClustering,
    AffinityPropagation,
    KMeans,
    AgglomerativeClustering
)
import hdbscan
from sklearn.metrics.pairwise import (
    # pairwise_distances_argmin_min,
    pairwise_distances,
    euclidean_distances,
    cosine_distances
)
from sklearn.metrics import(
    silhouette_score
)
import matplotlib as mpl
import matplotlib.patches as Patches
from matplotlib.patches import Polygon
from mpl_toolkits.mplot3d import Axes3D
import umap
from copy import deepcopy
import itertools
from scipy.spatial import (
    # distance,
    cKDTree,
    # KDTree
)
import math
# mpl.use('Agg')
from scipy import stats
from scipy.stats import (
    spearmanr,
    mannwhitneyu,
    # gaussian_kde,
    # kruskal
)
# from scipy.sparse import (
#     issparse,
#     lil_matrix,
#     csr_matrix
# )
from slugify import slugify
from decimal import *
import matplotlib.gridspec as gridspec
import pickle
import gzip
import shutil
import json

import pysam

from itertools import product
from rpy2.robjects.packages import importr
from rpy2.robjects import r as R
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

from .extra import *
# scikit_posthocs is currently not available in conda system.
# We will update it once it can be installed via conda.
# import scikit_posthocs as sp
from .scikit_posthocs import posthoc_conover

import ttUtils
from infoIntegrator import infoIntegrator

os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'


def set_figure_params(context='notebook',style='white',palette='deep',font='sans-serif',font_scale=1.1,color_codes=True,
                      dpi=80,dpi_save=150,figsize=[5.4, 4.8],rc=None):
    """ Set global parameters for figures. Modified from sns.set()
    Parameters
    ----------
    context : string or dict
        Plotting context parameters, see seaborn :func:`plotting_context
    style: `string`,optional (default: 'white')
        Axes style parameters, see seaborn :func:`axes_style`
    palette : string or sequence
        Color palette, see seaborn :func:`color_palette`
    font_scale: `float`, optional (default: 1.3)
        Separate scaling factor to independently scale the size of the font elements.        
    color_codes : `bool`, optional (default: True)
        If ``True`` and ``palette`` is a seaborn palette, remap the shorthand
        color codes (e.g. "b", "g", "r", etc.) to the colors from this palette.
    dpi: `int`,optional (default: 80)
        Resolution of rendered figures.
    dpi_save: `int`,optional (default: 150)
        Resolution of saved figures.
    rc: `dict`,optional (default: None)
        rc settings properties. Parameter mappings to override the values in the preset style.
        Please see https://matplotlib.org/tutorials/introductory/customizing.html#a-sample-matplotlibrc-file
    """
#     mpl.rcParams.update(mpl.rcParamsDefault)
    sns.set(context=context,style=style,palette=palette,font=font,font_scale=font_scale,color_codes=color_codes,
            rc={'figure.dpi':dpi,
                'savefig.dpi':dpi_save,
                'figure.figsize':figsize,
                'image.cmap': 'viridis',
                'lines.markersize':6,
                'legend.columnspacing':0.1,
                'legend.borderaxespad':0.1,
                'legend.handletextpad':0.1,
                'pdf.fonttype':42,})
    if(rc is not None):
        assert isinstance(rc,dict),"rc must be dict"  
        for key, value in rc.items():
            if key in mpl.rcParams.keys():
                mpl.rcParams[key] = value
            else:
                raise Exception("unrecognized property '%s'" % key)

def set_workdir(adata,workdir=None):
    """Set working directory.
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.           
    workdir: `float`, optional (default: None)
        Working directory. If it's not specified, a folder named 'stream_result' will be created under the current directory
    **kwargs: additional arguments to `Anndata` reading functions
   
    Returns
    -------
    updates `adata` with the following fields and create a new working directory if it doesn't existing.
    workdir: `str` (`adata.uns['workdir']`,dtype `str`)
        Working directory.  
    """       
    if(workdir==None):
        workdir = os.path.join(os.getcwd(), 'stream_result')
        print("Using default working directory.")
    if(not os.path.exists(workdir)):
        os.makedirs(workdir)
    adata.uns['workdir'] = workdir
    print('Saving results in: %s' % workdir)

def read(file_name,file_path=None,file_format=None,delimiter='\t',workdir=None, file_sample='barcodes.tsv',file_feature='genes.tsv',**kwargs):
    """Read gene expression matrix into anndata object.
    
    Parameters
    ----------
    file_name: `str`
        Input data file name.
    file_path: `str`, optional (default: None)
        File path. Empty string by default
    file_format: `str`, optional (default: None)
        File format. currently supported file formats: 'tsv','txt','tab','data','csv','mtx','h5ad','pklz','pkl'. If None, file_format will be inferred from the file extension.
    delimiter: `str`, optional (default: '\t')
        Delimiter to use.    
    workdir: `float`, optional (default: None)
        Working directory. If it's not specified, a folder named 'stream_result' will be created under the current directory
    file_sample:`str`, optional (default: 'barcodes.tsv')
        Sample file name. Only valid when file_format = 'mtx'
    file_feature:`str`, optional (default: 'genes.tsv')
        Feature file name. Only valid when file_format = 'mtx'
    **kwargs: additional arguments to `Anndata` reading functions
   
    Returns
    -------
    AnnData object
    """       
    if(file_format is None):
        file_format = get_extension(file_name)
    assert (file_format in ['tsv','txt','tab','data','csv','mtx','h5ad','pkl','pklz']),"file_format must be one of ['tsv','txt','tab','data','csv','mtx','h5ad','pkl','pklz']"
    if(file_path is None):
        file_path = ''
    _fp = lambda f:  os.path.join(file_path,f)     
    if(file_format in ['tsv','txt','tab','data']):
        adata = ad.read_text(_fp(file_name),delimiter=delimiter,**kwargs).T
        # adata.raw = adata        
    if(file_format == 'csv'):
        adata = ad.read_csv(_fp(file_name),delimiter=delimiter,**kwargs).T
        # adata.raw = adata
    if(file_format == 'mtx'):
        adata = ad.read_mtx(_fp(file_name),**kwargs).T 
        adata.X = adata.X
        adata.obs_names = pd.read_csv(_fp(file_sample), header=None)[0]
        features = pd.read_csv(_fp(file_feature), header=None, sep='\t')
        features.index = features.index.astype('str')
        adata.var = features
    if(file_format == 'h5ad'):
        adata = ad.read_h5ad(_fp(file_name),**kwargs)
    if(file_format == 'pkl'):
        if file_name.split('.')[-1]=='gz':
            f = gzip.open(_fp(file_name), 'rb')
            adata = pickle.load(f)
            f.close() 
        else:
            f = open(_fp(file_name), 'rb')
            adata = pickle.load(f)
            f.close()    
    if(file_format == 'pklz'):
        f = gzip.open(_fp(file_name), 'rb')
        adata = pickle.load(f)
        f.close()
        
    # if(not issparse(adata.X)):
    #     adata.X = csr_matrix(adata.X)
    
    if('workdir' not in adata.uns_keys()):
        set_workdir(adata,workdir=workdir)
    else:
        print("Working directory is already specified as '%s' " % adata.uns['workdir'])
        print("To change working directory, please run set_workdir(adata,workdir=new_directory)")
    return adata

def write(adata,file_name=None,file_path=None,file_format=None):
    """Write Anndate object to file
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix. 
    file_name: `str`, optional (default: None)
        File name. If it's not specified, a file named 'stream_result.pkl' will be created 
        under the working directory
    file_path: `str`, optional (default: '')
        File path. If it's not specified, it's set to working directory
    file_format: `str`, optional (default: None)
        File format. By default it's compressed pickle file. Currently two file formats are supported:
        'pklz': compressed pickle file
        'pkl': pickle file
    """

    if(file_name is None):
        file_name = 'stream_result.pkl'
    if(file_format is None):
        file_format = get_extension(file_name)
    assert (file_format in ['pkl','pklz']),"file_format must be one of ['pkl','pklz']"
    if(file_path is None):
        file_path = adata.uns['workdir']
    
    if(file_format == 'pklz'):
        f = gzip.open(os.path.join(file_path,file_name), 'wb')
        pickle.dump(adata, f, protocol=pickle.HIGHEST_PROTOCOL)
        f.close()  
    if(file_format == 'pkl'):
        f = open(os.path.join(file_path,file_name), 'wb')
        pickle.dump(adata, f, protocol=pickle.HIGHEST_PROTOCOL)
        f.close()            

def add_metadata(adata,file_name,delimiter='\t',file_path=''):
    """Add metadata.

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    file_path: `str`, optional (default: '')
        The file path of metadata file.
    file_name: `str`, optional (default: None)
        The file name of metadata file. 

    Returns
    -------
    updates `adata` with the following fields.
    label: `pandas.core.series.Series` (`adata.obs['label']`,dtype `str`)
        Array of #observations that stores the label of each cell.
    label_color: `pandas.core.series.Series` (`adata.obs['label_color']`,dtype `str`)
        Array of #observations that stores the color of each cell (hex color code).
    label_color: `dict` (`adata.uns['label_color']`,dtype `str`)
        Array of #observations that stores the color of each cell (hex color code). 

    updates `adata.obs` with additional columns in metadata file.
    """
    _fp = lambda f:  os.path.join(file_path,f)
    df_metadata = pd.read_csv(_fp(file_name),sep=delimiter,index_col=0)
    if('label' not in df_metadata.columns):
        print("No column 'label' found in metadata, \'unknown\' is used as the default cell labels")
        df_metadata['label'] = 'unknown'
    if('label_color' in df_metadata.columns):
        adata.uns['label_color'] = pd.Series(data=df_metadata.label_color.tolist(),index=df_metadata.label.tolist()).to_dict()
    else:
        print("No column 'label_color' found in metadata, random color is generated for each cell label")
        labels_unique = df_metadata['label'].unique()
        if(len(labels_unique)==1):
            adata.uns['label_color'] = {labels_unique[0]:'gray'}
        else:
            list_colors = sns.color_palette("hls",n_colors=len(labels_unique)).as_hex()
            adata.uns['label_color'] = {x:list_colors[i] for i,x in enumerate(labels_unique)}
        df_metadata['label_color'] = ''
        for x in labels_unique:
            id_cells = np.where(df_metadata['label']==x)[0]
            df_metadata.loc[df_metadata.index[id_cells],'label_color'] = adata.uns['label_color'][x]
    adata.obs = df_metadata.loc[adata.obs.index,:]
    return None


def add_cell_labels(adata,file_path='',file_name=None):
    """Add cell lables.

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    file_path: `str`, optional (default: '')
        The file path of cell label file.
    file_name: `str`, optional (default: None)
        The file name of cell label file. If file_name is not specified, 'unknown' is added as the label for all cells.
        

    Returns
    -------
    updates `adata` with the following fields.
    label: `pandas.core.series.Series` (`adata.obs['label']`,dtype `str`)
        Array of #observations that stores the label of each cell.
    """
    _fp = lambda f:  os.path.join(file_path,f)
    if(file_name!=None):
        df_labels = pd.read_csv(_fp(file_name),sep='\t',header=None,index_col=None,names=['label'],
                                dtype=str,compression= 'gzip' if file_name.split('.')[-1]=='gz' else None)
        df_labels['label'] = df_labels['label'].str.replace('/','-')        
        df_labels.index = adata.obs_names
        adata.obs['label'] = df_labels
    else:
        print('No cell label file is provided, \'unknown\' is used as the default cell labels')
        adata.obs['label'] = 'unknown'
    return None


def add_cell_colors(adata,file_path='',file_name=None):
    """Add cell colors.

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    file_path: `str`, optional (default: '')
        The file path of cell label color file.
    file_name: `str`, optional (default: None)
        The file name of cell label color file. If file_name is not specified, random color are generated for each cell label.
        
    Returns
    -------
    updates `adata` with the following fields.
    label_color: `pandas.core.series.Series` (`adata.obs['label_color']`,dtype `str`)
        Array of #observations that stores the color of each cell (hex color code).
    label_color: `dict` (`adata.uns['label_color']`,dtype `str`)
        Array of #observations that stores the color of each cell (hex color code).        
    """

    _fp = lambda f:  os.path.join(file_path,f)
    labels_unique = adata.obs['label'].unique()
    if(file_name!=None):
        df_colors = pd.read_csv(_fp(file_name),sep='\t',header=None,index_col=None,names=['label','color'],
                                dtype=str,compression= 'gzip' if file_name.split('.')[-1]=='gz' else None)
        df_colors['label'] = df_colors['label'].str.replace('/','-')   
        adata.uns['label_color'] = {df_colors.iloc[x,0]:df_colors.iloc[x,1] for x in range(df_colors.shape[0])}
    else:
        print('No cell color file is provided, random color is generated for each cell label')
        if(len(labels_unique)==1):
            adata.uns['label_color'] = {labels_unique[0]:'gray'}
        else:
            list_colors = sns.color_palette("hls",n_colors=len(labels_unique)).as_hex()
            adata.uns['label_color'] = {x:list_colors[i] for i,x in enumerate(labels_unique)}
    # df_cell_colors = adata.obs.copy()
    # df_cell_colors['label_color'] = ''
    # for x in labels_unique:
    #     id_cells = np.where(adata.obs['label']==x)[0]
    #     df_cell_colors.loc[df_cell_colors.index[id_cells],'label_color'] = adata.uns['label_color'][x]
    # adata.obs['label_color'] = df_cell_colors['label_color']
    return None

# def add_cell_colors(adata,file_path = None,file_name=None,key_label='label',key_color = 'label_color'):
#     labels_unique = adata.obs[key_label].unique()
#     if(file_name!=None):
#         df_colors = pd.read_csv(file_path+file_name,sep='\t',header=None,index_col=None,
#                                 compression= 'gzip' if file_name.split('.')[-1]=='gz' else None)

#         adata.uns[key_color] = {df_colors.iloc[x,0]:df_colors.iloc[x,1] for x in range(df_colors.shape[0])}
#     else:
#         list_colors = sns.color_palette("hls",n_colors=len(labels_unique)).as_hex()
#         adata.uns[key_color] = {x:list_colors[i] for i,x in enumerate(labels_unique)}
#     df_cell_colors = adata.obs.copy()
#     df_cell_colors[key_color] = ''
#     for x in labels_unique:
#         id_cells = np.where(adata.obs[key_label]==x)[0]
#         df_cell_colors.loc[df_cell_colors.index[id_cells],key_color] = adata.uns[key_color][x]
#     adata.obs[key_color] = df_cell_colors[key_color]
#     return None


def cal_qc(adata,expr_cutoff=1,assay='rna'):
    """Calculate quality control metrics.
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    expr_cutoff: `float`, optional (default: 1)
        Expression cutoff. If greater than expr_cutoff,the feature is considered 'expressed'
    assay: `str`, optional (default: 'rna')
            Choose from {{'rna','atac'}},case insensitive

    Returns
    -------
    updates `adata` with the following fields.
    n_counts: `pandas.Series` (`adata.var['n_counts']`,dtype `int`)
       The number of read count each gene has.
    n_cells: `pandas.Series` (`adata.var['n_cells']`,dtype `int`)
       The number of cells in which each gene is expressed.
    pct_cells: `pandas.Series` (`adata.var['pct_cells']`,dtype `float`)
       The percentage of cells in which each gene is expressed.
    n_counts: `pandas.Series` (`adata.obs['n_counts']`,dtype `int`)
       The number of read count each cell has.
    n_genes: `pandas.Series` (`adata.obs['n_genes']`,dtype `int`)
       The number of genes expressed in each cell.
    pct_genes: `pandas.Series` (`adata.obs['pct_genes']`,dtype `float`)
       The percentage of genes expressed in each cell.
    n_peaks: `pandas.Series` (`adata.obs['n_peaks']`,dtype `int`)
       The number of peaks expressed in each cell.
    pct_peaks: `pandas.Series` (`adata.obs['pct_peaks']`,dtype `int`)
       The percentage of peaks expressed in each cell.
    pct_mt: `pandas.Series` (`adata.obs['pct_mt']`,dtype `float`)
       the percentage of counts in mitochondrial genes
    """    
    
    assay = assay.lower()
    assert assay in ['rna','atac'], "`assay` must be chosen from ['rna','atac']"
    
    # if(not issparse(adata.X)):
    #     adata.X = csr_matrix(adata.X)

    n_counts = adata.X.sum(axis=0)
    adata.var['n_counts'] = n_counts
    n_cells = (adata.X>=expr_cutoff).sum(axis=0)
    adata.var['n_cells'] = n_cells 
    adata.var['pct_cells'] = n_cells/adata.shape[0]

    n_counts = adata.X.sum(axis=1)
    adata.obs['n_counts'] = n_counts
    n_features = (adata.X>=expr_cutoff).sum(axis=1)
    if(assay=='atac'):
        adata.obs['n_peaks'] = n_features  
        adata.obs['pct_peaks'] = n_features/adata.shape[1]
    if(assay=='rna'):
        adata.obs['n_genes'] = n_features
        adata.obs['pct_genes'] = n_features/adata.shape[1]
        r = re.compile("^MT-",flags=re.IGNORECASE)
        mt_genes = list(filter(r.match, adata.var_names))
        if(len(mt_genes)>0):
            n_counts_mt = adata[:,mt_genes].X.sum(axis=1)
            adata.obs['pct_mt'] = n_counts_mt/n_counts
        else:
            adata.obs['pct_mt'] = 0
    adata.uns['assay'] = assay

def plot_qc(adata,
            fig_size=(4,4),jitter=0.4,size=1,log_scale=None,hist_plot=None,
            pad=1.08,w_pad=None,h_pad=3,
            save_fig=False,fig_path=None,fig_name='qc.pdf',):
    
    """Plot QC metrics
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    fig_size: `tuple`, optional (default: (4,4))
        figure size.
    jitter: `float`, (default: 0.4)
        Amount of jitter (only along the categorical axis) to apply.
    size : `int`, (default: None)
        Size of the jitter points.
    log_scale: `list`, (default: None)
        Indices of plots to use log-transformation (0-indexed). 
        Index starts at 0 in the upper left corner and increases to the right.
    hist_plot: `list`, (default: None)
        Indices of plots to use histogram plot (0-indexed). 
        Index starts at 0 in the upper left corner and increases to the right.        
    pad: `float`, optional (default: 1.08)
        Padding between the figure edge and the edges of subplots, as a fraction of the font size.
    h_pad, w_pad: `float`, optional (default: None)
        Padding (height/width) between edges of adjacent subplots, as a fraction of the font size. Defaults to pad.
    save_fig: `bool`, optional (default: False)
        if True,save the figure.
    fig_path: `str`, optional (default: None)
        if save_fig is True, specify figure path. if None, adata.uns['workdir'] will be used.
    fig_name: `str`, optional (default: 'qc.pdf')
        if save_fig is True, specify figure name.

    Returns
    -------
    None

    """ 

    assert 'assay' in adata.uns_keys(),'please run `cal_qc(adata)` first'
    
    if(fig_path is None):
        fig_path = adata.uns['workdir']

    if(adata.uns['assay']=='rna'):
        feature = 'Genes'
        obs_list = ['n_counts','n_genes','pct_genes','pct_mt',] 
    if(adata.uns['assay']=='atac'):
        feature = 'Peaks'
        obs_list = ['n_counts','n_peaks','pct_peaks',]
    var_list = ['n_counts','n_cells','pct_cells']

    if(log_scale is None):
        log_scale = []
    for x in log_scale:
        if(x>=(len(obs_list)+len(var_list))):
            raise ValueError("Index '%s' out of range"  % (x))
            
    if(hist_plot is None):
        hist_plot = []
    for x in hist_plot:
        if(x>=(len(obs_list)+len(var_list))):
            raise ValueError("Index '%s' out of range"  % (x))            
  
    fig_ncol = len(obs_list)
    fig_nrow = 2    
    fig = plt.figure(figsize=(fig_size[0]*fig_ncol*1.2,fig_size[1]*fig_nrow*1.05))
    df_obs = adata.obs[obs_list].copy()
    df_var = adata.var[var_list].copy()
    
    for i in range(len(obs_list)+len(var_list)):
        ax_i = fig.add_subplot(fig_nrow,fig_ncol,i+1)
        if(i<fig_ncol):
            ann = obs_list[i]
            if(i in log_scale):
                df_obs['log10_'+ann] = np.log10(df_obs[ann]+1)
                ann = 'log10_'+ann 
            if(i in hist_plot):
                sc_i1 = sns.histplot(ax=ax_i, data=df_obs, y=ann)
                # sc_i1 = sns.distplot(df_obs[ann],vertical=True)
                ax_i.set_xticks([]) 
            else:
                sc_i1 = sns.violinplot(ax=ax_i,y=ann,data=df_obs,inner=None)
                sc_i2 = sns.stripplot(ax=ax_i,y=ann,data=df_obs,
                                     color='black',jitter=jitter,s=size) 
        else:
            ann = var_list[i%fig_ncol]
            if(i in log_scale):
                df_var['log10_'+ann] = np.log10(df_var[ann]+1)
                ann = 'log10_'+ann 
            if(i in hist_plot):
                sc_i1 = sns.histplot(ax=ax_i, data=df_var, y=ann)
                # sc_i1 = sns.distplot(df_var[ann],vertical=True)
                ax_i.set_xticks([])
            else:
                sc_i1 = sns.violinplot(ax=ax_i,y=ann,data=df_var,inner=None)
                sc_i2 = sns.stripplot(ax=ax_i,y=ann,data=df_var,
                                     color='black',jitter=jitter,s=size)  
        ax_i.set_title(ann)
        ax_i.locator_params(axis='y',nbins=6)
        ax_i.tick_params(axis="y",pad=-2)
        ax_i.spines['right'].set_visible(False)
        ax_i.spines['top'].set_visible(False)
        if(i==0):
            ax_i.set_ylabel('Cells',rotation=0, size='large',labelpad=50)
        elif(i==fig_ncol):
            ax_i.set_ylabel(feature,rotation=0, size='large',labelpad=50)
        else:
            ax_i.set_ylabel('')
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad)
    if(save_fig):
        plt.savefig(os.path.join(fig_path,fig_name),pad_inches=1,bbox_inches='tight')
        plt.close(fig)


def filter_features(adata,
                    min_n_cells = None, max_n_cells=None,
                    min_pct_cells = None, max_pct_cells=None,
                    min_n_counts = None, max_n_counts=None,
                    expr_cutoff = 1,assay=None):
    """Filter out features based on different metrics.

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    min_n_cells: `int`, optional (default: None)
        Minimum number of cells expressing one feature
    min_pct_cells: `float`, optional (default: None)
        Minimum percentage of cells expressing one feature
    min_n_counts: `int`, optional (default: None)
        Minimum number of read count for one feature
    expr_cutoff: `float`, optional (default: 1)
        Expression cutoff. If greater than expr_cutoff,the feature is considered 'expressed'
    assay: `str`, optional (default: 'rna')
            Choose from {{'rna','atac'}},case insensitive
            
    Returns
    -------
    updates `adata` with a subset of features that pass the filtering.   
    updates `adata` with the following fields if cal_qc() was not performed.
    n_counts: `pandas.Series` (`adata.var['n_counts']`,dtype `int`)
       The number of read count each gene has.
    n_cells: `pandas.Series` (`adata.var['n_cells']`,dtype `int`)
       The number of cells in which each gene is expressed.
    pct_cells: `pandas.Series` (`adata.var['pct_cells']`,dtype `float`)
       The percentage of cells in which each gene is expressed.
    """
    if(assay is None):
        if('assay' in adata.uns_keys()):
            assay = adata.uns['assay']
        else:
            raise Exception("Please either run 'st.cal_qc()' or specify the parameter`assay`")
    assay = assay.lower()
    assert assay in ['rna','atac'], "`assay` must be chosen from ['rna','atac']"
    
    if(assay=='rna'):
        feature = 'genes'
    if(assay=='atac'):
        feature = 'peaks'
    
    if('n_counts' in adata.var_keys()):
        n_counts = adata.var['n_counts']
    else:
        n_counts = np.sum(adata.X,axis=0).astype(int)
        adata.var['n_counts'] = n_counts
    if('n_cells' in adata.var_keys()):
        n_cells = adata.var['n_cells']
    else:
        n_cells = np.sum(adata.X>=expr_cutoff,axis=0).astype(int)
        adata.var['n_cells'] = n_cells   
    if('pct_cells' in adata.var_keys()): 
        pct_cells = adata.var['pct_cells']
    else:
        pct_cells = n_cells/adata.shape[0]
        adata.var['pct_cells'] = pct_cells
        
    if(sum(list(map(lambda x: x is None,[min_n_cells,min_pct_cells,min_n_counts,
                                         max_n_cells,max_pct_cells,max_n_counts,])))==6):
        print('No filtering')
    else:
        feature_subset = np.ones(len(adata.var_names),dtype=bool)
        if(min_n_cells!=None):
            print('Filter '+feature+' based on min_n_cells')
            feature_subset = (n_cells>=min_n_cells) & feature_subset
        if(max_n_cells!=None):
            print('Filter '+feature+' based on max_n_cells')
            feature_subset = (n_cells<=max_n_cells) & feature_subset
        if(min_pct_cells!=None):
            print('Filter '+feature+' based on min_pct_cells')
            feature_subset = (pct_cells>=min_pct_cells) & feature_subset
        if(max_pct_cells!=None):
            print('Filter '+feature+' based on max_pct_cells')
            feature_subset = (pct_cells<=max_pct_cells) & feature_subset
        if(min_n_counts!=None):
            print('Filter '+feature+' based on min_n_counts')
            feature_subset = (n_counts>=min_n_counts) & feature_subset 
        if(max_n_counts!=None):
            print('Filter '+feature+' based on max_n_counts')
            feature_subset = (n_counts<=max_n_counts) & feature_subset 
        adata._inplace_subset_var(feature_subset)
        print('After filtering out low-expressed '+feature+': ')
        print(str(adata.shape[0])+' cells, ' + str(adata.shape[1])+' '+feature)
    return None

def filter_cells(adata,
                 min_n_features = None, max_n_features = None,
                 min_pct_features = None, max_pct_features = None,
                 min_n_counts=None, max_n_counts = None,
                 expr_cutoff = 1,assay=None):
    """Filter out cells based on different metrics.

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    min_n_features: `int`, optional (default: None)
        Minimum number of features expressed
    min_pct_features: `float`, optional (default: None)
        Minimum percentage of features expressed
    min_n_counts: `int`, optional (default: None)
        Minimum number of read count for one cell
    expr_cutoff: `float`, optional (default: 1)
        Expression cutoff. If greater than expr_cutoff,the gene is considered 'expressed'
    assay: `str`, optional (default: 'rna')
            Choose from {{'rna','atac'}},case insensitive    
    Returns
    -------
    updates `adata` with a subset of cells that pass the filtering.      
    updates `adata` with the following fields if cal_qc() was not performed.
    n_counts: `pandas.Series` (`adata.obs['n_counts']`,dtype `int`)
       The number of read count each cell has.
    n_genes: `pandas.Series` (`adata.obs['n_genes']`,dtype `int`)
       The number of genes expressed in each cell.
    pct_genes: `pandas.Series` (`adata.obs['pct_genes']`,dtype `float`)
       The percentage of genes expressed in each cell.
    n_peaks: `pandas.Series` (`adata.obs['n_peaks']`,dtype `int`)
       The number of peaks expressed in each cell.
    pct_peaks: `pandas.Series` (`adata.obs['pct_peaks']`,dtype `int`)
       The percentage of peaks expressed in each cell.
    """
    if(assay is None):
        if('assay' in adata.uns_keys()):
            assay = adata.uns['assay']
        else:
            raise Exception("Please either run 'st.cal_qc()' or specify the parameter`assay`")
    assay = assay.lower()
    assert assay in ['rna','atac'], "`assay` must be chosen from ['rna','atac']"
            
    if('n_counts' in adata.obs_keys()):
        n_counts = adata.obs['n_counts']
    else:
        n_counts = np.sum(adata.X,axis=1).astype(int)
        adata.obs['n_counts'] = n_counts
    if(assay == 'rna'):
        if('n_genes' in adata.obs_keys()):
            n_features = adata.obs['n_genes']
        else:
            n_features = np.sum(adata.X>=expr_cutoff,axis=1).astype(int)
            adata.obs['n_genes'] = n_features       
        if('pct_genes' in adata.obs_keys()):
            pct_features = adata.obs['pct_genes']
        else:
            pct_features = n_features/adata.shape[1]
            adata.obs['pct_genes'] = pct_features
        feature = 'genes'
    if(assay == 'atac'):
        if('n_peaks' in adata.obs_keys()):
            n_features = adata.obs['n_peaks']
        else:
            n_features = np.sum(adata.X>=expr_cutoff,axis=1).astype(int)
            adata.obs['n_peaks'] = n_features       
        if('pct_peaks' in adata.obs_keys()):
            pct_features = adata.obs['pct_peaks']
        else:
            pct_features = n_features/adata.shape[1]
            adata.obs['pct_peaks'] = pct_features
        feature = 'peaks'            

    if(sum(list(map(lambda x: x is None,[min_n_features,min_pct_features,min_n_counts,
                                         max_n_features,max_pct_features,max_n_counts,])))==6):
        print('No filtering')    
    else:
        cell_subset = np.ones(len(adata.obs_names),dtype=bool)
        if(min_n_features!=None):
            print('filter cells based on min_n_features')
            cell_subset = (n_features>=min_n_features) & cell_subset
        if(max_n_features!=None):
            print('filter cells based on max_n_features')
            cell_subset = (n_features<=max_n_features) & cell_subset
        if(min_pct_features!=None):
            print('filter cells based on min_pct_features')
            cell_subset = (pct_features>=min_pct_features) & cell_subset
        if(max_pct_features!=None):
            print('filter cells based on max_pct_features')
            cell_subset = (pct_features<=max_pct_features) & cell_subset
        if(min_n_counts!=None):
            print('filter cells based on min_n_counts')
            cell_subset = (n_counts>=min_n_counts) & cell_subset 
        if(max_n_counts!=None):
            print('filter cells based on max_n_counts')
            cell_subset = (n_counts<=max_n_counts) & cell_subset 
        adata._inplace_subset_obs(cell_subset)
        print('after filtering out low-quality cells: ')
        print(str(adata.shape[0])+' cells, ' + str(adata.shape[1])+' '+feature)
    return None

def log_transform(adata,base=2):
    """Logarithmize gene expression.

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    base: `int`, optional (default: 2)
        The base used to calculate logarithm

    Returns
    -------
    updates `adata` with the following fields.
    X: `numpy.ndarray` (`adata.X`)
        Store #observations × #var_genes logarithmized data matrix.
    """

    adata.X = np.log2(adata.X+1)/np.log2(base)
    return None


# def normalize_per_cell(adata):
#     """Normalize gene expression based on library size.

#     Parameters
#     ----------
#     adata: AnnData
#         Annotated data matrix.

#     Returns
#     -------
#     updates `adata` with the following fields.
#     X: `numpy.ndarray` (`adata.X`)
#         Store #observations × #var_genes normalized data matrix.
#     """
#     adata.X = (np.divide(adata.X.T,adata.X.sum(axis=1)).T)*1e6

def normalize(adata,method='lib_size'):
    """Normalize count matrix.

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    method: `str`, optional (default: 'lib_size')
        Choose from {{'lib_size','tf_idf'}}
        Method used for dimension reduction.
        'lib_size': Total-count normalize (library-size correct)
        'tf_idf': TF-IDF (term frequency–inverse document frequency) transformation
    Returns
    -------
    updates `adata` with the following fields.
    X: `numpy.ndarray` (`adata.X`)
        Store #observations × #var_genes normalized data matrix.
    """
    if(method not in ['lib_size','tf_idf']):
        raise ValueError("unrecognized method '%s'" % method)
    if(method == 'lib_size'):
        adata.X = (np.divide(adata.X.T,adata.X.sum(axis=1)).T)*1e6
    if(method == 'tf_idf'):
        adata.X = cal_tf_idf(adata.X)

def remove_mt_genes(adata):
    """remove mitochondrial genes.

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.

    Returns
    -------
    updates `adata` with a subset of genes that excluded mitochondrial genes.
    """        

    r = re.compile("^MT-",flags=re.IGNORECASE)
    mt_genes = list(filter(r.match, adata.var_names))
    if(len(mt_genes)>0):
        print('remove mitochondrial genes:')
        print(mt_genes)
        gene_subset = ~adata.var_names.isin(mt_genes)
        adata._inplace_subset_var(gene_subset)


def select_variable_genes(adata,loess_frac=0.01,percentile=95,n_genes = None,n_jobs = multiprocessing.cpu_count(),
                          save_fig=False,fig_name='std_vs_means.pdf',fig_path=None,fig_size=(4,4),
                          pad=1.08,w_pad=None,h_pad=None):

    """Select the most variable genes.

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    loess_frac: `float`, optional (default: 0.1)
        Between 0 and 1. The fraction of the data used when estimating each y-value in LOWESS function.
    percentile: `int`, optional (default: 95)
        Between 0 and 100. Specify the percentile to select genes.Genes are ordered based on its distance from the fitted curve.
    n_genes: `int`, optional (default: None)
        Specify the number of selected genes. Genes are ordered based on its distance from the fitted curve.
    n_jobs: `int`, optional (default: all available cpus)
        The number of parallel jobs to run when calculating the distance from each gene to the fitted curve
    save_fig: `bool`, optional (default: False)
        if True,save the figure.
    fig_size: `tuple`, optional (default: (4,4))
        figure size.
    fig_path: `str`, optional (default: '')
        if empty, adata.uns['workdir'] will be used.
    fig_name: `str`, optional (default: 'std_vs_means.pdf')
        if save_fig is True, specify figure name.
    pad: `float`, optional (default: 1.08)
        Padding between the figure edge and the edges of subplots, as a fraction of the font size.
    h_pad, w_pad: `float`, optional (default: None)
        Padding (height/width) between edges of adjacent subplots, as a fraction of the font size. Defaults to pad.

    Returns
    -------
    updates `adata` with the following fields.
    var_genes: `numpy.ndarray` (`adata.obsm['var_genes']`)
        Store #observations × #var_genes data matrix used for subsequent dimension reduction.
    var_genes: `pandas.core.indexes.base.Index` (`adata.uns['var_genes']`)
        The selected variable gene names.
    """

    # if(not issparse(adata.X)):
    #     adata.X = csr_matrix(adata.X) 

    if(fig_path is None):
        fig_path = adata.uns['workdir']  
    fig_size = mpl.rcParams['figure.figsize'] if fig_size is None else fig_size
    mean_genes = np.mean(adata.X,axis=0)
    std_genes = np.std(adata.X,ddof=1,axis=0)
    loess_fitted = lowess(std_genes,mean_genes,return_sorted=False,frac=loess_frac)
    residuals = std_genes - loess_fitted
    XP = np.column_stack((np.sort(mean_genes),loess_fitted[np.argsort(mean_genes)]))
    mat_p = np.column_stack((mean_genes,std_genes))
    with multiprocessing.Pool(processes=n_jobs) as pool:
        dist_point_to_curve = pool.starmap(project_point_to_curve_distance,[(XP,mat_p[i,]) for i in range(XP.shape[0])])
    mat_sign = np.ones(XP.shape[0])
    mat_sign[np.where(residuals<0)[0]] = -1
    dist_point_to_curve = np.array(dist_point_to_curve)*mat_sign
    if(n_genes is None):
        cutoff = np.percentile(dist_point_to_curve,percentile)
        id_var_genes = np.where(dist_point_to_curve>cutoff)[0]
        id_non_var_genes = np.where(residuals<=cutoff)[0]
    else:
        id_var_genes = np.argsort(dist_point_to_curve)[::-1][:n_genes]
        id_non_var_genes = np.argsort(dist_point_to_curve)[::-1][n_genes:]
 
    adata.obsm['var_genes'] = adata.X[:,id_var_genes].copy()
    adata.uns['var_genes'] = adata.var_names[id_var_genes]
    print(str(len(id_var_genes))+' variable genes are selected')
    ###plotting
    fig = plt.figure(figsize=fig_size)      
    plt.scatter(mean_genes[id_non_var_genes], std_genes[id_non_var_genes],s=5,alpha=0.2,zorder=1,c='#6baed6')
    plt.scatter(mean_genes[id_var_genes], std_genes[id_var_genes],s=5,alpha=0.9,zorder=2,c='#EC4E4E')
    plt.plot(np.sort(mean_genes), loess_fitted[np.argsort(mean_genes)],linewidth=3,zorder=3,c='#3182bd')
    plt.xlabel('mean value')
    plt.ylabel('standard deviation')
    plt.locator_params(axis='x',nbins=5)
    plt.locator_params(axis='y',nbins=5)
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad)
    if(save_fig):
        plt.savefig(os.path.join(fig_path,fig_name),pad_inches=1,bbox_inches='tight')
        plt.close(fig)
    return None

def select_gini_genes(adata,loess_frac=0.1,percentile=95,n_genes = None,
                          save_fig=False,fig_name='gini_vs_max.pdf',fig_path=None,fig_size=(4,4),
                          pad=1.08,w_pad=None,h_pad=None):

    """Select high gini genes for rare cell types.
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    loess_frac: `float`, optional (default: 0.1)
        Between 0 and 1. The fraction of the data used when estimating each y-value in LOWESS function.
    percentile: `int`, optional (default: 95)
        Between 0 and 100. Specify the percentile to select genes.Genes are ordered based on the residuals.
    n_genes: `int`, optional (default: None)
        Specify the number of selected genes. Genes are ordered based on the residuals.
    save_fig: `bool`, optional (default: False)
        if True,save the figure.
    fig_size: `tuple`, optional (default: (4,4))
        figure size.
    fig_path: `str`, optional (default: '')
        if empty, adata.uns['workdir'] will be used.
    fig_name: `str`, optional (default: 'gini_vs_max.pdf')
        if save_fig is True, specify figure name.
    pad: `float`, optional (default: 1.08)
        Padding between the figure edge and the edges of subplots, as a fraction of the font size.
    h_pad, w_pad: `float`, optional (default: None)
        Padding (height/width) between edges of adjacent subplots, as a fraction of the font size. Defaults to pad.

    Returns
    -------
    updates `adata` with the following fields.
    gini: `pandas.core.series.Series` (`adata.var['gini']`,dtype `float`)
        Gini coefficients for all genes.
    gini_genes: `pandas.core.indexes.base.Index` (`adata.uns['gini_genes']`)
        The selected high gini genes.
    """

    if(fig_path is None):
        fig_path = adata.uns['workdir']
    fig_size = mpl.rcParams['figure.figsize'] if fig_size is None else fig_size
    if('gini' not in adata.var.columns):
        gini_values = np.array([gini(adata[:,x].X) for x in adata.var_names])
        adata.var['gini'] = gini_values
    gini_values = adata.var['gini']
    max_genes = np.max(adata.X,axis=0)
    loess_fitted = lowess(gini_values,max_genes,return_sorted=False,frac=loess_frac)
    residuals = gini_values - loess_fitted
    if(n_genes is None):
        cutoff = np.percentile(residuals,percentile)
        id_gini_genes = np.where(residuals>cutoff)[0]
    else:
        id_gini_genes = np.argsort(residuals)[::-1][:n_genes] 
    
    adata.uns['gini_genes'] = adata.var_names[id_gini_genes]
    print(str(len(id_gini_genes))+' gini genes are selected')
    ###plotting
    fig = plt.figure(figsize=fig_size)      
    plt.scatter(max_genes, gini_values,s=5,alpha=0.2,zorder=1,c='#6baed6')
    plt.scatter(max_genes[id_gini_genes], gini_values[id_gini_genes],s=5,alpha=0.9,zorder=2,c='#EC4E4E')
    plt.plot(np.sort(max_genes), loess_fitted[np.argsort(max_genes)],linewidth=3,zorder=3,c='#3182bd')
    plt.xlabel('max gene expression')
    plt.ylabel('Gini coefficient')
    plt.locator_params(axis='x',nbins=5)
    plt.locator_params(axis='y',nbins=5)
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad)
    if(save_fig):
        plt.savefig(os.path.join(fig_path,fig_name),pad_inches=1,bbox_inches='tight')
        plt.close(fig)
    return None


def select_top_principal_components(adata,feature=None,n_pc = 15,max_pc = 100,first_pc = False,use_precomputed=True,
                                    save_fig=False,fig_name='top_pcs.pdf',fig_path=None,fig_size=(4,4),
                                    pad=1.08,w_pad=None,h_pad=None):
    """Select top principal components.
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    feature: `str`, optional (default: None)
        Choose from {{'var_genes'}}
        Features used for pricipal component analysis
        If None, all the genes will be used.
        IF 'var_genes', the most variable genes obtained from select_variable_genes() will be used.
    n_pc: `int`, optional (default: 15)
        The number of selected principal components.
    max_pc: `int`, optional (default: 100)
        The maximum number of principal components used for variance Ratio plot.
    first_pc: `bool`, optional (default: False)
        If True, the first principal component will be included
    use_precomputed: `bool`, optional (default: True)
        If True, the PCA results from previous computing will be used
    save_fig: `bool`, optional (default: False)
        if True,save the figure.
    fig_size: `tuple`, optional (default: (4,4))
        figure size.
    fig_path: `str`, optional (default: None)
        if None, adata.uns['workdir'] will be used.
    fig_name: `str`, optional (default: 'top_pcs.pdf')
        if save_fig is True, specify figure name.
    pad: `float`, optional (default: 1.08)
        Padding between the figure edge and the edges of subplots, as a fraction of the font size.
    h_pad, w_pad: `float`, optional (default: None)
        Padding (height/width) between edges of adjacent subplots, as a fraction of the font size. Defaults to pad.

    Returns
    -------
    updates `adata` with the following fields.
    pca: `numpy.ndarray` (`adata.obsm['pca']`)
        Store #observations × n_components data matrix after pca. Number of components to keep is min(#observations,#variables)
    top_pcs: `numpy.ndarray` (`adata.obsm['top_pcs']`)
        Store #observations × n_pc data matrix used for subsequent dimension reduction.
    top_pcs: `sklearn.decomposition.PCA` (`adata.uns['top_pcs']`)
        Store pca object.
    pca_variance_ratio: `numpy.ndarray` (`adata.uns['pca_variance_ratio']`)
        Percentage of variance explained by each of the selected components.
    """
    
    if(fig_path is None):
        fig_path = adata.uns['workdir']
    fig_size = mpl.rcParams['figure.figsize'] if fig_size is None else fig_size  
    if(feature is None):
        feature = 'all'
    if(feature not in ['all','var_genes']):
        raise ValueError("unrecognized feature '%s'" % feature)

    if(use_precomputed and ('pca' in adata.obsm_keys())):
        print('Importing precomputed principal components')
        X_pca = adata.obsm['pca']
        pca_variance_ratio = adata.uns['pca_variance_ratio']
    else:
        sklearn_pca = sklearnPCA(n_components=max_pc,svd_solver='arpack',random_state=42)
        if(feature == 'var_genes'):
            print('using top variable genes ...')
            trans = sklearn_pca.fit(adata.obsm['var_genes'])
            X_pca = trans.transform(adata.obsm['var_genes'])
            #X_pca = sklearn_pca.fit_transform(adata.obsm['var_genes'])
            pca_variance_ratio = trans.explained_variance_ratio_
            adata.obsm['pca'] = X_pca
            adata.uns['pca_variance_ratio'] = pca_variance_ratio                
        else:
            print('using all the features ...')
            trans = sklearn_pca.fit(adata.X)
            X_pca = trans.transform(adata.X)            
#             X_pca = sklearn_pca.fit_transform(adata.X)
            pca_variance_ratio = trans.explained_variance_ratio_
            adata.obsm['pca'] = X_pca
            adata.uns['pca_variance_ratio'] = pca_variance_ratio  
        adata.uns['top_pcs'] = trans
    if(first_pc):
        adata.obsm['top_pcs'] = X_pca[:,0:(n_pc)]
    else:
        #discard the first Principal Component
        adata.obsm['top_pcs'] = X_pca[:,1:(n_pc+1)]
    print(str(n_pc) + ' PCs are selected')

    if('params' not in adata.uns_keys()):
        adata.uns['params'] = dict()
    adata.uns['params']['select_top_principal_components'] = {'feature':feature,'first_pc':first_pc,'n_pc':n_pc}

    ##plotting
    fig = plt.figure(figsize=fig_size)
    plt.plot(range(max_pc),pca_variance_ratio[:max_pc])
    if(first_pc):
        plt.axvline(n_pc,c='red',ls = '--')
    else:
        plt.axvline(1,c='red',ls = '--')
        plt.axvline(n_pc+1,c='red',ls = '--')
    plt.xlabel('Principal Component')
    plt.ylabel('Variance Ratio')
    plt.locator_params(axis='x',nbins=5)
    plt.locator_params(axis='y',nbins=5)
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad)
    if(save_fig):
        plt.savefig(os.path.join(fig_path,fig_name),pad_inches=1,bbox_inches='tight')
        plt.close(fig)
    return None

def dimension_reduction(adata,n_neighbors=50, nb_pct = None,n_components = 3,n_jobs = 1,
                        feature='var_genes',method = 'se',eigen_solver=None):

    """Perform dimension reduction.
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    n_neighbors: `int`, optional (default: 50)
        The number of neighbor cells used for manifold learning (only valid when 'mlle','se', or 'umap' is specified).
    nb_pct: `float`, optional (default: None)
        The percentage of neighbor cells (when sepcified, it will overwrite n_neighbors).
    n_components: `int`, optional (default: 3)
        Number of components to keep.
    n_jobs: `int`, optional (default: 1)
        The number of parallel jobs to run.
    feature: `str`, optional (default: 'var_genes')
        Choose from {{'var_genes','top_pcs','all'}}
        Feature used for dimension reduction.
        'var_genes': most variable genes
        'top_pcs': top principal components
        'all': all available features (genes)
    method: `str`, optional (default: 'se')
        Choose from {{'se','mlle','umap','pca'}}
        Method used for dimension reduction.
        'se': Spectral embedding algorithm
        'mlle': Modified locally linear embedding algorithm
        'umap': Uniform Manifold Approximation and Projection
        'pca': Principal component analysis
    eigen_solver: `str`, optional (default: None)
        For 'mlle', choose from {{'arpack', 'dense'}}
        For 'se', choose from {{'arpack', 'lobpcg', or 'amg'}}
        The eigenvalue decomposition strategy to use
   
    Returns
    -------
    updates `adata` with the following fields.
    
    X_dr : `numpy.ndarray` (`adata.obsm['X_dr']`)
        A #observations × n_components data matrix after dimension reduction.
    X_mlle : `numpy.ndarray` (`adata.obsm['X_mlle']`)
        Store #observations × n_components data matrix after mlle.
    X_se : `numpy.ndarray` (`adata.obsm['X_se']`)
        Store #observations × n_components data matrix after spectral embedding.    
    X_umap : `numpy.ndarray` (`adata.obsm['X_umap']`)
        Store #observations × n_components data matrix after umap.
    X_pca : `numpy.ndarray` (`adata.obsm['X_pca']`)
        Store #observations × n_components data matrix after pca.
    trans_mlle : `sklearn.manifold.locally_linear.LocallyLinearEmbedding` (`adata.uns['trans_mlle']`)
        Store mlle object
    trans_se : `sklearn.manifold.spectral_embedding_.SpectralEmbedding` (`adata.uns['trans_se']`)
        Store se object
    trans_umap : `umap.UMAP` (`adata.uns['trans_umap']`)
        Store umap object
    trans_pca : `sklearn.decomposition.PCA` (`adata.uns['trans_pca']`)
        Store pca object 
    """

    if(feature not in ['var_genes','top_pcs','all']):
        raise ValueError("unrecognized feature '%s'" % feature)
    if(method not in ['mlle','se','umap','pca']):
        raise ValueError("unrecognized method '%s'" % method)
    if(feature == 'var_genes'):
        input_data = adata.obsm['var_genes']
    if(feature == 'top_pcs'):
        input_data = adata.obsm['top_pcs']
    if(feature == 'all'):
        input_data = adata.X
    print('feature ' + feature + ' is being used ...')
    print(str(n_jobs)+' cpus are being used ...')
    if(nb_pct!=None):
        n_neighbors = int(np.around(input_data.shape[0]*nb_pct))

    if(method == 'mlle'):
        np.random.seed(2)
        if(eigen_solver==None):
            if(input_data.shape[0]<=2000):
                reducer = LocallyLinearEmbedding(n_neighbors=n_neighbors, 
                                                     n_components=n_components,
                                                     n_jobs = n_jobs,
                                                     method = 'modified',eigen_solver = 'dense',random_state=10,
                                                     neighbors_algorithm = 'kd_tree')
            else:
                reducer = LocallyLinearEmbedding(n_neighbors=n_neighbors, 
                                                     n_components=n_components,
                                                     n_jobs = n_jobs,
                                                     method = 'modified',eigen_solver = 'arpack',random_state=10,
                                                     neighbors_algorithm = 'kd_tree')
                
        else:
            reducer = LocallyLinearEmbedding(n_neighbors=n_neighbors, 
                                                 n_components=n_components,
                                                 n_jobs = n_jobs,
                                                 method = 'modified',eigen_solver = eigen_solver,random_state=10,
                                                 neighbors_algorithm = 'kd_tree')        
        trans = reducer.fit(input_data)
        adata.uns['trans_mlle'] = trans
        adata.obsm['X_mlle'] = trans.embedding_
        adata.obsm['X_dr'] = trans.embedding_
    if(method == 'se'):
        np.random.seed(2)
        reducer = SpectralEmbedding(n_neighbors=n_neighbors, 
                                         n_components=n_components,
                                         n_jobs = n_jobs,
                                         eigen_solver = eigen_solver,random_state=10)
        trans = reducer.fit(input_data)
        adata.uns['trans_se'] = trans
        adata.obsm['X_se'] = trans.embedding_
        adata.obsm['X_dr'] = trans.embedding_
    if(method == 'umap'):
        reducer = umap.UMAP(n_neighbors=n_neighbors,n_components=n_components,random_state=42)
        trans = reducer.fit(input_data)
        adata.uns['trans_umap'] = trans
        adata.obsm['X_umap'] = trans.embedding_
        adata.obsm['X_dr'] = trans.embedding_
    if(method == 'pca'):
        reducer = sklearnPCA(n_components=n_components,svd_solver='arpack',random_state=42)
        trans = reducer.fit(input_data)
        adata.uns['trans_pca'] = trans
        adata.obsm['X_pca'] = trans.transform(input_data) 
        adata.obsm['X_dr'] = adata.obsm['X_pca']
    if('params' not in adata.uns_keys()):
        adata.uns['params'] = dict()
    adata.uns['params']['dimension_reduction'] = {'feature':feature,'method':method,'n_components':n_components}
    return None


def plot_dimension_reduction(adata, n_components=None, comp1=0, comp2=1, comp3=2, color=None, key_graph='epg',
                             fig_size=None, fig_ncol=3, fig_legend_ncol=1, fig_legend_order=None,
                             vmin=None, vmax=None, alpha=0.8,
                             pad=1.08, w_pad=None, h_pad=None,
                             show_text=False, show_graph=False,
                             save_fig=False, fig_path=None, fig_name='dimension_reduction.pdf',
                             plotly=False, gif_name='tmp.gif'):
    """Plot the manifold where the graph is learned

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    n_components: `int`, optional (default: 3)
        Number of components to be plotted.
    comp1: `int`, optional (default: 0)
        Component used for x axis.
    comp2: `int`, optional (default: 1)
        Component used for y axis.
    comp3: `int`, optional (default: 2)
        Component used for z axis. (only valid when n_components=3)
    color: `list` optional (default: None)
        Column names of observations (adata.obs.columns) or variable names(adata.var_names). A list of names to be plotted.
    key_graph: `str`, optional (default: None):
        Choose from {{'epg','seed_epg','ori_epg'}}
        Specify gragh to be plotted.
        'epg' current elastic principal graph
        'seed_epg' seed structure used for elastic principal graph learning, which is obtained by running seed_elastic_principal_graph()
        'ori_epg' original elastic principal graph, which is obtained by running elastic_principal_graph()
    fig_size: `tuple`, optional (default: None)
        figure size.
    fig_ncol: `int`, optional (default: 3)
        the number of columns of the figure panel
    fig_legend_order: `dict`,optional (default: None)
        Specified order for the appearance of the annotation keys.Only valid for ategorical variable
        e.g. fig_legend_order = {'ann1':['a','b','c'],'ann2':['aa','bb','cc']}
    fig_legend_ncol: `int`, optional (default: 1)
        The number of columns that the legend has.
    vmin,vmax: `float`, optional (default: None)
        The min and max values are used to normalize continuous values. If None, the respective min and max of continuous values is used.
    alpha: `float`, optional (default: 0.8)
        0.0 transparent through 1.0 opaque
    pad: `float`, optional (default: 1.08)
        Padding between the figure edge and the edges of subplots, as a fraction of the font size.
    h_pad, w_pad: `float`, optional (default: None)
        Padding (height/width) between edges of adjacent subplots, as a fraction of the font size. Defaults to pad.
    show_text: `bool`, optional (default: False)
        If True, node state label will be shown
    show_graph: `bool`, optional (default: False)
        If True, the learnt principal graph will be shown
    save_fig: `bool`, optional (default: False)
        if True,save the figure.
    fig_path: `str`, optional (default: None)
        if save_fig is True, specify figure path. if None, adata.uns['workdir'] will be used.
    fig_name: `str`, optional (default: 'dimension_reduction.pdf')
        if save_fig is True, specify figure name.
    plotly: `bool`, optional (default: False)
        if True, plotly will be used to make interactive plots
    Returns
    -------
    None

    """

    if (fig_path is None):
        fig_path = adata.uns['workdir']
    fig_size = mpl.rcParams['figure.figsize'] if fig_size is None else fig_size

    if (n_components == None):
        n_components = min(3, adata.obsm['X_dr'].shape[1])
    if n_components not in [2, 3]:
        raise ValueError("n_components should be 2 or 3")
    if (n_components > adata.obsm['X_dr'].shape[1]):
        print(
            "`n_components` is greater than the available dimension. It is corrected to %s" % adata.obsm['X_dr'].shape[
                1])
        n_components = adata.obsm['X_dr'].shape[1]

    if (color is None):
        color = ['label']
    ###remove duplicate keys
    color = list(dict.fromkeys(color))
    dict_ann = dict()
    for ann in color:
        if (ann in adata.obs.columns):
            dict_ann[ann] = adata.obs[ann]
        elif (ann in adata.var_names):
            dict_ann[ann] = adata.obs_vector(ann)
        else:
            raise ValueError("could not find '%s' in `adata.obs.columns` and `adata.var_names`" % (ann))

    df_plot = pd.DataFrame(index=adata.obs.index, data=adata.obsm['X_dr'],
                           columns=['Dim' + str(x + 1) for x in range(adata.obsm['X_dr'].shape[1])])
    for ann in color:
        df_plot[ann] = dict_ann[ann]
    df_plot_shuf = df_plot.sample(frac=1, random_state=100)

    legend_order = {ann: np.unique(df_plot_shuf[ann]) for ann in color if is_string_dtype(df_plot_shuf[ann])}
    if (fig_legend_order is not None):
        if (not isinstance(fig_legend_order, dict)):
            raise TypeError("`fig_legend_order` must be a dictionary")
        for ann in fig_legend_order.keys():
            if (ann in legend_order.keys()):
                legend_order[ann] = fig_legend_order[ann]
            else:
                print("'%s' is ignored for ordering legend labels due to incorrect name or data type" % ann)

    if (show_graph or show_text):
        assert (all(np.isin(['epg', 'flat_tree'], adata.uns_keys()))), '''graph is not learnt yet. 
        please first run: `st.seed_elastic_principal_graph` and `st.elastic_principal_graph` to learn graph'''
        assert (key_graph in ['epg', 'seed_epg', 'ori_epg']), "key_graph must be one of ['epg','seed_epg','ori_epg']"
        if (fig_path is None):
            fig_path = adata.uns['workdir']
        if (key_graph == 'epg'):
            epg = adata.uns['epg']
            flat_tree = adata.uns['flat_tree']
        else:
            epg = adata.uns[key_graph.split('_')[0] + '_epg']
            flat_tree = adata.uns[key_graph.split('_')[0] + '_flat_tree']
        ft_node_pos = nx.get_node_attributes(flat_tree, 'pos')
        ft_node_label = nx.get_node_attributes(flat_tree, 'label')
        epg_node_pos = nx.get_node_attributes(epg, 'pos')

    if (plotly):
        for ann in color:
            if (n_components == 3):
                fig = px.scatter_3d(df_plot_shuf, x='Dim' + str(comp1 + 1), y='Dim' + str(comp2 + 1),
                                    z='Dim' + str(comp3 + 1), color=ann,
                                    opacity=alpha,
                                    color_continuous_scale=px.colors.sequential.Viridis,
                                    color_discrete_map=adata.uns[
                                        ann + '_color'] if ann + '_color' in adata.uns_keys() else {})
                fig.update_traces(marker=dict(size=2))
                if (show_graph):
                    for edge_i in flat_tree.edges():
                        branch_i_pos = np.array([epg_node_pos[i] for i in flat_tree.edges[edge_i]['nodes']])
                        edge_i_label = flat_tree.nodes[edge_i[0]]['label'] + '_' + flat_tree.nodes[edge_i[1]]['label']
                        curve_i = pd.DataFrame(branch_i_pos, columns=np.arange(branch_i_pos.shape[1]))
                        fig.add_trace(go.Scatter3d(x=curve_i[comp1],
                                                   y=curve_i[comp2],
                                                   z=curve_i[comp3],
                                                   mode='lines',
                                                   line=dict(color='black', width=3),
                                                   name=edge_i_label,
                                                   showlegend=True if is_string_dtype(df_plot[ann]) else False))
                if (show_text):
                    fig.add_trace(go.Scatter3d(x=np.array(list(ft_node_pos.values()))[:, comp1],
                                               y=np.array(list(ft_node_pos.values()))[:, comp2],
                                               z=np.array(list(ft_node_pos.values()))[:, comp3],
                                               mode='markers+text',
                                               opacity=1,
                                               marker=dict(size=4, color='#767070'),
                                               text=[ft_node_label[x] for x in ft_node_pos.keys()],
                                               textposition="bottom center",
                                               name='states',
                                               showlegend=True if is_string_dtype(df_plot[ann]) else False))
                fig.update_layout(legend={'itemsizing': 'constant'}, width=500, height=500,
                                  scene=dict(aspectmode='cube'))
            else:
                fig = px.scatter(df_plot_shuf, x='Dim' + str(comp1 + 1), y='Dim' + str(comp2 + 1), color=ann,
                                 opacity=alpha,
                                 color_continuous_scale=px.colors.sequential.Viridis,
                                 color_discrete_map=adata.uns[
                                     ann + '_color'] if ann + '_color' in adata.uns_keys() else {})
                if (show_graph):
                    for edge_i in flat_tree.edges():
                        branch_i_pos = np.array([epg_node_pos[i] for i in flat_tree.edges[edge_i]['nodes']])[:,
                                       [comp1, comp2]]
                        edge_i_label = flat_tree.nodes[edge_i[0]]['label'] + '_' + flat_tree.nodes[edge_i[1]]['label']
                        curve_i = pd.DataFrame(branch_i_pos, columns=np.arange(branch_i_pos.shape[1]))
                        fig.add_trace(go.Scatter(x=curve_i[comp1],
                                                 y=curve_i[comp2],
                                                 mode='lines',
                                                 line=dict(color='black', width=3),
                                                 name=edge_i_label,
                                                 showlegend=True if is_string_dtype(df_plot[ann]) else False))
                if (show_text):
                    fig.add_trace(go.Scatter(x=np.array(list(ft_node_pos.values()))[:, comp1],
                                             y=np.array(list(ft_node_pos.values()))[:, comp2],
                                             mode='markers+text',
                                             opacity=1,
                                             marker=dict(size=1.5 * mpl.rcParams['lines.markersize'], color='#767070'),
                                             text=[ft_node_label[x] for x in ft_node_pos.keys()],
                                             textposition="bottom center",
                                             name='states',
                                             showlegend=True if is_string_dtype(df_plot[ann]) else False))

                fig.update_layout(legend={'itemsizing': 'constant'}, width=500, height=500)
            fig.show(renderer="notebook")

    else:
        if (len(color) < fig_ncol):
            fig_ncol = len(color)
        fig_nrow = int(np.ceil(len(color) / fig_ncol))
        fig = plt.figure(figsize=(fig_size[0] * fig_ncol * 1.05, fig_size[1] * fig_nrow))
        axList = []
        from matplotlib import animation
        def init():
            for i, ann in enumerate(color):
                if (n_components == 3):
                    if (is_string_dtype(df_plot[ann])):
                        ax_i = fig.add_subplot(fig_nrow, fig_ncol, i + 1)
                        sc_i = sns.scatterplot(ax=ax_i,
                                               x="Dim1", y="Dim2",
                                               hue=ann, hue_order=legend_order[ann],
                                               data=df_plot_shuf,
                                               alpha=alpha, linewidth=0,
                                               palette=adata.uns[ann + '_color'] \
                                                   if (ann + '_color' in adata.uns_keys()) and (
                                                           set(adata.uns[ann + '_color'].keys()) >= set(
                                                       np.unique(df_plot_shuf[ann]))) \
                                                   else None
                                               )
                        colors_sns = sc_i.get_children()[0].get_facecolors()
                        if (ann + '_color' not in adata.uns_keys()):
                            colors_sns_scaled = (255 * colors_sns).astype(int)
                            adata.uns[ann + '_color'] = {df_plot_shuf[ann][i]: '#%02x%02x%02x' % (
                            colors_sns_scaled[i][0], colors_sns_scaled[i][1], colors_sns_scaled[i][2])
                                                         for i in np.unique(df_plot_shuf[ann], return_index=True)[1]}
                        legend_sns, labels_sns = ax_i.get_legend_handles_labels()
                        ax_i.remove()
                        ax_i = fig.add_subplot(fig_nrow, fig_ncol, i + 1, projection='3d')
                        axList.append(ax_i)
                        ax_i.scatter3D(df_plot_shuf['Dim' + str(comp1 + 1)], df_plot_shuf['Dim' + str(comp2 + 1)],
                                       df_plot_shuf['Dim' + str(comp3 + 1)],
                                       c=colors_sns, alpha=alpha, linewidth=0,
                                       zorder=-1, s=3)
                        ax_i.legend(legend_sns, labels_sns, loc='upper center', ncol=3, fontsize=8, markerscale=0.5)
                    else:
                        ax_i = fig.add_subplot(fig_nrow, fig_ncol, i + 1, projection='3d')
                        vmin_i = df_plot[ann].min() if vmin is None else vmin
                        vmax_i = df_plot[ann].max() if vmax is None else vmax
                        sc_i = ax_i.scatter3D(df_plot_shuf['Dim' + str(comp1 + 1)],
                                              df_plot_shuf['Dim' + str(comp2 + 1)],
                                              df_plot_shuf['Dim' + str(comp3 + 1)],
                                              # c=df_plot_shuf[ann],vmin=vmin_i,vmax=vmax_i,alpha=alpha,linewidth=0,
                                              c='black', vmin=vmin_i, vmax=vmax_i, alpha=alpha, linewidth=0,
                                              zorder=-1)

                        # cbar = plt.colorbar(sc_i, ax=ax_i, pad=0.1, fraction=0.05, aspect=30)
                        # cbar.solids.set_edgecolor("face")
                        # cbar.ax.locator_params(nbins=5)
                    if (show_graph):
                        for edge_i in flat_tree.edges():
                            branch_i_pos = np.array([epg_node_pos[i] for i in flat_tree.edges[edge_i]['nodes']])
                            curve_i = pd.DataFrame(branch_i_pos, columns=np.arange(branch_i_pos.shape[1]))
                            ax_i.plot(curve_i[comp1], curve_i[comp2], curve_i[comp3], c='black', zorder=5)
                    if (show_text):
                        for node_i in flat_tree.nodes():
                            ### matplotlib has zcorder issues for 3D plot.To avoid the incorrect order, here the node point is removed for 3D
                            # ax_i.scatter(ft_node_pos[node_i][comp1],ft_node_pos[node_i][comp2],ft_node_pos[node_i][comp3],
                            #              color='#767070',s=1.5*(mpl.rcParams['lines.markersize']**2),zorder=10)
                            ax_i.text(ft_node_pos[node_i][comp1], ft_node_pos[node_i][comp2],
                                      ft_node_pos[node_i][comp3],
                                      ft_node_label[node_i],
                                      color='red', fontsize=0.4 * mpl.rcParams['font.size'],
                                      ha='left', va='bottom',
                                      zorder=10)
                    ax_i.set_xlabel('PCoA' + str(comp1 + 1))
                    ax_i.set_ylabel('PCoA' + str(comp2 + 1))
                    ax_i.set_zlabel('PCoA' + str(comp3 + 1))
                    ax_i.set_xticklabels([])
                    ax_i.set_yticklabels([])
                    ax_i.set_zticklabels([])
                    # ax_i.legend(handles=legendElementList, loc='upper center', ncol=3, fontsize=20)
                else:
                    ax_i = fig.add_subplot(fig_nrow, fig_ncol, i + 1)
                    if (is_string_dtype(df_plot[ann])):
                        sc_i = sns.scatterplot(ax=ax_i,
                                               x='Dim' + str(comp1 + 1), y='Dim' + str(comp2 + 1),
                                               hue=ann, hue_order=legend_order[ann],
                                               data=df_plot_shuf,
                                               alpha=alpha, linewidth=0,
                                               palette=adata.uns[ann + '_color'] \
                                                   if (ann + '_color' in adata.uns_keys()) and (
                                                           set(adata.uns[ann + '_color'].keys()) >= set(
                                                       np.unique(df_plot_shuf[ann]))) \
                                                   else None
                                               )
                        legend_handles, legend_labels = ax_i.get_legend_handles_labels()
                        ax_i.legend(handles=legend_handles, labels=legend_labels,
                                    bbox_to_anchor=(1, 0.5), loc='center left', ncol=fig_legend_ncol,
                                    frameon=False,
                                    )
                        if (ann + '_color' not in adata.uns_keys()):
                            colors_sns = sc_i.get_children()[0].get_facecolors()
                            colors_sns_scaled = (255 * colors_sns).astype(int)
                            adata.uns[ann + '_color'] = {df_plot_shuf[ann][i]: '#%02x%02x%02x' % (
                            colors_sns_scaled[i][0], colors_sns_scaled[i][1], colors_sns_scaled[i][2])
                                                         for i in np.unique(df_plot_shuf[ann], return_index=True)[1]}
                        ### remove legend title
                        # ax_i.get_legend().texts[0].set_text("")
                    else:
                        vmin_i = df_plot[ann].min() if vmin is None else vmin
                        vmax_i = df_plot[ann].max() if vmax is None else vmax
                        sc_i = ax_i.scatter(df_plot_shuf['Dim' + str(comp1 + 1)], df_plot_shuf['Dim' + str(comp2 + 1)],
                                            c=df_plot_shuf[ann], vmin=vmin_i, vmax=vmax_i, alpha=alpha)
                        cbar = plt.colorbar(sc_i, ax=ax_i, pad=0.01, fraction=0.05, aspect=40)
                        cbar.solids.set_edgecolor("face")
                        cbar.ax.locator_params(nbins=5)
                    if (show_graph):
                        for edge_i in flat_tree.edges():
                            branch_i_pos = np.array([epg_node_pos[i] for i in flat_tree.edges[edge_i]['nodes']])
                            curve_i = pd.DataFrame(branch_i_pos, columns=np.arange(branch_i_pos.shape[1]))
                            ax_i.plot(curve_i[comp1], curve_i[comp2], c='black')
                    if (show_text):
                        for node_i in flat_tree.nodes():
                            ax_i.scatter(ft_node_pos[node_i][comp1], ft_node_pos[node_i][comp2],
                                         color='#767070', s=1.5 * (mpl.rcParams['lines.markersize'] ** 2))
                            ax_i.text(ft_node_pos[node_i][comp1], ft_node_pos[node_i][comp2], ft_node_label[node_i],
                                      color='black', fontsize=0.9 * mpl.rcParams['font.size'],
                                      ha='left', va='bottom')
                    ax_i.set_xlabel('PCoA' + str(comp1 + 1))
                    ax_i.set_ylabel('PCoA' + str(comp2 + 1))
                    ax_i.set_zlabel('PCoA' + str(comp3 + 1))
                    ax_i.set_title(ann)
                plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad)

        def rotate(angle):
            for ax in axList:
                ax.view_init(elev=10, azim=angle)
        rot_animation = animation.FuncAnimation(fig, rotate, init_func=init, frames=np.arange(0, 362, 2), interval=50)
        rot_animation.save(gif_name, dpi=200, writer='imagemagick')
        if (save_fig):
            plt.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')
            plt.close(fig)


def plot_branches(adata,n_components = None,comp1=0,comp2=1,comp3=2,key_graph='epg',
                  fig_size=None,
                  pad=1.08,w_pad=None,h_pad=None,
                  show_text=False,
                  save_fig=False,fig_path=None,fig_name='branches.pdf',
                  plotly=False):    
    """Plot branches. The branches contain all the nodes learnt from ElPiGraph
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    n_components: `int`, optional (default: 3)
        Number of components to be plotted.
    comp1: `int`, optional (default: 0)
        Component used for x axis.
    comp2: `int`, optional (default: 1)
        Component used for y axis.
    comp3: `int`, optional (default: 2)
        Component used for z axis. (only valid when n_components=3)
    key_graph: `str`, optional (default: None): 
        Choose from {{'epg','seed_epg','ori_epg'}}
        Specify gragh to be plotted.
        'epg' current elastic principal graph
        'seed_epg' seed structure used for elastic principal graph learning, which is obtained by running seed_elastic_principal_graph()
        'ori_epg' original elastic principal graph, which is obtained by running elastic_principal_graph()
    fig_size: `tuple`, optional (default: None)
        figure size.
    pad: `float`, optional (default: 1.08)
        Padding between the figure edge and the edges of subplots, as a fraction of the font size.
    h_pad, w_pad: `float`, optional (default: None)
        Padding (height/width) between edges of adjacent subplots, as a fraction of the font size. Defaults to pad.
    show_text: `bool`, optional (default: False)
        If True, node state label will be shown
    save_fig: `bool`, optional (default: False)
        if True,save the figure.
    fig_path: `str`, optional (default: None)
        if save_fig is True, specify figure path. if None, adata.uns['workdir'] will be used.
    fig_name: `str`, optional (default: 'branches.pdf')
        if save_fig is True, specify figure name.
    plotly: `bool`, optional (default: False)
        if True, plotly will be used to make interactive plots       

    Returns
    -------
    None

    """

    if(fig_path is None):
        fig_path = adata.uns['workdir']
    fig_size = mpl.rcParams['figure.figsize'] if fig_size is None else fig_size

    if(n_components==None):
        n_components = min(3,adata.obsm['X_dr'].shape[1])
    if n_components not in [2,3]:
        raise ValueError("n_components should be 2 or 3")

    assert (all(np.isin(['epg','flat_tree'],adata.uns_keys()))),'''graph is not learnt yet. 
    please first run: `st.seed_elastic_principal_graph` and `st.elastic_principal_graph` to learn graph'''
    assert (key_graph in ['epg','seed_epg','ori_epg']),"key_graph must be one of ['epg','seed_epg','ori_epg']"    
    if(fig_path is None):
        fig_path = adata.uns['workdir']
    if(key_graph=='epg'):
        epg = adata.uns['epg']
        flat_tree = adata.uns['flat_tree']
    else:
        epg = adata.uns[key_graph.split('_')[0]+'_epg']
        flat_tree = adata.uns[key_graph.split('_')[0]+'_flat_tree']
    ft_node_pos = nx.get_node_attributes(flat_tree,'pos')
    ft_node_label = nx.get_node_attributes(flat_tree,'label')
    epg_node_pos = nx.get_node_attributes(epg,'pos')
        
    if(plotly):
        if(n_components==3): 
            fig = go.Figure(data=go.Scatter3d())
            for edge_i in flat_tree.edges():
                branch_i_pos = np.array([epg_node_pos[i] for i in flat_tree.edges[edge_i]['nodes']])
                edge_i_label = flat_tree.nodes[edge_i[0]]['label'] +'_'+flat_tree.nodes[edge_i[1]]['label']
                curve_i = pd.DataFrame(branch_i_pos,columns=np.arange(branch_i_pos.shape[1]))
                fig.add_trace(go.Scatter3d(x=curve_i[comp1], 
                                           y=curve_i[comp2], 
                                           z=curve_i[comp3],
                                           mode='markers+lines',
                                           marker=dict(size=3,color='black'),
                                           line=dict(color='black', width=3),
                                           name=edge_i_label))
            if(show_text):
                fig.add_trace(go.Scatter3d(x=np.array(list(epg_node_pos.values()))[:,comp1], 
                                           y=np.array(list(epg_node_pos.values()))[:,comp2], 
                                           z=np.array(list(epg_node_pos.values()))[:,comp3],
                                           mode='text',
                                           opacity=1,
                                           text=[x for x in epg_node_pos.keys()],
                                           textposition="bottom center",))
                
            fig.update_layout(legend= {'itemsizing': 'constant'},width=500,height=500, 
                              scene = dict(xaxis_title='Dim'+str(comp1+1),yaxis_title='Dim'+str(comp2+1),zaxis_title='Dim'+str(comp3+1),
                                           aspectmode='cube')) 
        else:
            fig = go.Figure(data=go.Scatter())
            for edge_i in flat_tree.edges():
                branch_i_pos = np.array([epg_node_pos[i] for i in flat_tree.edges[edge_i]['nodes']])[:,[comp1,comp2]]
                edge_i_label = flat_tree.nodes[edge_i[0]]['label'] +'_'+flat_tree.nodes[edge_i[1]]['label']
                curve_i = pd.DataFrame(branch_i_pos,columns=np.arange(branch_i_pos.shape[1]))
                fig.add_trace(go.Scatter(x=curve_i[comp1], 
                                         y=curve_i[comp2],
                                         mode='markers+lines',
                                         marker=dict(size=mpl.rcParams['lines.markersize'],color='black'),                                         
                                         line=dict(color='black', width=3),
                                         name=edge_i_label)),
            if(show_text):
                fig.add_trace(go.Scatter(x=np.array(list(epg_node_pos.values()))[:,comp1], 
                                           y=np.array(list(epg_node_pos.values()))[:,comp2], 
                                           mode='text',
                                           opacity=1,
                                           text=[x for x in epg_node_pos.keys()],
                                           textposition="bottom center"),)
            fig.update_layout(legend= {'itemsizing': 'constant'},width=500,height=500, 
                              xaxis_title='Dim'+str(comp1+1),yaxis_title='Dim'+str(comp2+1))         
        fig.show(renderer="notebook")
            
    else:
        fig = plt.figure(figsize=(fig_size[0],fig_size[1]))
        if(n_components==3):
            ax_i = fig.add_subplot(1,1,1,projection='3d')
            for edge_i in flat_tree.edges():
                branch_i_pos = np.array([epg_node_pos[i] for i in flat_tree.edges[edge_i]['nodes']])
                edge_i_label = flat_tree.nodes[edge_i[0]]['label'] +'_'+flat_tree.nodes[edge_i[1]]['label']
                curve_i = pd.DataFrame(branch_i_pos,columns=range(branch_i_pos.shape[1]))
                ax_i.plot(curve_i[comp1],curve_i[comp2],curve_i[comp3],marker='o',c = 'black',
                          ms=0.8*mpl.rcParams['lines.markersize'])
            if(show_text):
                for node_i in epg.nodes():
                    ax_i.text(epg_node_pos[node_i][comp1],epg_node_pos[node_i][comp2],epg_node_pos[node_i][comp3],node_i,
                              color='black',fontsize=0.8*mpl.rcParams['font.size'],
                              ha='left', va='bottom',
                              zorder=10)                
            ax_i.set_xlabel('Dim'+str(comp1+1),labelpad=-5,rotation=-15)
            ax_i.set_ylabel('Dim'+str(comp2+1),labelpad=0,rotation=45)
            ax_i.set_zlabel('Dim'+str(comp3+1),labelpad=5,rotation=90)
            ax_i.locator_params(axis='x',nbins=4)
            ax_i.locator_params(axis='y',nbins=4)
            ax_i.locator_params(axis='z',nbins=4)
            ax_i.tick_params(axis="x",pad=-4)
            ax_i.tick_params(axis="y",pad=-1)
            ax_i.tick_params(axis="z",pad=3.5)
        else:  
            ax_i = fig.add_subplot(1,1,1)
            for edge_i in flat_tree.edges():
                branch_i_pos = np.array([epg_node_pos[i] for i in flat_tree.edges[edge_i]['nodes']])
                edge_i_label = flat_tree.nodes[edge_i[0]]['label'] +'_'+flat_tree.nodes[edge_i[1]]['label']
                curve_i = pd.DataFrame(branch_i_pos,columns=range(branch_i_pos.shape[1]))
                ax_i.plot(curve_i[comp1],curve_i[comp2],marker='o',c = 'black',
                          ms=0.8*mpl.rcParams['lines.markersize'])
            if(show_text):
                for node_i in epg.nodes():
                    ax_i.text(epg_node_pos[node_i][comp1],epg_node_pos[node_i][comp2],node_i,
                              color='black',fontsize=0.8*mpl.rcParams['font.size'],
                              ha='left', va='bottom')  
            ax_i.set_xlabel('Dim'+str(comp1+1),labelpad=2)
            ax_i.set_ylabel('Dim'+str(comp2+1),labelpad=-6)
            ax_i.locator_params(axis='x',nbins=5)
            ax_i.locator_params(axis='y',nbins=5)
            ax_i.tick_params(axis="x",pad=-1)
            ax_i.tick_params(axis="y",pad=-3)
        plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad)
        if(save_fig):
            plt.savefig(os.path.join(fig_path,fig_name),pad_inches=1,bbox_inches='tight')
            plt.close(fig)

def switch_to_low_dimension(adata,n_components=2):
    """Switch to low dimension space, in which the preliminary structure will be learnt.    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix
    n_components: `int`, optional (default: 2)
        Number of components used to infer the initial elastic principal graph.

    Returns
    -------
    adata_low: AnnData
        Annotated data matrix used in low dimensional space
    """
    if('X_dr' not in adata.obsm_keys()):
        raise ValueError('Please run dimension reduction first')
    if(n_components > adata.obsm['X_dr'].shape[1]):
        raise ValueError("`n_components` should be smaller than %s" %adata.obsm['X_dr'].shape[1])
    adata_low = adata.copy()
    adata_low.obsm['X_dr'] = adata.obsm['X_dr'][:,:n_components]
    adata_low.obsm['X_dr_ori'] = adata.obsm['X_dr']
    return adata_low


def infer_initial_structure(adata_low,nb_min=5):
    """Infer the initial node positions and edges. It helps infer the initial structure used in high-dimensional space
    
    Parameters
    ----------
    adata_low: AnnData
        Annotated data matrix used in low dimensional space
    nb_min: `int`, optional (default: 2)
        Minimum number of neighbour cells when mapping elastic principal graph from low-dimension to high-dimension.
        if the number of cells within one node is greater than nb_min, these cells will be used to calculate the new position of this node in high dimensional space
        if the number of cells within one node is less than or equal to nb_min, then nb_min nearest neighbor cells of this node will be used to calculate its new position in high dimensional space    

    Returns
    -------
    init_nodes_pos: `array`, shape = [n_nodes,n_dimension], optional (default: `None`)
        initial node positions
    init_edges: `array`, shape = [n_edges,2], optional (default: `None`)
        initial edges, all the initial nodes should be included in the tree structure
    """    
    n_components = adata_low.obsm['X_dr'].shape[1]
    epg_low = adata_low.uns['epg']
    kdtree=cKDTree(adata_low.obsm['X_dr'])
    dict_nodes_pos = dict()
    nx.set_node_attributes(epg_low,values=False,name='inferred_by_knn')
    for node_i in epg_low.nodes():
        ids = np.where(adata_low.obs['node'] == node_i)[0]
        if(ids.shape[0]<=nb_min):
            ids = kdtree.query(nx.get_node_attributes(epg_low,'pos')[node_i],k=nb_min)[1]
            print('Node '+ str(node_i) +' is calculated using ' + str(nb_min) + 'nearest neighbor cells')
            epg_low.nodes[node_i]['inferred_by_knn'] = True
        dict_nodes_pos[node_i] = np.concatenate((nx.get_node_attributes(epg_low,'pos')[node_i],
                                                 np.mean(adata_low.obsm['X_dr_ori'][ids,n_components:],axis=0)))
        init_nodes_pos = np.array(list(dict_nodes_pos.values()))
        init_edges = epg_low.edges()
    return init_nodes_pos,init_edges

def seed_elastic_principal_graph(adata,init_nodes_pos=None,init_edges=None,clustering='kmeans',damping=0.75,pref_perc=50,n_clusters=10,max_n_clusters=200,n_neighbors=50, nb_pct=None,use_vis=False):
    
    """Seeding the initial elastic principal graph.
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    init_nodes_pos: `array`, shape = [n_nodes,n_dimension], optional (default: `None`)
        initial node positions
    init_edges: `array`, shape = [n_edges,2], optional (default: `None`)
        initial edges, all the initial nodes should be included in the tree structure
    clustering: `str`, optional (default: 'kmeans')
        Choose from {{'ap','kmeans','sc'}}
        clustering method used to infer the initial nodes.
        'ap' affinity propagation
        'kmeans' K-Means clustering
        'sc' spectral clustering
    damping: `float`, optional (default: 0.75)
        Damping factor (between 0.5 and 1) for affinity propagation.
    pref_perc: `int`, optional (default: 50)
        Preference percentile (between 0 and 100). The percentile of the input similarities for affinity propagation.
    n_clusters: `int`, optional (default: 10)
        Number of clusters (only valid once 'clustering' is specificed as 'sc' or 'kmeans').
    max_n_clusters: `int`, optional (default: 200)
        The allowed maximum number of clusters for 'ap'.
    n_neighbors: `int`, optional (default: 50)
        The number of neighbor cells used for spectral clustering.
    nb_pct: `float`, optional (default: None)
        The percentage of neighbor cells (when sepcified, it will overwrite n_neighbors).
    use_vis: `bool`, optional (default: False)
        If True, the manifold learnt from st.plot_visualization_2D() will replace the manifold learnt from st.dimension_reduction().
        The principal graph will be learnt in the new manifold.

    Returns
    -------
    updates `adata` with the following fields.

    adata.obs: `pandas.core.frame.DataFrame` (`adata.obs`)
        Update adata.obs with adding the columns of current root_node_pseudotime and removing the previous ones.  
    clustering: `pandas.core.series.Series` (`adata.obs['clustering']`,dtype `str`)
        Array of dim (number of samples) that stores the clustering labels ('0', '1', …) for each cell.      
    epg : `networkx.classes.graph.Graph` (`adata.uns['epg']`)
        Elastic principal graph structure. It contains node attributes ('pos')
    flat_tree : `networkx.classes.graph.Graph` (`adata.uns['flat_tree']`)
        An abstract of elastic principle graph structure by only keeping leaf nodes and branching nodes. 
        It contains node attribtutes ('pos','label') and edge attributes ('nodes','id','len','color').
    seed_epg : `networkx.classes.graph.Graph` (`adata.uns['epg']`)
        Store seeded elastic principal graph structure
    seed_flat_tree : `networkx.classes.graph.Graph` (`adata.uns['flat_tree']`)
        Store seeded flat_tree

    Notes
    -------
    The default procedure is fast and good enough when seeding structure in low-dimensional space.

    when seeding structure in high-dimensional space, it's strongly recommended that using 'infer_initial_structure' to get the initial node positions and edges

    """

    print('Seeding initial elastic principal graph...')
    if(use_vis):
        if('X_vis' in adata.obsm_keys()):
            adata.obsm['X_dr'] = adata.obsm['X_vis']
        else:
            raise ValueError("Please run `st.plot_visualization_2D()` first")
    input_data = adata.obsm['X_dr']
    if(nb_pct!=None):
        n_neighbors = int(np.around(input_data.shape[0]*nb_pct))    
    if(init_nodes_pos is None):
        print('Clustering...')
        if(clustering=='ap'):
            print('Affinity propagation ...')
            ap = AffinityPropagation(damping=damping,random_state=42,preference=np.percentile(-euclidean_distances(input_data,squared=True),pref_perc)).fit(input_data)
            # ap = AffinityPropagation(damping=damping).fit(input_data)
            if(ap.cluster_centers_.shape[0]>max_n_clusters):
                print('The number of clusters is ' + str(ap.cluster_centers_.shape[0]))
                print('Too many clusters are generated, please lower pref_perc or increase damping and retry it')
                return
            cluster_labels = ap.labels_
            init_nodes_pos = ap.cluster_centers_
            epg_nodes_pos = init_nodes_pos  
        elif(clustering=='sc'):
            print('Spectral clustering ...')
            sc = SpectralClustering(n_clusters=n_clusters,affinity='nearest_neighbors',n_neighbors=n_neighbors,
                                    eigen_solver='arpack',random_state=42).fit(input_data)
            cluster_labels = sc.labels_ 
            init_nodes_pos = np.empty((0,input_data.shape[1])) #cluster centers
            for x in np.unique(cluster_labels):
                id_cells = np.array(range(input_data.shape[0]))[cluster_labels==x]
                init_nodes_pos = np.vstack((init_nodes_pos,np.median(input_data[id_cells,:],axis=0)))
            epg_nodes_pos = init_nodes_pos
        elif(clustering=='kmeans'):
            print('K-Means clustering ...')
            kmeans = KMeans(n_clusters=n_clusters,init='k-means++',random_state=42).fit(input_data)
            cluster_labels = kmeans.labels_
            init_nodes_pos = kmeans.cluster_centers_
            epg_nodes_pos = init_nodes_pos     
        else:
            print("'"+clustering+"'" + ' is not supported')
        adata.obs[clustering] = ['cluster '+str(x) for x in cluster_labels]
    else:
        epg_nodes_pos = init_nodes_pos
        print('Setting initial nodes...')
    print('The number of initial nodes is ' + str(epg_nodes_pos.shape[0]))

    if(init_edges is None):
        #Minimum Spanning Tree
        print('Calculatng minimum spanning tree...')
        D=pairwise_distances(epg_nodes_pos)
        G=nx.from_numpy_matrix(D)
        mst=nx.minimum_spanning_tree(G)
        epg_edges = np.array(mst.edges())
    else:
        print('Setting initial edges...')
        epg_edges = init_edges


    #store graph information and update adata
    epg=nx.Graph()
    epg.add_nodes_from(range(epg_nodes_pos.shape[0]))
    epg.add_edges_from(epg_edges)
    dict_nodes_pos = {i:x for i,x in enumerate(epg_nodes_pos)}
    nx.set_node_attributes(epg,values=dict_nodes_pos,name='pos')
    dict_branches = extract_branches(epg)
    flat_tree = construct_flat_tree(dict_branches)
    nx.set_node_attributes(flat_tree,values={x:dict_nodes_pos[x] for x in flat_tree.nodes()},name='pos')
    adata.uns['epg'] = deepcopy(epg)
    adata.uns['flat_tree'] = deepcopy(flat_tree)
    adata.uns['seed_epg'] = deepcopy(epg)
    adata.uns['seed_flat_tree'] = deepcopy(flat_tree)  
    project_cells_to_epg(adata)
    calculate_pseudotime(adata)
    if('params' not in adata.uns_keys()):
        adata.uns['params'] = dict()
    adata.uns['params']['epg'] = {'use_vis':use_vis}
    print('Number of initial branches: ' + str(len(dict_branches))) 


def elastic_principal_graph(adata,epg_n_nodes = 50,incr_n_nodes=30,epg_lambda=0.02,epg_mu=0.1,epg_trimmingradius='Inf',
                            epg_finalenergy = 'Penalized',epg_alpha=0.02,epg_beta=0.0,epg_n_processes=1,use_vis=False,
                            save_fig=False,fig_name='ElPiGraph_analysis.pdf',fig_path=None,fig_size=(8,8), **kwargs):
    """Elastic principal graph learning.
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    epg_n_nodes: `int`, optional (default: 50)
        Number of nodes for elastic principal graph.
    incr_n_nodes: `int`, optional (default: 30)
        Incremental number of nodes for elastic principal graph when epg_n_nodes is not big enough.
    epg_lambda: `float`, optional (default: 0.02)
        lambda parameter used to compute the elastic energy.
    epg_mu: `float`, optional (default: 0.1)
        mu parameter used to compute the elastic energy.
    epg_trimmingradius: `float`, optional (default: 'Inf')  
        maximal distance from a node to the points it controls in the embedding.
    epg_finalenergy: `str`, optional (default: 'Penalized')
        indicate the final elastic energy associated with the configuration.
    epg_alpha: `float`, optional (default: 0.02)
        alpha parameter of the penalized elastic energy.
    epg_beta: `float`, optional (default: 0.0)
        beta parameter of the penalized elastic energy.
    epg_n_processes: `int`, optional (default: 1)
        the number of processes to use.
    use_vis: `bool`, optional (default: False)
        Only valid when st.seed_elastic_principal_graph() is not performed.
        If True, the manifold learnt from st.plot_visualization_2D() will replace the manifold learnt from st.dimension_reduction().
        The principal graph will be learnt in the new manifold.
    save_fig: `bool`, optional (default: False)
        if True,save the figure.
    fig_size: `tuple`, optional (default: (8,8))
        figure size.
    fig_path: `str`, optional (default: None)
        if None, adata.uns['workdir'] will be used.
    fig_name: `str`, optional (default: 'ElPigraph_analysis.pdf')
        if save_fig is True, specify figure name.
    **kwargs: additional arguments to `ElPiGraph.computeElasticPrincipalTree`
   
    Returns
    -------
    updates `adata` with the following fields.
    
    adata.obs: `pandas.core.frame.DataFrame` (`adata.obs`)
        Update adata.obs with adding the columns of current root_node_pseudotime and removing the previous ones.
    epg : `networkx.classes.graph.Graph` (`adata.uns['epg']`)
        Elastic principal graph structure. It contains node attributes ('pos')
    ori_epg : `networkx.classes.graph.Graph` (`adata.uns['ori_epg']`)
        Store original elastic principal graph structure
    epg_obj : `rpy2.rinterface.ListSexpVector` (`adata.uns['epg_obj']`)
        R object of elastic principal graph learning.
    ori_epg_obj : `rpy2.rinterface.ListSexpVector` (`adata.uns['ori_epg_obj']`)
        Store original R object of elastic principal graph learning.
    flat_tree : `networkx.classes.graph.Graph` (`adata.uns['flat_tree']`)
        An abstract of elastic principle graph structure by only keeping leaf nodes and branching nodes. 
        It contains node attribtutes ('pos','label') and edge attributes ('nodes','id','len','color').
    ori_flat_tree : `networkx.classes.graph.Graph` (`adata.uns['flat_tree']`)
        Store original flat_tree
    """
    def markCardioZeroNode(epg):
        indList = [n for n in epg.nodes]
        posList = [epg.nodes[n]['pos'] for n in epg.nodes]
        posSer = pd.Series(posList, index=indList)
        disSer = posSer.apply(np.linalg.norm)
        markDict = {disSer.idxmin():{'cardioZero':''}}
        nx.set_node_attributes(epg, markDict)
        return epg
    if(fig_path is None):
        fig_path = adata.uns['workdir']

    if('seed_epg' in adata.uns_keys()):
        epg = adata.uns['seed_epg']
        dict_nodes_pos = nx.get_node_attributes(epg,'pos')
        init_nodes_pos = np.array(list(dict_nodes_pos.values()))
        init_edges = np.array(list(epg.edges())) 
        R_init_edges = init_edges + 1
        if((init_nodes_pos.shape[0]+incr_n_nodes)>=epg_n_nodes):
            print('epg_n_nodes is too small. It is corrected to the initial number of nodes plus incr_n_nodes')
            epg_n_nodes = init_nodes_pos.shape[0]+incr_n_nodes
    else:
        print('No initial structure is seeded')
        if(use_vis):
            if('X_vis' in adata.obsm_keys()):
                adata.obsm['X_dr'] = adata.obsm['X_vis']
            else:
                raise ValueError("Please run `st.plot_visualization_2D()` first")
        if('params' not in adata.uns_keys()):
            adata.uns['params'] = dict()
        adata.uns['params']['epg'] = {'use_vis':use_vis}
        init_nodes_pos = robjects.NULL
        R_init_edges = robjects.NULL

    input_data = adata.obsm['X_dr']
        
    ElPiGraph = importr('ElPiGraph.R')
    pandas2ri.activate()
    # print('Learning elastic principal graph...')
    R.pdf(os.path.join(fig_path,fig_name))
    epg_obj = ElPiGraph.computeElasticPrincipalTree(X=input_data,
                                                    NumNodes = epg_n_nodes, 
                                                    InitNodePositions = init_nodes_pos,
                                                    InitEdges=R_init_edges,
                                                    Lambda=epg_lambda, Mu=epg_mu,
                                                    TrimmingRadius= epg_trimmingradius,
                                                    FinalEnergy = epg_finalenergy,
                                                    alpha = epg_alpha,
                                                    beta = epg_beta,                                                    
                                                    Do_PCA=False,CenterData=True,
                                                    n_cores = epg_n_processes,
                                                    nReps=1,
                                                    ProbPoint=1.0,
                                                    MaxNumberOfIterations=5,
                                                    **kwargs)

    R('dev.off()')

    epg_nodes_pos = np.array(epg_obj[0].rx2('NodePositions'))
    epg_edges = np.array((epg_obj[0].rx2('Edges')).rx2('Edges'),dtype=int)-1

    #store graph information and update adata
    epg=nx.Graph()
    epg.add_nodes_from(range(epg_nodes_pos.shape[0]))
    epg.add_edges_from(epg_edges)
    dict_nodes_pos = {i:x for i,x in enumerate(epg_nodes_pos)}
    nx.set_node_attributes(epg,values=dict_nodes_pos,name='pos')
    isTwoBranches = pfTreeIsLine(epg)
    if isTwoBranches:
        epg = markCardioZeroNode(epg)
    dict_branches = extract_branches(epg)
    flat_tree = construct_flat_tree(dict_branches)
    nx.set_node_attributes(flat_tree,values={x:dict_nodes_pos[x] for x in flat_tree.nodes()},name='pos')
    adata.uns['epg'] = deepcopy(epg)
    adata.uns['ori_epg'] = deepcopy(epg)
    adata.uns['epg_obj'] = deepcopy(epg_obj)    
    adata.uns['ori_epg_obj'] = deepcopy(epg_obj)
    adata.uns['flat_tree'] = deepcopy(flat_tree)
    project_cells_to_epg(adata)
    calculate_pseudotime(adata)
    print('Number of branches after learning elastic principal graph: ' + str(len(dict_branches)))
    if('params' not in adata.uns_keys()):
        adata.uns['params'] = dict()
    adata.uns['params']['epg'].update({'epg_alpha':epg_alpha,'epg_lambda':epg_lambda,'epg_mu':epg_mu,'epg_trimmingradius':epg_trimmingradius,'epg_beta':epg_beta})


def prune_elastic_principal_graph(adata,epg_collapse_mode = 'PointNumber',epg_collapse_par = 5,   
                                  epg_lambda=None,epg_mu=None,epg_trimmingradius=None,
                                  epg_finalenergy = 'base',epg_alpha=None,epg_beta=None,epg_n_processes=1,reset=False,**kwargs): 
    """Prune the learnt elastic principal graph by filtering out 'trivial' branches.
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    epg_collapse_mode: `str`, optional (default: 'PointNumber')
        The mode used to prune the graph.
        Choose from {{'PointNumber','PointNumber_Extrema','PointNumber_Leaves','EdgesNumber','EdgesLength'}}
        'PointNumber': branches with less than epg_collapse_par points (points projected on the extreme points are not considered) are removed
        'PointNumber_Extrema', branches with less than epg_collapse_par (points projected on the extreme points are not considered) are removed
        'PointNumber_Leaves', branches with less than epg_collapse_par points (points projected on non-leaf extreme points are not considered) are removed
        'EdgesNumber', branches with less than epg_collapse_par edges are removed
        'EdgesLength', branches shorter than epg_collapse_par are removed        
    epg_collapse_par: `float`, optional (default: 5)
        The paramter used to control different modes.
    epg_lambda: `float`, optional (default: None)
        lambda parameter used to compute the elastic energy. By default using the same `epg_lambda` from `elastic_principal_graph()`
    epg_mu: `float`, optional (default: None)
        mu parameter used to compute the elastic energy. By default using the same `epg_mu` from `elastic_principal_graph()`
    epg_trimmingradius: `float`, optional (default: None)  
        maximal distance from a node to the points it controls in the embedding. By default using the same `epg_trimmingradius` from `elastic_principal_graph()`
    epg_finalenergy: `str`, optional (default: 'Penalized')
        indicate the final elastic energy associated with the configuration.
    epg_alpha: `float`, optional (default: None)
        alpha parameter of the penalized elastic energy. By default using the same `epg_alpha` from `elastic_principal_graph()`
    epg_beta: `float`, optional (default: None)
        beta parameter of the penalized elastic energy. By default using the same `epg_beta` from `elastic_principal_graph()`
    epg_n_processes: `int`, optional (default: 1)
        The number of processes to use.
    reset: `bool`, optional (default: False)
        If true, reset the current elastic principal graph to the initial elastic principal graph (i.e. the graph obtained from running 'elastic_principal_graph')
    **kwargs: additional arguments to `ElPiGraph.CollapseBrances`
   
    Returns
    -------
    updates `adata` with the following fields.

    adata.obs: `pandas.core.frame.DataFrame` (`adata.obs`)
        Update adata.obs with adding the columns of current root_node_pseudotime and removing the previous ones.    
    epg : `networkx.classes.graph.Graph` (`adata.uns['epg']`)
        Elastic principal graph structure. It contains node attributes ('pos')
    epg_obj : `rpy2.rinterface.ListSexpVector` (`adata.uns['epg_obj']`)
        R object of elastic principal graph learning.
    flat_tree : `networkx.classes.graph.Graph` (`adata.uns['flat_tree']`)
        An abstract of elastic principle graph structure by only keeping leaf nodes and branching nodes. 
        It contains node attribtutes ('pos','label') and edge attributes ('nodes','id','len','color').
    """

    print('Collasping small branches ...')
    if(epg_alpha is None):
        epg_alpha = adata.uns['params']['epg']['epg_alpha']
    if(epg_beta is None):
        epg_beta = adata.uns['params']['epg']['epg_beta']
    if(epg_lambda is None):
        epg_lambda = adata.uns['params']['epg']['epg_lambda']
    if(epg_mu is None):
        epg_mu = adata.uns['params']['epg']['epg_mu']   
    if(epg_trimmingradius is None):
        epg_trimmingradius = adata.uns['params']['epg']['epg_trimmingradius']

    ElPiGraph = importr('ElPiGraph.R')
    pandas2ri.activate()
    if(reset):
        epg_obj = adata.uns['ori_epg_obj']
        epg = adata.uns['ori_epg']
    else:
        epg_obj = adata.uns['epg_obj']
        epg = adata.uns['epg']
    if(len(extract_branches(epg))<3):
        print("No branching points are detected. This step is skipped")
        return

    input_data = adata.obsm['X_dr']
    epg_obj_collapse = ElPiGraph.CollapseBrances(X = input_data, TargetPG = epg_obj[0], Mode = epg_collapse_mode, ControlPar = epg_collapse_par, **kwargs)

    init_nodes_pos = np.array(epg_obj_collapse.rx2('Nodes'))
    init_edges = np.array(epg_obj_collapse.rx2('Edges')) - 1    
    epg_n_nodes = init_nodes_pos.shape[0]
    epg_obj = ElPiGraph.computeElasticPrincipalTree(X=input_data,
                                                    NumNodes = epg_n_nodes, 
                                                    InitNodePositions = init_nodes_pos,
                                                    InitEdges=init_edges + 1,
                                                    Lambda=epg_lambda, Mu=epg_mu,
                                                    TrimmingRadius= epg_trimmingradius,
                                                    FinalEnergy = epg_finalenergy,
                                                    alpha = epg_alpha,
                                                    beta = epg_beta,                                                    
                                                    Do_PCA=False,CenterData=False,
                                                    drawAccuracyComplexity = False, drawEnergy = False,drawPCAView = False,
                                                    n_cores = epg_n_processes,
                                                    nReps=1,
                                                    ProbPoint=1.0)
    
    epg_nodes_pos = np.array(epg_obj[0].rx2('NodePositions'))
    epg_edges = np.array((epg_obj[0].rx2('Edges')).rx2('Edges'),dtype=int)-1    
        
    #store graph information and update adata
    epg=nx.Graph()
    epg.add_nodes_from(range(epg_nodes_pos.shape[0]))
    epg.add_edges_from(epg_edges)
    dict_nodes_pos = {i:x for i,x in enumerate(epg_nodes_pos)}
    nx.set_node_attributes(epg,values=dict_nodes_pos,name='pos')
    dict_branches = extract_branches(epg)
    flat_tree = construct_flat_tree(dict_branches)
    nx.set_node_attributes(flat_tree,values={x:dict_nodes_pos[x] for x in flat_tree.nodes()},name='pos')
    adata.uns['epg'] = deepcopy(epg)
    adata.uns['epg_obj'] = deepcopy(epg_obj)
    adata.uns['flat_tree'] = deepcopy(flat_tree)
    project_cells_to_epg(adata)
    calculate_pseudotime(adata)
    print('Number of branches after pruning ElPiGraph: ' + str(len(dict_branches)))


def optimize_branching(adata,incr_n_nodes=30,epg_maxsteps=50,mode=2,                                  
                       epg_lambda=None,epg_mu=None,epg_trimmingradius=None,
                       epg_finalenergy = 'base',epg_alpha=None,epg_beta=None,epg_n_processes=1,reset=False,**kwargs):
    """Optimize branching node by expanding the nodes around a branching point.
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    incr_n_nodes: `int`, optional (default: 30)
        Incremental number of nodes for elastic principal graph.       
    epg_maxsteps: `float`, optional (default: 50)
        The maximum number of iteration steps .
    mode: `int`, optional (default: 2)
        The energy computation mode.
    epg_lambda: `float`, optional (default: None)
        lambda parameter used to compute the elastic energy. By default using the same `epg_lambda` from `elastic_principal_graph()`
    epg_mu: `float`, optional (default: None)
        mu parameter used to compute the elastic energy. By default using the same `epg_mu` from `elastic_principal_graph()`
    epg_trimmingradius: `float`, optional (default: None)  
        maximal distance from a node to the points it controls in the embedding. By default using the same `epg_trimmingradius` from `elastic_principal_graph()`
    epg_finalenergy: `str`, optional (default: 'Penalized')
        indicate the final elastic energy associated with the configuration.
    epg_alpha: `float`, optional (default: 0.02)
        alpha parameter of the penalized elastic energy. By default using the same `epg_alpha` from `elastic_principal_graph()`
    epg_beta: `float`, optional (default: 0.0)
        beta parameter of the penalized elastic energy. By default using the same `epg_beta` from `elastic_principal_graph()`
    epg_n_processes: `int`, optional (default: 1)
        The number of processes to use.
    reset: `bool`, optional (default: False)
        If true, reset the current elastic principal graph to the initial elastic principal graph (i.e. the graph obtained from running 'elastic_principal_graph')
    **kwargs: additional arguments to `ElPiGraph.CollapseBrances`
   
    Returns
    -------
    updates `adata` with the following fields.

    adata.obs: `pandas.core.frame.DataFrame` (`adata.obs`)
        Update adata.obs with adding the columns of current root_node_pseudotime and removing the previous ones.        
    epg : `networkx.classes.graph.Graph` (`adata.uns['epg']`)
        Elastic principal graph structure. It contains node attributes ('pos')
    epg_obj : `rpy2.rinterface.ListSexpVector` (`adata.uns['epg_obj']`)
        R object of elastic principal graph learning.
    flat_tree : `networkx.classes.graph.Graph` (`adata.uns['flat_tree']`)
        An abstract of elastic principle graph structure by only keeping leaf nodes and branching nodes. 
        It contains node attribtutes ('pos','label') and edge attributes ('nodes','id','len','color').
    """

    print('Optimizing branching...')
    if(epg_alpha is None):
        epg_alpha = adata.uns['params']['epg']['epg_alpha']
    if(epg_beta is None):
        epg_beta = adata.uns['params']['epg']['epg_beta']
    if(epg_lambda is None):
        epg_lambda = adata.uns['params']['epg']['epg_lambda']
    if(epg_mu is None):
        epg_mu = adata.uns['params']['epg']['epg_mu']   
    if(epg_trimmingradius is None):
        epg_trimmingradius = adata.uns['params']['epg']['epg_trimmingradius']

    ElPiGraph = importr('ElPiGraph.R')
    pandas2ri.activate()
    if(reset):
        epg_obj = adata.uns['ori_epg_obj']
        epg = adata.uns['ori_epg']
    else:
        epg_obj = adata.uns['epg_obj']
        epg = adata.uns['epg']
    if(len(extract_branches(epg))<3):
        print("No branching points are detected. This step is skipped")
        return

    input_data = adata.obsm['X_dr']

    dict_nodes_pos = nx.get_node_attributes(epg,'pos')
    init_nodes_pos = np.array(list(dict_nodes_pos.values()))
    init_edges = np.array(list(epg.edges()))  
    epg_n_nodes = init_nodes_pos.shape[0] + incr_n_nodes
    
    epg_obj = ElPiGraph.fineTuneBR(X=input_data,
                                    MaxSteps = epg_maxsteps,
                                    Mode = 2,
                                    NumNodes = epg_n_nodes, 
                                    InitNodePositions = init_nodes_pos,
                                    InitEdges=init_edges + 1,
                                    Lambda=epg_lambda, Mu=epg_mu,
                                    TrimmingRadius= epg_trimmingradius,
                                    FinalEnergy = epg_finalenergy,
                                    alpha = epg_alpha,
                                    beta = epg_beta,                                                    
                                    Do_PCA=False,CenterData=False,
                                    drawAccuracyComplexity = False, drawEnergy = False,drawPCAView = False,
                                    n_cores = epg_n_processes,
                                    nReps=1,
                                    ProbPoint=1.0,
                                    **kwargs)
    
    epg_nodes_pos = np.array(epg_obj[0].rx2('NodePositions'))
    epg_edges = np.array((epg_obj[0].rx2('Edges')).rx2('Edges'),dtype=int)-1    
        
    #store graph information and update adata
    epg=nx.Graph()
    epg.add_nodes_from(range(epg_nodes_pos.shape[0]))
    epg.add_edges_from(epg_edges)
    dict_nodes_pos = {i:x for i,x in enumerate(epg_nodes_pos)}
    nx.set_node_attributes(epg,values=dict_nodes_pos,name='pos')
    dict_branches = extract_branches(epg)
    flat_tree = construct_flat_tree(dict_branches)
    nx.set_node_attributes(flat_tree,values={x:dict_nodes_pos[x] for x in flat_tree.nodes()},name='pos')
    adata.uns['epg'] = deepcopy(epg)
    adata.uns['epg_obj'] = deepcopy(epg_obj)
    adata.uns['flat_tree'] = deepcopy(flat_tree)
    project_cells_to_epg(adata)
    calculate_pseudotime(adata)
    print('Number of branches after optimizing branching: ' + str(len(dict_branches)))   


def shift_branching(adata,epg_shift_mode = 'NodeDensity',epg_shift_radius = 0.05,epg_shift_max=5,                             
                   epg_lambda=None,epg_mu=None,epg_trimmingradius=None,
                   epg_finalenergy = 'base',epg_alpha=None,epg_beta=None,epg_n_processes=1,reset=False,**kwargs):
    """Move branching node to the area with higher density.
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    epg_shift_mode: `str`, optional (default: 'NodeDensity')
        The mode used to shift the branching nodes.
        Choose from {{'NodePoints','NodeDensity'}}
    epg_shift_radius: `float`, optional (default: 0.05)
        The radius used when computing point density if epg_shift_mode = 'NodeDensity'.
    epg_shift_max: `float`, optional (default: 5)
        The maxium distance (defined as the number of edges) to consider when exploring the neighborhood of branching point
    epg_lambda: `float`, optional (default: None)
        lambda parameter used to compute the elastic energy. By default using the same `epg_lambda` from `elastic_principal_graph()`
    epg_mu: `float`, optional (default: None)
        mu parameter used to compute the elastic energy. By default using the same `epg_mu` from `elastic_principal_graph()`
    epg_trimmingradius: `float`, optional (default: None)  
        maximal distance from a node to the points it controls in the embedding. By default using the same `epg_trimmingradius` from `elastic_principal_graph()`
    epg_finalenergy: `str`, optional (default: 'Penalized')
        indicate the final elastic energy associated with the configuration.
    epg_alpha: `float`, optional (default: None)
        alpha parameter of the penalized elastic energy. By default using the same `epg_alpha` from `elastic_principal_graph()`
    epg_beta: `float`, optional (default: None)
        beta parameter of the penalized elastic energy. By default using the same `epg_beta` from `elastic_principal_graph()`
    epg_n_processes: `int`, optional (default: 1)
        The number of processes to use.
    reset: `bool`, optional (default: False)
        If true, reset the current elastic principal graph to the initial elastic principal graph (i.e. the graph obtained from running 'elastic_principal_graph')
    **kwargs: additional arguments to `ElPiGraph.CollapseBrances`
   
    Returns
    -------
    updates `adata` with the following fields.

    adata.obs: `pandas.core.frame.DataFrame` (`adata.obs`)
        Update adata.obs with adding the columns of current root_node_pseudotime and removing the previous ones.        
    epg : `networkx.classes.graph.Graph` (`adata.uns['epg']`)
        Elastic principal graph structure. It contains node attributes ('pos')
    epg_obj : `rpy2.rinterface.ListSexpVector` (`adata.uns['epg_obj']`)
        R object of elastic principal graph learning.
    flat_tree : `networkx.classes.graph.Graph` (`adata.uns['flat_tree']`)
        An abstract of elastic principle graph structure by only keeping leaf nodes and branching nodes. 
        It contains node attribtutes ('pos','label') and edge attributes ('nodes','id','len','color').
    """

    print('Shifting branching point to denser area ...')
    if(epg_alpha is None):
        epg_alpha = adata.uns['params']['epg']['epg_alpha']
    if(epg_beta is None):
        epg_beta = adata.uns['params']['epg']['epg_beta']
    if(epg_lambda is None):
        epg_lambda = adata.uns['params']['epg']['epg_lambda']
    if(epg_mu is None):
        epg_mu = adata.uns['params']['epg']['epg_mu']   
    if(epg_trimmingradius is None):
        epg_trimmingradius = adata.uns['params']['epg']['epg_trimmingradius']

    ElPiGraph = importr('ElPiGraph.R')
    pandas2ri.activate()
    if(reset):
        epg_obj = adata.uns['ori_epg_obj']
        epg = adata.uns['ori_epg']
    else:
        epg_obj = adata.uns['epg_obj']
        epg = adata.uns['epg']
    if(len(extract_branches(epg))<3):
        print("No branching points are detected. This step is skipped")
        return

    input_data = adata.obsm['X_dr']
    epg_obj_shift = ElPiGraph.ShiftBranching(X = input_data, 
                                           TargetPG = epg_obj[0], 
                                           TrimmingRadius = epg_trimmingradius,                       
                                           SelectionMode = epg_shift_mode, 
                                           DensityRadius = epg_shift_radius,
                                           MaxShift = epg_shift_max,
                                           **kwargs)

    init_nodes_pos = np.array(epg_obj_shift.rx2('NodePositions'))
    init_edges = np.array(epg_obj_shift.rx2('Edges')) - 1  
    epg_n_nodes = init_nodes_pos.shape[0]
 
    epg_obj = ElPiGraph.computeElasticPrincipalTree(X=input_data,
                                                    NumNodes = epg_n_nodes, 
                                                    InitNodePositions = init_nodes_pos,
                                                    InitEdges=init_edges + 1,
                                                    Lambda=epg_lambda, Mu=epg_mu,
                                                    TrimmingRadius= epg_trimmingradius,
                                                    FinalEnergy = epg_finalenergy,
                                                    alpha = epg_alpha,
                                                    beta = epg_beta,                                                    
                                                    Do_PCA=False,CenterData=False,
                                                    drawAccuracyComplexity = False, drawEnergy = False,drawPCAView = False,
                                                    n_cores = epg_n_processes,
                                                    nReps=1,
                                                    ProbPoint=1.0)
    
    epg_nodes_pos = np.array(epg_obj[0].rx2('NodePositions'))
    epg_edges = np.array((epg_obj[0].rx2('Edges')).rx2('Edges'),dtype=int)-1   
                    
    #store graph information and update adata
    epg=nx.Graph()
    epg.add_nodes_from(range(epg_nodes_pos.shape[0]))
    epg.add_edges_from(epg_edges)
    dict_nodes_pos = {i:x for i,x in enumerate(epg_nodes_pos)}
    nx.set_node_attributes(epg,values=dict_nodes_pos,name='pos')
    dict_branches = extract_branches(epg)
    flat_tree = construct_flat_tree(dict_branches)
    nx.set_node_attributes(flat_tree,values={x:dict_nodes_pos[x] for x in flat_tree.nodes()},name='pos')
    adata.uns['epg'] = deepcopy(epg)
    adata.uns['epg_obj'] = deepcopy(epg_obj)
    adata.uns['flat_tree'] = deepcopy(flat_tree)
    project_cells_to_epg(adata)
    calculate_pseudotime(adata)
    print('Number of branches after shifting branching: ' + str(len(dict_branches)))


def extend_elastic_principal_graph(adata,epg_ext_mode = 'QuantDists',epg_ext_par = 0.5,epg_trimmingradius=None,reset=False,**kwargs):
    """Extend the leaves of elastic principal graph with additional nodes.
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    epg_ext_mode: `str`, optional (default: 'QuantDists')
        The mode used to extend the leaves.
        Choose from {{'QuantDists','QuantCentroid','WeigthedCentroid'}}
        'QuantCentroid':for each leaf node, the extreme points are ordered by their distance from the node and the centroid of the points further than epg_ext_par is returned.
        'WeigthedCentroid':for each leaf node, a weight is computed for each points by raising the distance to the epg_ext_par power. Larger epg_ext_par results in a bigger influence of points further than the node
        'QuantDists':for each leaf node, the extreme points are ordered by their distance from the node and the 100*epg_ext_par th percentile of the points farther than epg_ext_par is returned
    epg_ext_par: `float`, optional (default: 0.5)
        The paramter used to control different modes.
    epg_trimmingradius: `float`, optional (default: None)  
        maximal distance from a node to the points it controls in the embedding. By default using the same `epg_trimmingradius` from `elastic_principal_graph()`
    reset: `bool`, optional (default: False)
        If true, reset the current elastic principal graph to the initial elastic principal graph (i.e. the graph obtained from running 'elastic_principal_graph')
    **kwargs: additional arguments to `ElPiGraph.CollapseBrances`
   
    Returns
    -------
    updates `adata` with the following fields.

    adata.obs: `pandas.core.frame.DataFrame` (`adata.obs`)
        Update adata.obs with adding the columns of current root_node_pseudotime and removing the previous ones.        
    epg : `networkx.classes.graph.Graph` (`adata.uns['epg']`)
        Elastic principal graph structure. It contains node attributes ('pos')
    epg_obj : `rpy2.rinterface.ListSexpVector` (`adata.uns['epg_obj']`)
        R object of elastic principal graph learning.
    flat_tree : `networkx.classes.graph.Graph` (`adata.uns['flat_tree']`)
        An abstract of elastic principle graph structure by only keeping leaf nodes and branching nodes. 
        It contains node attribtutes ('pos','label') and edge attributes ('nodes','id','len','color').
    """

    print('Extending leaves with additional nodes ...') 
    if(epg_trimmingradius is None):
        epg_trimmingradius = adata.uns['params']['epg']['epg_trimmingradius']

    ElPiGraph = importr('ElPiGraph.R')
    pandas2ri.activate()
    if(reset):
        epg_obj = adata.uns['ori_epg_obj']
        epg = adata.uns['ori_epg']
    else:
        epg_obj = adata.uns['epg_obj']
        epg = adata.uns['epg']

    input_data = adata.obsm['X_dr']
    epg_obj_extend = ElPiGraph.ExtendLeaves(X = input_data, 
                                          TargetPG = epg_obj[0],
                                          TrimmingRadius = epg_trimmingradius,
                                          Mode = epg_ext_mode, 
                                          ControlPar = epg_ext_par,
                                          PlotSelected = False,
                                          **kwargs)
    epg_nodes_pos = np.array(epg_obj_extend.rx2('NodePositions'))
    epg_edges = np.array((epg_obj_extend.rx2('Edges')).rx2('Edges'),dtype=int)-1   
        
    #store graph information and update adata
    epg=nx.Graph()
    epg.add_nodes_from(range(epg_nodes_pos.shape[0]))
    epg.add_edges_from(epg_edges)
    dict_nodes_pos = {i:x for i,x in enumerate(epg_nodes_pos)}
    nx.set_node_attributes(epg,values=dict_nodes_pos,name='pos')
    dict_branches = extract_branches(epg)
    flat_tree = construct_flat_tree(dict_branches)
    nx.set_node_attributes(flat_tree,values={x:dict_nodes_pos[x] for x in flat_tree.nodes()},name='pos')
    adata.uns['epg'] = deepcopy(epg)
    adata.uns['flat_tree'] = deepcopy(flat_tree)
#     adata.uns['epg_obj'] = deepcopy(epg_obj_extend)
    project_cells_to_epg(adata)
    calculate_pseudotime(adata)
    print('Number of branches after extending leaves: ' + str(len(dict_branches)))    


def plot_flat_tree(adata,color=None,dist_scale=1,
                   fig_size=None,fig_ncol=3,fig_legend_ncol=1,fig_legend_order = None,
                   vmin=None,vmax=None,alpha=0.8,
                   pad=1.08,w_pad=None,h_pad=None,
                   show_text=False,show_graph=False,
                   save_fig=False,fig_path=None,fig_name='flat_tree.pdf',
                   plotly=False):  
    """Plot flat tree based on a modified version of the force-directed layout Fruchterman-Reingold algorithm.
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    color: `list` optional (default: None)
        Column names of observations (adata.obs.columns) or variable names(adata.var_names). A list of names to be plotted. 
    dist_scale: `float`,optional (default: 1)
        Scaling factor to scale the distance from cells to tree branches 
        (by default, it keeps the same distance as in original manifold)
    fig_size: `tuple`, optional (default: None)
        figure size.
    fig_legend_order: `dict`,optional (default: None)
        Specified order for the appearance of the annotation keys.Only valid for ategorical variable  
        e.g. fig_legend_order = {'ann1':['a','b','c'],'ann2':['aa','bb','cc']}
    fig_legend_ncol: `int`, optional (default: 1)
        The number of columns that the legend has.
    vmin,vmax: `float`, optional (default: None)
        The min and max values are used to normalize continuous values. If None, the respective min and max of continuous values is used.
    alpha: `float`, optional (default: 0.8)
        0.0 transparent through 1.0 opaque
    pad: `float`, optional (default: 1.08)
        Padding between the figure edge and the edges of subplots, as a fraction of the font size.
    h_pad, w_pad: `float`, optional (default: None)
        Padding (height/width) between edges of adjacent subplots, as a fraction of the font size. Defaults to pad.
    show_text: `bool`, optional (default: False)
        If True, node state label will be shown
    show_graph: `bool`, optional (default: False)
        If True, the learnt principal graph will be shown
    save_fig: `bool`, optional (default: False)
        if True,save the figure.
    fig_path: `str`, optional (default: None)
        if save_fig is True, specify figure path. if None, adata.uns['workdir'] will be used.
    fig_name: `str`, optional (default: 'flat_tree.pdf')
        if save_fig is True, specify figure name.
    plotly: `bool`, optional (default: False)
        if True, plotly will be used to make interactive plots.

    Returns
    -------
    updates `adata` with the following fields.
    X_spring: `numpy.ndarray` (`adata.obsm['X_spring']`)
        Store #observations × 2 coordinates of cells in flat tree.
    """

    if(fig_path is None):
        fig_path = adata.uns['workdir']
    fig_size = mpl.rcParams['figure.figsize'] if fig_size is None else fig_size
       
    if(color is None):
        color = ['label']
    ###remove duplicate keys
    color = list(dict.fromkeys(color))     

    dict_ann = dict()
    for ann in color:
        if(ann in adata.obs.columns):
            dict_ann[ann] = adata.obs[ann]
        elif(ann in adata.var_names):
            dict_ann[ann] = adata.obs_vector(ann)
        else:
            raise ValueError("could not find '%s' in `adata.obs.columns` and `adata.var_names`"  % (ann))

    ## add the positions of flat tree's nodes
    add_flat_tree_node_pos(adata)
    flat_tree = adata.uns['flat_tree']
    ## add the positions of cells on flat tre
    add_flat_tree_cell_pos(adata,dist_scale)
    
    ft_node_pos = nx.get_node_attributes(flat_tree,'pos_spring')
    ft_node_label = nx.get_node_attributes(flat_tree,'label')
    
    df_plot = pd.DataFrame(index=adata.obs.index,data = adata.obsm['X_spring'],columns=['FlatTree'+str(x+1) for x in range(adata.obsm['X_spring'].shape[1])])
    for ann in color:
        df_plot[ann] = dict_ann[ann]
    df_plot_shuf = df_plot.sample(frac=1,random_state=100)
    
    legend_order = {ann:np.unique(df_plot_shuf[ann]) for ann in color if is_string_dtype(df_plot_shuf[ann])}
    if(fig_legend_order is not None):
        if(not isinstance(fig_legend_order, dict)):
            raise TypeError("`fig_legend_order` must be a dictionary")
        for ann in fig_legend_order.keys():
            if(ann in legend_order.keys()):
                legend_order[ann] = fig_legend_order[ann]
            else:
                print("'%s' is ignored for ordering legend labels due to incorrect name or data type" % ann)
        
    if(plotly):
        for ann in color:
            fig = px.scatter(df_plot_shuf, x='FlatTree1', y='FlatTree2',color=ann,
                                 opacity=alpha,width=500,height=500,
                                 color_continuous_scale=px.colors.sequential.Viridis,
                                 color_discrete_map=adata.uns[ann+'_color'] if ann+'_color' in adata.uns_keys() else {})
            if(show_graph):
                for edge_i in flat_tree.edges():
                    branch_i_pos = np.array([ft_node_pos[i] for i in edge_i])
                    branch_i = pd.DataFrame(branch_i_pos,columns=range(branch_i_pos.shape[1]))
                    edge_i_label = flat_tree.nodes[edge_i[0]]['label'] +'_'+flat_tree.nodes[edge_i[1]]['label']
                    fig.add_trace(go.Scatter(x=branch_i[0], 
                                               y=branch_i[1],
                                               mode='lines',
                                               line=dict(color='black', width=3),
                                               name=edge_i_label))
            if(show_text):
                fig.add_trace(go.Scatter(x=np.array(list(ft_node_pos.values()))[:,0], 
                                           y=np.array(list(ft_node_pos.values()))[:,1], 
                                           mode='markers+text',
#                                            mode='text',
                                           opacity=1,
                                           marker=dict(size=1.5*mpl.rcParams['lines.markersize'],color='#767070'),
                                           text=[ft_node_label[x] for x in ft_node_pos.keys()],
                                           textposition="bottom center",
                                           name='states'),)

        fig.update_layout(legend= {'itemsizing': 'constant'},width=500,height=500) 
        fig.show(renderer="notebook")
    else:
        if(len(color)<fig_ncol):
            fig_ncol=len(color)
        fig_nrow = int(np.ceil(len(color)/fig_ncol))
        fig = plt.figure(figsize=(fig_size[0]*fig_ncol*1.05,fig_size[1]*fig_nrow))
        for i,ann in enumerate(color):
            ax_i = fig.add_subplot(fig_nrow,fig_ncol,i+1)
            if(is_string_dtype(df_plot[ann])):
                sc_i=sns.scatterplot(ax=ax_i,
                                    x='FlatTree1', y='FlatTree2',
                                    hue=ann,hue_order = legend_order[ann],
                                    data=df_plot_shuf,
                                    alpha=alpha,linewidth=0,
                                    palette= adata.uns[ann+'_color'] \
                                            if (ann+'_color' in adata.uns_keys()) and (set(adata.uns[ann+'_color'].keys()) >= set(np.unique(df_plot_shuf[ann]))) \
                                            else None)
                legend_handles, legend_labels = ax_i.get_legend_handles_labels()
                ax_i.legend(handles=legend_handles, labels=legend_labels,
                            bbox_to_anchor=(1, 0.5), loc='center left', ncol=fig_legend_ncol,
                            frameon=False,
                            )
                if(ann+'_color' not in adata.uns_keys()):
                    colors_sns = sc_i.get_children()[0].get_facecolors()
                    colors_sns_scaled = (255*colors_sns).astype(int)
                    adata.uns[ann+'_color'] = {df_plot_shuf[ann][i]:'#%02x%02x%02x' % (colors_sns_scaled[i][0], colors_sns_scaled[i][1], colors_sns_scaled[i][2])
                                               for i in np.unique(df_plot_shuf[ann],return_index=True)[1]}
                ### remove legend title
                # ax_i.get_legend().texts[0].set_text("")
            else:
                vmin_i = df_plot[ann].min() if vmin is None else vmin
                vmax_i = df_plot[ann].max() if vmax is None else vmax
                sc_i = ax_i.scatter(df_plot_shuf['FlatTree1'], df_plot_shuf['FlatTree2'],
                                    c=df_plot_shuf[ann],vmin=vmin_i,vmax=vmax_i,alpha=alpha)
                cbar = plt.colorbar(sc_i,ax=ax_i, pad=0.01, fraction=0.05, aspect=40)
                cbar.solids.set_edgecolor("face")
                cbar.ax.locator_params(nbins=5)
            if(show_graph):
                for edge_i in flat_tree.edges():
                    branch_i_pos = np.array([ft_node_pos[i] for i in edge_i])
                    branch_i = pd.DataFrame(branch_i_pos,columns=range(branch_i_pos.shape[1]))
                    ax_i.plot(branch_i[0],branch_i[1],c = 'black',alpha=0.8)
            if(show_text):
                for node_i in flat_tree.nodes():
                    ax_i.text(ft_node_pos[node_i][0],ft_node_pos[node_i][1],ft_node_label[node_i],
                              color='black',fontsize=0.9*mpl.rcParams['font.size'],
                               ha='left', va='bottom')  
            ax_i.set_xlabel("FlatTree1",labelpad=2)
            ax_i.set_ylabel("FlatTree2",labelpad=0)
            ax_i.locator_params(axis='x',nbins=5)
            ax_i.locator_params(axis='y',nbins=5)
            ax_i.tick_params(axis="x",pad=-1)
            ax_i.tick_params(axis="y",pad=-3)
            ax_i.set_title(ann)
            plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad)
        if(save_fig):
            plt.savefig(os.path.join(fig_path,fig_name),pad_inches=1,bbox_inches='tight')
            plt.close(fig)            

def plot_visualization_2D(adata,method='umap',n_neighbors=50, nb_pct=None,perplexity=30.0,color=None,use_precomputed=True,
                          fig_size=None,fig_ncol=3,fig_legend_ncol=1,fig_legend_order = None,
                          vmin=None,vmax=None,alpha=0.8,
                          pad=1.08,w_pad=None,h_pad=None,
                          save_fig=False,fig_path=None,fig_name='visualization_2D.pdf',
                          plotly=False):

    """ Visualize the results in 2D plane
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix. 
    method: `str`, optional (default: 'umap')
        Choose from {{'umap','tsne'}}
        Method used for visualization.
        'umap': Uniform Manifold Approximation and Projection      
        'tsne': t-Distributed Stochastic Neighbor Embedding
    n_neighbors: `int`, optional (default: 50)
        The number of neighbor cells (only valid when 'umap' is specified).
    nb_pct: `float`, optional (default: None)
        The percentage of neighbor cells (when sepcified, it will overwrite n_neighbors).
    perplexity: `float`, optional (default: 30.0)
        The perplexity used for tSNE. 
    color: `list` optional (default: None)
        Column names of observations (adata.obs.columns) or variable names(adata.var_names). A list of names to be plotted.    
    use_precomputed: `bool`, optional (default: True)
        If True, the visualization coordinates from previous computing will be used
    fig_size: `tuple`, optional (default: None)
        figure size.
    fig_legend_order: `dict`,optional (default: None)
        Specified order for the appearance of the annotation keys.Only valid for ategorical variable  
        e.g. fig_legend_order = {'ann1':['a','b','c'],'ann2':['aa','bb','cc']}
    fig_legend_ncol: `int`, optional (default: 1)
        The number of columns that the legend has.
    vmin,vmax: `float`, optional (default: None)
        The min and max values are used to normalize continuous values. If None, the respective min and max of continuous values is used.
    alpha: `float`, optional (default: 0.8)
        0.0 transparent through 1.0 opaque
    pad: `float`, optional (default: 1.08)
        Padding between the figure edge and the edges of subplots, as a fraction of the font size.
    h_pad, w_pad: `float`, optional (default: None)
        Padding (height/width) between edges of adjacent subplots, as a fraction of the font size. Defaults to pad.
    save_fig: `bool`, optional (default: False)
        if True,save the figure.
    fig_path: `str`, optional (default: None)
        if save_fig is True, specify figure path. if None, adata.uns['workdir'] will be used.
    fig_name: `str`, optional (default: 'visualization_2D.pdf')
        if save_fig is True, specify figure name.
    plotly: `bool`, optional (default: False)
        if True, plotly will be used to make interactive plots

    Returns
    -------
    updates `adata` with the following fields. (Depending on `method`)
    vis_trans_umap : `umap.UMAP` (`adata.uns['vis_trans_umap']`)
        Store umap object
    vis_trans_tsne : `sklearn.manifold._t_sne.TSNE` (`adata.uns['vis_trans_tsne']`)
        Store tsne object
    X_vis_umap: `numpy.ndarray` (`adata.obsm['X_vis_umap']`)
        Store #observations × 2 umap data matrix. 
    X_vis_tsne: `numpy.ndarray` (`adata.obsm['X_vis_tsne']`)
        Store #observations × 2 tsne data matrix. 
    X_vis : `numpy.ndarray` (`adata.obsm['X_vis']`)
        A #observations × 2 data matrix after visualization.   
    """    
    if(fig_path is None):
        fig_path = adata.uns['workdir']
    fig_size = mpl.rcParams['figure.figsize'] if fig_size is None else fig_size

    if(method not in ['umap','tsne']):
        raise ValueError("unrecognized method '%s'" % method)
    if(color is None):
        color = ['label']
    ###remove duplicate keys
    color = list(dict.fromkeys(color))

    dict_ann = dict()
    for ann in color:
        if(ann in adata.obs.columns):
            dict_ann[ann] = adata.obs[ann]
        elif(ann in adata.var_names):
            dict_ann[ann] = adata.obs_vector(ann)
        else:
            raise ValueError("could not find '%s' in `adata.obs.columns` and `adata.var_names`"  % (ann))
    input_data = adata.obsm['X_dr']
    if(nb_pct!=None):
        n_neighbors = int(np.around(input_data.shape[0]*nb_pct)) 
    if(method == 'umap'):       
        if(use_precomputed and ('X_vis_umap' in adata.obsm_keys())):
            print('Importing precomputed umap visualization ...')
            embedding = adata.obsm['X_vis_umap']
        else:
            reducer = umap.UMAP(n_neighbors=n_neighbors,n_components=2,random_state=42)
            trans = reducer.fit(input_data)
            embedding = trans.embedding_
            adata.uns['vis_trans_umap'] = trans
            adata.obsm['X_vis_umap'] = embedding
    if(method == 'tsne'):
        if(use_precomputed and ('X_vis_tsne' in adata.obsm_keys())):
            print('Importing precomputed tsne visualization ...')
            embedding = adata.obsm['X_vis_tsne']
        else:
            reducer = TSNE(n_components=2, init='pca',perplexity=perplexity, random_state=42)
            trans = reducer.fit(input_data)
            embedding = trans.embedding_
            adata.uns['vis_trans_tsne'] = trans
            adata.obsm['X_vis_tsne'] = embedding
    adata.obsm['X_vis'] = embedding
    if('params' not in adata.uns_keys()):
        adata.uns['params'] = dict()
    adata.uns['params']['plot_visualization_2D'] = {'method':method}
    
    df_plot = pd.DataFrame(index=adata.obs.index,data = embedding,columns=[method.upper()+str(x) for x in [1,2]])
    for ann in color:
        df_plot[ann] = dict_ann[ann]
    df_plot_shuf = df_plot.sample(frac=1,random_state=100)
    
    legend_order = {ann:np.unique(df_plot_shuf[ann]) for ann in color if is_string_dtype(df_plot_shuf[ann])}
    if(fig_legend_order is not None):
        if(not isinstance(fig_legend_order, dict)):
            raise TypeError("`fig_legend_order` must be a dictionary")
        for ann in fig_legend_order.keys():
            if(ann in legend_order.keys()):
                legend_order[ann] = fig_legend_order[ann]
            else:
                print("'%s' is ignored for ordering legend labels due to incorrect name or data type" % ann)

    if(plotly):
        for ann in color:
            fig = px.scatter(df_plot_shuf, x='Dim'+str(comp1+1), y='Dim'+str(comp2+1),color=ann,
                                opacity=alpha,width=500,height=500,
                                color_continuous_scale=px.colors.sequential.Viridis,
                                color_discrete_map=adata.uns[ann+'_color'] if ann+'_color' in adata.uns_keys() else {})
            fig.update_layout(legend= {'itemsizing': 'constant'}) 
            fig.show(renderer="notebook")
    else:
        if(len(color)<fig_ncol):
            fig_ncol=len(color)
        fig_nrow = int(np.ceil(len(color)/fig_ncol))
        fig = plt.figure(figsize=(fig_size[0]*fig_ncol*1.05,fig_size[1]*fig_nrow))
        for i,ann in enumerate(color):
            ax_i = fig.add_subplot(fig_nrow,fig_ncol,i+1)
            if(is_string_dtype(df_plot[ann])):
                sc_i=sns.scatterplot(ax=ax_i,
                                    x=method.upper()+'1', y=method.upper()+'2', 
                                    hue=ann,hue_order = legend_order[ann],
                                    data=df_plot_shuf,
                                    alpha=alpha,linewidth=0,
                                    palette= adata.uns[ann+'_color'] \
                                            if (ann+'_color' in adata.uns_keys()) and (set(adata.uns[ann+'_color'].keys()) >= set(np.unique(df_plot_shuf[ann]))) \
                                            else None
                                    )
                legend_handles, legend_labels = ax_i.get_legend_handles_labels()
                ax_i.legend(handles=legend_handles, labels=legend_labels,
                            bbox_to_anchor=(1, 0.5), loc='center left', ncol=fig_legend_ncol,
                            frameon=False,
                            )
                if(ann+'_color' not in adata.uns_keys()):
                    colors_sns = sc_i.get_children()[0].get_facecolors()
                    colors_sns_scaled = (255*colors_sns).astype(int)
                    adata.uns[ann+'_color'] = {df_plot_shuf[ann][i]:'#%02x%02x%02x' % (colors_sns_scaled[i][0], colors_sns_scaled[i][1], colors_sns_scaled[i][2])
                                                for i in np.unique(df_plot_shuf[ann],return_index=True)[1]}
                ### remove legend title
                # ax_i.get_legend().texts[0].set_text("")
            else:
                vmin_i = df_plot[ann].min() if vmin is None else vmin
                vmax_i = df_plot[ann].max() if vmax is None else vmax
                sc_i = ax_i.scatter(df_plot_shuf[method.upper()+'1'], df_plot_shuf[method.upper()+'2'],
                                    c=df_plot_shuf[ann],vmin=vmin_i,vmax=vmax_i,alpha=alpha)
                cbar = plt.colorbar(sc_i,ax=ax_i, pad=0.01, fraction=0.05, aspect=40)
                cbar.solids.set_edgecolor("face")
                cbar.ax.locator_params(nbins=5)                    
            ax_i.set_xlabel(method.upper()+'1')
            ax_i.set_ylabel(method.upper()+'2',labelpad=2)
            ax_i.get_xaxis().set_ticks([])
            ax_i.get_yaxis().set_ticks([])
            ax_i.set_title(ann)
#             plt.subplots_adjust(hspace=hspace,wspace=wspace)
        plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad)
        if(save_fig):
            plt.savefig(os.path.join(fig_path,fig_name),pad_inches=1,bbox_inches='tight')
            plt.close(fig)

def plot_stream_sc(adata,root='S0',color=None,dist_scale=1,dist_pctl=95,preference=None,
                   fig_size=(7,4.5),fig_legend_ncol=1,fig_legend_order = None,
                   vmin=None,vmax=None,alpha=0.8,
                   pad=1.08,w_pad=None,h_pad=None,
                   show_text=True,show_graph=True,
                   save_fig=False,fig_path=None,fig_format='pdf',
                   plotly=False): 
    """Generate stream plot at single cell level (aka, subway map plots)
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    root: `str`, optional (default: 'S0'): 
        The starting node
    color: `list` optional (default: None)
        Column names of observations (adata.obs.columns) or variable names(adata.var_names). A list of names to be plotted. 
    dist_scale: `float`,optional (default: 1)
        Scaling factor to scale the distance from cells to tree branches 
        (by default, it keeps the same distance as in original manifold)
    dist_pctl: `int`, optional (default: 95)
        Percentile of cells' distances from branches (between 0 and 100) used for calculating the distances between branches.
    preference: `list`, optional (default: None): 
        The preference of nodes. The branch with speficied nodes are preferred and put on the top part of stream plot. 
        The higher ranks the node have, the closer to the top the branch with that node is.
    fig_size: `tuple`, optional (default: (7,4.5))
        figure size.
    fig_legend_order: `dict`,optional (default: None)
        Specified order for the appearance of the annotation keys.Only valid for ategorical variable  
        e.g. fig_legend_order = {'ann1':['a','b','c'],'ann2':['aa','bb','cc']}
    fig_legend_ncol: `int`, optional (default: 1)
        The number of columns that the legend has.
    vmin,vmax: `float`, optional (default: None)
        The min and max values are used to normalize continuous values. If None, the respective min and max of continuous values is used.
    alpha: `float`, optional (default: 0.8)
        0.0 transparent through 1.0 opaque
    pad: `float`, optional (default: 1.08)
        Padding between the figure edge and the edges of subplots, as a fraction of the font size.
    h_pad, w_pad: `float`, optional (default: None)
        Padding (height/width) between edges of adjacent subplots, as a fraction of the font size. Defaults to pad.
    show_text: `bool`, optional (default: False)
        If True, node state label will be shown
    show_graph: `bool`, optional (default: False)
        If True, the learnt principal graph will be shown
    save_fig: `bool`, optional (default: False)
        if True,save the figure.
    fig_path: `str`, optional (default: None)
        if save_fig is True, specify figure path. if None, adata.uns['workdir'] will be used.
    fig_format: `str`, optional (default: 'pdf')
        if save_fig is True, specify figure format.
    plotly: `bool`, optional (default: False)
        if True, plotly will be used to make interactive plots

    Returns
    -------
    updates `adata` with the following fields.
    X_stream_root: `numpy.ndarray` (`adata.obsm['X_stream_root']`)
        Store #observations × 2 coordinates of cells in subwaymap plot.
    stream_root: `dict` (`adata.uns['stream_root']`)
        Store the coordinates of nodes ('nodes') and edges ('edges') in subwaymap plot.
    """

    if(fig_path is None):
        fig_path = adata.uns['workdir']
    fig_size = mpl.rcParams['figure.figsize'] if fig_size is None else fig_size

    if(color is None):
        color = ['label']
    ###remove duplicate keys
    color = list(dict.fromkeys(color))     

    dict_ann = dict()
    for ann in color:
        if(ann in adata.obs.columns):
            dict_ann[ann] = adata.obs[ann]
        elif(ann in adata.var_names):
            dict_ann[ann] = adata.obs_vector(ann)
        else:
            raise ValueError("could not find '%s' in `adata.obs.columns` and `adata.var_names`"  % (ann))
    
    flat_tree = adata.uns['flat_tree']
    ft_node_label = nx.get_node_attributes(flat_tree,'label')
    label_to_node = {value: key for key,value in nx.get_node_attributes(flat_tree,'label').items()}    
    if(root not in label_to_node.keys()):
        raise ValueError("There is no root '%s'" % root)  
            
    add_stream_sc_pos(adata,root=root,dist_scale=dist_scale,dist_pctl=dist_pctl,preference=preference)
    stream_nodes = adata.uns['stream_'+root]['nodes']
    stream_edges = adata.uns['stream_'+root]['edges']

    df_plot = pd.DataFrame(index=adata.obs.index,data = adata.obsm['X_stream_'+root],
                           columns=['pseudotime','dist'])
    for ann in color:
        df_plot[ann] = dict_ann[ann]
    df_plot_shuf = df_plot.sample(frac=1,random_state=100)

    legend_order = {ann:np.unique(df_plot_shuf[ann]) for ann in color if is_string_dtype(df_plot_shuf[ann])}
    if(fig_legend_order is not None):
        if(not isinstance(fig_legend_order, dict)):
            raise TypeError("`fig_legend_order` must be a dictionary")
        for ann in fig_legend_order.keys():
            if(ann in legend_order.keys()):
                legend_order[ann] = fig_legend_order[ann]
            else:
                print("'%s' is ignored for ordering legend labels due to incorrect name or data type" % ann)        

    if(plotly):
        for ann in color:
            fig = px.scatter(df_plot_shuf, x='pseudotime', y='dist',color=ann,
                                 opacity=alpha,
                                 color_continuous_scale=px.colors.sequential.Viridis,
                                 color_discrete_map=adata.uns[ann+'_color'] if ann+'_color' in adata.uns_keys() else {})
            if(show_graph):
                for edge_i in stream_edges.keys():
                    branch_i_pos = stream_edges[edge_i]
                    branch_i = pd.DataFrame(branch_i_pos,columns=range(branch_i_pos.shape[1]))
                    for ii in np.arange(start=0,stop=branch_i.shape[0],step=2):
                        if(branch_i.iloc[ii,0]==branch_i.iloc[ii+1,0]):
                            fig.add_trace(go.Scatter(x=branch_i.iloc[[ii,ii+1],0], 
                                                       y=branch_i.iloc[[ii,ii+1],1],
                                                       mode='lines',
                                                       opacity=0.8,
                                                       line=dict(color='#767070', width=3),
                                                       showlegend=False))
                        else:
                            fig.add_trace(go.Scatter(x=branch_i.iloc[[ii,ii+1],0], 
                                                       y=branch_i.iloc[[ii,ii+1],1],
                                                       mode='lines',
                                                       line=dict(color='black', width=3),
                                                       showlegend=False))
            if(show_text):
                fig.add_trace(go.Scatter(x=np.array(list(stream_nodes.values()))[:,0], 
                                           y=np.array(list(stream_nodes.values()))[:,1],
                                           mode='text',
                                           opacity=1,
                                           marker=dict(size=1.5*mpl.rcParams['lines.markersize'],color='#767070'),
                                           text=[ft_node_label[x] for x in stream_nodes.keys()],
                                           textposition="bottom center",
                                           name='states',
                                           showlegend=False),)
            fig.update_layout(legend= {'itemsizing': 'constant'},
                              xaxis={'showgrid': False,'zeroline': False,},
                              yaxis={'visible':False},
                              width=800,height=500) 
            fig.show(renderer="notebook")            
    else:
        for i,ann in enumerate(color):
            fig = plt.figure(figsize=(fig_size[0],fig_size[1]))
            ax_i = fig.add_subplot(1,1,1)
            if(is_string_dtype(df_plot[ann])):
                sc_i=sns.scatterplot(ax=ax_i,
                                    x='pseudotime', y='dist',
                                    hue=ann,hue_order = legend_order[ann],
                                    data=df_plot_shuf,
                                    alpha=alpha,linewidth=0,
                                    palette= adata.uns[ann+'_color'] \
                                            if (ann+'_color' in adata.uns_keys()) and (set(adata.uns[ann+'_color'].keys()) >= set(np.unique(df_plot_shuf[ann]))) \
                                            else None                                    
                                    )
                legend_handles, legend_labels = ax_i.get_legend_handles_labels()
                ax_i.legend(handles=legend_handles, labels=legend_labels,
                            bbox_to_anchor=(1, 0.5), loc='center left', ncol=fig_legend_ncol,
                            frameon=False,
                            )
                if(ann+'_color' not in adata.uns_keys()):
                    colors_sns = sc_i.get_children()[0].get_facecolors()
                    colors_sns_scaled = (255*colors_sns).astype(int)
                    adata.uns[ann+'_color'] = {df_plot_shuf[ann][i]:'#%02x%02x%02x' % (colors_sns_scaled[i][0], colors_sns_scaled[i][1], colors_sns_scaled[i][2])
                                               for i in np.unique(df_plot_shuf[ann],return_index=True)[1]}
                ### remove legend title
                # ax_i.get_legend().texts[0].set_text("")
            else:
                vmin_i = df_plot[ann].min() if vmin is None else vmin
                vmax_i = df_plot[ann].max() if vmax is None else vmax
                sc_i = ax_i.scatter(df_plot_shuf['pseudotime'], df_plot_shuf['dist'],
                                    c=df_plot_shuf[ann],vmin=vmin_i,vmax=vmax_i,alpha=alpha)
                cbar = plt.colorbar(sc_i,ax=ax_i, pad=0.01, fraction=0.05, aspect=40)
                cbar.solids.set_edgecolor("face")
                cbar.ax.locator_params(nbins=5)
            if(show_graph):
                for edge_i in stream_edges.keys():
                    branch_i_pos = stream_edges[edge_i]
                    branch_i = pd.DataFrame(branch_i_pos,columns=range(branch_i_pos.shape[1]))
                    for ii in np.arange(start=0,stop=branch_i.shape[0],step=2):
                        if(branch_i.iloc[ii,0]==branch_i.iloc[ii+1,0]):
                            ax_i.plot(branch_i.iloc[[ii,ii+1],0],branch_i.iloc[[ii,ii+1],1],
                                      c = '#767070',alpha=0.8)
                        else:
                            ax_i.plot(branch_i.iloc[[ii,ii+1],0],branch_i.iloc[[ii,ii+1],1],
                                      c = 'black',alpha=1)
            if(show_text):
                for node_i in flat_tree.nodes():
                    ax_i.text(stream_nodes[node_i][0],stream_nodes[node_i][1],ft_node_label[node_i],
                              color='black',fontsize=0.9*mpl.rcParams['font.size'],
                               ha='left', va='bottom')  
            ax_i.set_xlabel("pseudotime",labelpad=2)
            ax_i.spines['left'].set_visible(False) 
            ax_i.spines['right'].set_visible(False)
            ax_i.spines['top'].set_visible(False) 
            ax_i.get_yaxis().set_visible(False)
            ax_i.locator_params(axis='x',nbins=8)
            ax_i.tick_params(axis="x",pad=-1)
            annots = arrowed_spines(ax_i, locations=('bottom right',),
                                    lw=ax_i.spines['bottom'].get_linewidth()*1e-5)
            ax_i.set_title(ann)
            plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad)
            if(save_fig):
                file_path_S = os.path.join(fig_path,root)
                if(not os.path.exists(file_path_S)):
                    os.makedirs(file_path_S) 
                plt.savefig(os.path.join(file_path_S,'stream_sc_' + slugify(ann) + '.' + fig_format),pad_inches=1,bbox_inches='tight')
                plt.close(fig)

def plot_stream(adata,root='S0',color = None,preference=None,dist_scale=0.9,
                factor_num_win=10,factor_min_win=2.0,factor_width=2.5,factor_nrow=200,factor_ncol=400,
                log_scale = False,factor_zoomin=100.0,
                fig_size=(7,4.5),fig_legend_order=None,fig_legend_ncol=1,
                fig_colorbar_aspect=30,
                vmin=None,vmax=None,
                pad=1.08,w_pad=None,h_pad=None,
                save_fig=False,fig_path=None,fig_format='pdf'):  
    """Generate stream plot at density level
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    root: `str`, optional (default: 'S0'): 
        The starting node
    color: `list` optional (default: None)
        Column names of observations (adata.obs.columns) or variable names(adata.var_names). A list of names to be plotted. 
    preference: `list`, optional (default: None): 
        The preference of nodes. The branch with speficied nodes are preferred and put on the top part of stream plot. 
        The higher ranks the node have, the closer to the top the branch with that node is.
    dist_scale: `float`,optional (default: 0.9)
        Scaling factor. It controls the width of STREAM plot branches. The smaller, the thinner the branch will be.
    factor_num_win: `int`, optional (default: 10)
        Number of sliding windows used for making stream plot. It controls the smoothness of STREAM plot.
    factor_min_win: `float`, optional (default: 2.0)
        The minimum number of sliding windows. It controls the resolution of STREAM plot. The window size is calculated based on shortest branch. (suggested range: 1.5~3.0)
    factor_width: `float`, optional (default: 2.5)
        The ratio between length and width of stream plot. 
    factor_nrow: `int`, optional (default: 200)
        The number of rows in the array used to plot continuous values 
    factor_ncol: `int`, optional (default: 400)
        The number of columns in the array used to plot continuous values
    log_scale: `bool`, optional (default: False)
        If True,the number of cells (the width) is logarithmized when drawing stream plot.
    factor_zoomin: `float`, optional (default: 100.0)
        If log_scale is True, the factor used to zoom in the thin branches
    fig_size: `tuple`, optional (default: (7,4.5))
        figure size.
    fig_legend_order: `dict`,optional (default: None)
        Specified order for the appearance of the annotation keys.Only valid for ategorical variable  
        e.g. fig_legend_order = {'ann1':['a','b','c'],'ann2':['aa','bb','cc']}
    fig_legend_ncol: `int`, optional (default: 1)
        The number of columns that the legend has.
    vmin,vmax: `float`, optional (default: None)
        The min and max values are used to normalize continuous values. If None, the respective min and max of continuous values is used.
    pad: `float`, optional (default: 1.08)
        Padding between the figure edge and the edges of subplots, as a fraction of the font size.
    h_pad, w_pad: `float`, optional (default: None)
        Padding (height/width) between edges of adjacent subplots, as a fraction of the font size. Defaults to pad.
    save_fig: `bool`, optional (default: False)
        if True,save the figure.
    fig_path: `str`, optional (default: None)
        if save_fig is True, specify figure path. if None, adata.uns['workdir'] will be used.
    fig_format: `str`, optional (default: 'pdf')
        if save_fig is True, specify figure format.

    Returns
    -------
    None

    """

    if(fig_path is None):
        fig_path = adata.uns['workdir']
    fig_size = mpl.rcParams['figure.figsize'] if fig_size is None else fig_size

    if(color is None):
        color = ['label']
    ###remove duplicate keys
    color = list(dict.fromkeys(color))     

    dict_ann = dict()
    for ann in color:
        if(ann in adata.obs.columns):
            dict_ann[ann] = adata.obs[ann]
        elif(ann in adata.var_names):
            dict_ann[ann] = adata.obs_vector(ann)
        else:
            raise ValueError("could not find '%s' in `adata.obs.columns` and `adata.var_names`"  % (ann))

    flat_tree = adata.uns['flat_tree']
    ft_node_label = nx.get_node_attributes(flat_tree,'label')
    label_to_node = {value: key for key,value in nx.get_node_attributes(flat_tree,'label').items()}    
    if(root not in label_to_node.keys()):
        raise ValueError("There is no root '%s'" % root)  

    if(preference!=None):
        preference_nodes = [label_to_node[x] for x in preference]
    else:
        preference_nodes = None

    legend_order = {ann:np.unique(dict_ann[ann]) for ann in color if is_string_dtype(dict_ann[ann])}
    if(fig_legend_order is not None):
        if(not isinstance(fig_legend_order, dict)):
            raise TypeError("`fig_legend_order` must be a dictionary")
        for ann in fig_legend_order.keys():
            if(ann in legend_order.keys()):
                legend_order[ann] = fig_legend_order[ann]
            else:
                print("'%s' is ignored for ordering legend labels due to incorrect name or data type" % ann)

    dict_plot = dict()
    
    list_string_type = [k for k,v in dict_ann.items() if is_string_dtype(v)]
    if(len(list_string_type)>0):
        dict_verts,dict_extent = \
        cal_stream_polygon_string(adata,dict_ann,root=root,preference=preference,dist_scale=dist_scale,
                                  factor_num_win=factor_num_win,factor_min_win=factor_min_win,factor_width=factor_width,
                                  log_scale=log_scale,factor_zoomin=factor_zoomin)  
        dict_plot['string'] = [dict_verts,dict_extent]

    list_numeric_type = [k for k,v in dict_ann.items() if is_numeric_dtype(v)]
    if(len(list_numeric_type)>0):
        verts,extent,ann_order,dict_ann_df,dict_im_array = \
        cal_stream_polygon_numeric(adata,dict_ann,root=root,preference=preference,dist_scale=dist_scale,
                                   factor_num_win=factor_num_win,factor_min_win=factor_min_win,factor_width=factor_width,
                                   factor_nrow=factor_nrow,factor_ncol=factor_ncol,
                                   log_scale=log_scale,factor_zoomin=factor_zoomin)     
        dict_plot['numeric'] = [verts,extent,ann_order,dict_ann_df,dict_im_array]
        
    for ann in color:  
        if(is_string_dtype(dict_ann[ann])):
            if(not ((ann+'_color' in adata.uns_keys()) and (set(adata.uns[ann+'_color'].keys()) >= set(np.unique(dict_ann[ann]))))):
                ### a hacky way to generate colors from seaborn
                tmp = pd.DataFrame(index=adata.obs_names,
                                   data=np.random.rand(adata.shape[0], 2))
                tmp[ann] = dict_ann[ann]
                fig = plt.figure(figsize=fig_size)
                sc_i=sns.scatterplot(x=0,y=1,hue=ann,data=tmp,linewidth=0)
                colors_sns = sc_i.get_children()[0].get_facecolors()
                plt.close(fig)
                colors_sns_scaled = (255*colors_sns).astype(int)
                adata.uns[ann+'_color'] = {tmp[ann][i]:'#%02x%02x%02x' % (colors_sns_scaled[i][0], colors_sns_scaled[i][1], colors_sns_scaled[i][2])
                                           for i in np.unique(tmp[ann],return_index=True)[1]}
            dict_palette = adata.uns[ann+'_color']

            verts = dict_plot['string'][0][ann]
            extent = dict_plot['string'][1][ann]
            xmin = extent['xmin']
            xmax = extent['xmax']
            ymin = extent['ymin'] - (extent['ymax'] - extent['ymin'])*0.1
            ymax = extent['ymax'] + (extent['ymax'] - extent['ymin'])*0.1            
            
            fig = plt.figure(figsize=fig_size)
            ax = fig.add_subplot(1,1,1)
            legend_labels = []
            for ann_i in legend_order[ann]:
                legend_labels.append(ann_i)
                verts_cell = verts[ann_i]
                polygon = Polygon(verts_cell,closed=True,color=dict_palette[ann_i],alpha=0.8,lw=0)
                ax.add_patch(polygon)
            ax.legend(legend_labels,bbox_to_anchor=(1.03, 0.5), loc='center left', ncol=fig_legend_ncol,frameon=False,  
                      columnspacing=0.4,
                      borderaxespad=0.2,
                      handletextpad=0.3,)        
        else:
            verts = dict_plot['numeric'][0] 
            extent = dict_plot['numeric'][1]
            ann_order = dict_plot['numeric'][2]
            dict_ann_df = dict_plot['numeric'][3]  
            dict_im_array = dict_plot['numeric'][4]
            xmin = extent['xmin']
            xmax = extent['xmax']
            ymin = extent['ymin'] - (extent['ymax'] - extent['ymin'])*0.1
            ymax = extent['ymax'] + (extent['ymax'] - extent['ymin'])*0.1

            #clip parts according to determined polygon
            fig = plt.figure(figsize=fig_size)
            ax = fig.add_subplot(1,1,1)
            for ann_i in ann_order:
                vmin_i = dict_ann_df[ann].loc[ann_i,:].min() if vmin is None else vmin
                vmax_i = dict_ann_df[ann].loc[ann_i,:].max() if vmax is None else vmax
                im = ax.imshow(dict_im_array[ann][ann_i],interpolation='bicubic',
                               extent=[xmin,xmax,ymin,ymax],vmin=vmin_i,vmax=vmax_i,aspect='auto') 
                verts_cell = verts[ann_i]
                clip_path = Polygon(verts_cell, facecolor='none', edgecolor='none', closed=True)
                ax.add_patch(clip_path)
                im.set_clip_path(clip_path)
                cbar = plt.colorbar(im, ax=ax, pad=0.04, fraction=0.02, aspect=fig_colorbar_aspect)
                cbar.ax.locator_params(nbins=5)  
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        ax.set_xlabel("pseudotime",labelpad=2)
        ax.spines['left'].set_visible(False) 
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False) 
        ax.get_yaxis().set_visible(False)
        ax.locator_params(axis='x',nbins=8)
        ax.tick_params(axis="x",pad=-1)
        annots = arrowed_spines(ax, locations=('bottom right',),
                                lw=ax.spines['bottom'].get_linewidth()*1e-5)
        ax.set_title(ann)
        plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad)             
        if(save_fig):
            file_path_S = os.path.join(fig_path,root)
            if(not os.path.exists(file_path_S)):
                os.makedirs(file_path_S) 
            plt.savefig(os.path.join(file_path_S,'stream_' + slugify(ann) + '.' + fig_format),pad_inches=1,bbox_inches='tight')
            plt.close(fig)

def detect_transition_markers(adata,marker_list=None,cutoff_spearman=0.4, cutoff_logfc = 0.25, percentile_expr=95, n_jobs = 1,min_num_cells=5,
                               use_precomputed=True, root='S0',preference=None):

    """Detect transition markers along one branch.
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    marker_list: `list`, optional (default: None): 
        A list of candidate markers to be scanned. Instead of scanning all available genes/peaks/kmers/motifs, this will limit the scanning to a specific list of genes/peaks/kmers/motifs
        If none, all available features (genes/peaks/kmers/motifs) will be scanned.
    cutoff_spearman: `float`, optional (default: 0.4)
        Between 0 and 1. The cutoff used for Spearman's rank correlation coefficient.
    cutoff_logfc: `float`, optional (default: 0.25)
        The log2-transformed fold change cutoff between cells around start and end node.
    percentile_expr: `int`, optional (default: 95)
        Between 0 and 100. Between 0 and 100. Specify the percentile of marker expression greater than 0 to filter out some extreme marker expressions. 
    min_num_cells: `int`, optional (default: 5)
        The minimum number of cells in which markers are expressed.
    n_jobs: `int`, optional (default: 1)
        The number of parallel jobs to run when scaling the marker expressions .
    use_precomputed: `bool`, optional (default: True)
        If True, the previously computed scaled marker expression will be used
    root: `str`, optional (default: 'S0'): 
        The starting node
    preference: `list`, optional (default: None): 
        The preference of nodes. The branch with speficied nodes are preferred and will be dealt with first. The higher ranks the node have, The earlier the branch with that node will be analyzed. 
        This will help generate the consistent results as shown in subway map and stream plot.

    Returns
    -------
    updates `adata` with the following fields.
    scaled_marker_expr: `list` (`adata.uns['scaled_marker_expr']`)
        Scaled marker expression for marker marker detection.    
    transition_markers: `dict` (`adata.uns['transition_markers']`)
        Transition markers for each branch deteced by STREAM.
    """

    file_path = os.path.join(adata.uns['workdir'],'transition_markers')
    if(not os.path.exists(file_path)):
        os.makedirs(file_path)    

    if(marker_list is None):
        print('Scanning all features ...')
        marker_list = adata.var_names.tolist()
    else:
        print('Scanning the specified marker list ...')
        ###remove duplicate keys
        marker_list = list(dict.fromkeys(marker_list)) 
        for marker in marker_list:
            if(marker not in adata.var_names):
                raise ValueError("could not find '%s' in `adata.var_names`"  % (marker))
    
    flat_tree = adata.uns['flat_tree']
    dict_node_state = nx.get_node_attributes(flat_tree,'label')
    df_marker_detection = adata.obs.copy()
    df_marker_detection.rename(columns={"label": "CELL_LABEL", "branch_lam": "lam"},inplace = True)
    if(use_precomputed and ('scaled_marker_expr' in adata.uns_keys())):
        print('Importing precomputed scaled marker expression matrix ...')
        results = adata.uns['scaled_marker_expr']    
        df_results = pd.DataFrame(results).T 
        if(all(np.isin(marker_list,df_results.columns.tolist()))):
            input_markers_expressed = marker_list
            df_scaled_marker_expr = df_results[input_markers_expressed]
        else:
            raise ValueError("Not all markers in `marker_list` can be found in precomputed scaled marker expression matrix. Please set `use_precomputed=False`")
        
    else:
        input_markers = marker_list
        df_sc = pd.DataFrame(index= adata.obs_names.tolist(),
                             data = adata[:,input_markers].X,
                             columns=input_markers)
        #exclude markers that are expressed in fewer than min_num_cells cells
        #min_num_cells = max(5,int(round(df_marker_detection.shape[0]*0.001)))
        # print('Minimum number of cells expressing markers: '+ str(min_num_cells))

        print("Filtering out markers that are expressed in less than " + str(min_num_cells) + " cells ...")
        input_markers_expressed = np.array(input_markers)[np.where((df_sc[input_markers]>0).sum(axis=0)>min_num_cells)[0]].tolist()
        df_marker_detection[input_markers_expressed] = df_sc[input_markers_expressed].copy()

        print(str(n_jobs)+' cpus are being used ...')
        params = [(df_marker_detection,x,percentile_expr) for x in input_markers_expressed]
        pool = multiprocessing.Pool(processes=n_jobs)
        results = pool.map(scale_marker_expr,params)
        pool.close()
        adata.uns['scaled_marker_expr'] = results
        df_scaled_marker_expr = pd.DataFrame(results).T

    print(str(len(input_markers_expressed)) + ' markers are being scanned ...')
    df_marker_detection[input_markers_expressed] = df_scaled_marker_expr
    #### TG (Transition markers) along each branch
    dict_tg_edges = dict()
    dict_label_node = {value: key for key,value in nx.get_node_attributes(flat_tree,'label').items()}
    if(preference!=None):
        preference_nodes = [dict_label_node[x] for x in preference]
    else:
        preference_nodes = None
    root_node = dict_label_node[root]
    bfs_edges = bfs_edges_modified(flat_tree,root_node,preference=preference_nodes)
#     all_branches = np.unique(df_marker_detection['branch_id']).tolist()
    for edge_i in bfs_edges:
        edge_i_alias = (dict_node_state[edge_i[0]],dict_node_state[edge_i[1]])
        if edge_i in nx.get_edge_attributes(flat_tree,'id').values():
            df_cells_edge_i = deepcopy(df_marker_detection[df_marker_detection.branch_id==edge_i])
            df_cells_edge_i['lam_ordered'] = df_cells_edge_i['lam']
        else:
            df_cells_edge_i = deepcopy(df_marker_detection[df_marker_detection.branch_id==(edge_i[1],edge_i[0])])
            df_cells_edge_i['lam_ordered'] = flat_tree.edges[edge_i]['len'] - df_cells_edge_i['lam']
        df_cells_edge_i_sort = df_cells_edge_i.sort_values(['lam_ordered'])
        df_stat_pval_qval = pd.DataFrame(columns = ['stat','logfc','pval','qval'],dtype=float)
        for markername in input_markers_expressed:
            id_initial = range(0,int(df_cells_edge_i_sort.shape[0]*0.2))
            id_final = range(int(df_cells_edge_i_sort.shape[0]*0.8),int(df_cells_edge_i_sort.shape[0]*1))
            values_initial = df_cells_edge_i_sort.iloc[id_initial,:][markername]
            values_final = df_cells_edge_i_sort.iloc[id_final,:][markername]
            diff_initial_final = abs(values_final.mean() - values_initial.mean())
            if(diff_initial_final>0):
                logfc = np.log2(max(values_final.mean(),values_initial.mean())/(min(values_final.mean(),values_initial.mean())+diff_initial_final/1000.0))
            else:
                logfc = 0
            if(logfc>cutoff_logfc):
                df_stat_pval_qval.loc[markername] = np.nan
                df_stat_pval_qval.loc[markername,['stat','pval']] = spearmanr(df_cells_edge_i_sort.loc[:,markername],\
                                                                            df_cells_edge_i_sort.loc[:,'lam_ordered'])
                df_stat_pval_qval.loc[markername,'logfc'] = logfc
        if(df_stat_pval_qval.shape[0]==0):
            print('No Transition markers are detected in branch ' + edge_i_alias[0]+'_'+edge_i_alias[1])
        else:
            p_values = df_stat_pval_qval['pval']
            q_values = multipletests(p_values, method='fdr_bh')[1]
            df_stat_pval_qval['qval'] = q_values
            dict_tg_edges[edge_i_alias] = df_stat_pval_qval[(abs(df_stat_pval_qval.stat)>=cutoff_spearman)].sort_values(['qval'])
            dict_tg_edges[edge_i_alias].to_csv(os.path.join(file_path,'transition_markers_'+ edge_i_alias[0]+'_'+edge_i_alias[1] + '.tsv'),sep = '\t',index = True)
    adata.uns['transition_markers'] = dict_tg_edges   


def plot_transition_markers(adata,num_markers = 15,
                            save_fig=False,fig_path=None,fig_size=None):
    if(fig_path is None):
        fig_path = os.path.join(adata.uns['workdir'],'transition_markers')
    if(not os.path.exists(fig_path)):
        os.makedirs(fig_path)

    dict_tg_edges = adata.uns['transition_markers']
    flat_tree = adata.uns['flat_tree']
    colors = sns.color_palette("Set1", n_colors=8, desat=0.8)
    for edge_i in dict_tg_edges.keys():

        df_tg_edge_i = deepcopy(dict_tg_edges[edge_i])
        df_tg_edge_i = df_tg_edge_i.iloc[:num_markers,:]

        stat = df_tg_edge_i.stat[::-1]
        qvals = df_tg_edge_i.qval[::-1]

        pos = np.arange(df_tg_edge_i.shape[0])-1
        bar_colors = np.tile(colors[4],(len(stat),1))
        id_neg = np.arange(len(stat))[np.array(stat<0)]
        bar_colors[id_neg]=colors[2]

        if(fig_size is None):
            fig_size = (12,np.ceil(0.4*len(stat)))
        fig = plt.figure(figsize=fig_size)
        ax = fig.add_subplot(1,1,1, adjustable='box')
        ax.barh(pos,stat,align='center',height=0.8,tick_label=[''],color = bar_colors)
        ax.set_xlabel('Spearman correlation coefficient')
        ax.set_title("branch " + edge_i[0]+'_'+edge_i[1])
        ax.spines['left'].set_visible(False) 
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False) 
        ax.set_ylim((min(pos)-1,max(pos)+1))

        rects = ax.patches
        for i,rect in enumerate(rects):
            if(stat[i]>0):
                alignment = {'horizontalalignment': 'left', 'verticalalignment': 'center'}
                ax.text(rect.get_x()+rect.get_width()+0.01, rect.get_y() + rect.get_height()/2.0, \
                        qvals.index[i],**alignment)
                ax.text(rect.get_x()+0.01, rect.get_y()+rect.get_height()/2.0,
                        "{:.2E}".format(Decimal(str(qvals[i]))),size=0.8*mpl.rcParams['font.size'],
                        color='black',**alignment)
            else:
                alignment = {'horizontalalignment': 'right', 'verticalalignment': 'center'}
                ax.text(rect.get_x()+rect.get_width()-0.01, rect.get_y()+rect.get_height()/2.0, \
                        qvals.index[i],**alignment)
                ax.text(rect.get_x()-0.01, rect.get_y()+rect.get_height()/2.0,
                        "{:.2E}".format(Decimal(str(qvals[i]))),size=0.8*mpl.rcParams['font.size'],
                        color='w',**alignment)
        plt.tight_layout()
        if(save_fig):        
            plt.savefig(os.path.join(fig_path,'transition_markers_'+ edge_i[0]+'_'+edge_i[1]+'.pdf'),\
                        pad_inches=1,bbox_inches='tight')
            plt.close(fig)


def detect_de_markers(adata,marker_list=None,cutoff_zscore=1,cutoff_logfc = 0.25,percentile_expr=95,n_jobs = 1,min_num_cells=5,
                      use_precomputed=True, root='S0',preference=None):

    """Detect differentially expressed markers between different sub-branches.
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    marker_list: `list`, optional (default: None): 
        A list of candidate markers to be scanned. Instead of scanning all available genes/peaks/kmers/motifs, this will limit the scanning to a specific list of genes/peaks/kmers/motifs
        If none, all available features (genes/peaks/kmers/motifs) will be scanned.
    cutoff_zscore: `float`, optional (default: 1)
        The z-score cutoff used for Mann–Whitney U test.
    cutoff_logfc: `float`, optional (default: 0.25)
        The log-transformed fold change cutoff between a pair of branches.
    percentile_expr: `int`, optional (default: 95)
        Between 0 and 100. Between 0 and 100. Specify the percentile of marker expression greater than 0 to filter out some extreme marker expressions. 
    n_jobs: `int`, optional (default: 1)
        The number of parallel jobs to run when scaling the marker expressions .
    min_num_cells: `int`, optional (default: 5)
        The minimum number of cells in which markers are expressed.
    use_precomputed: `bool`, optional (default: True)
        If True, the previously computed scaled marker expression will be used
    root: `str`, optional (default: 'S0'): 
        The starting node
    preference: `list`, optional (default: None): 
        The preference of nodes. The branch with speficied nodes are preferred and will be dealt with first. The higher ranks the node have, The earlier the branch with that node will be analyzed. 
        This will help generate the consistent results as shown in subway map and stream plot.

    Returns
    -------
    updates `adata` with the following fields.
    scaled_marker_expr: `list` (`adata.uns['scaled_marker_expr']`)
        Scaled marker expression for marker marker detection.    
    de_markers_greater: `dict` (`adata.uns['de_markers_greater']`)
        DE(differentially expressed) markers for each pair of branches deteced by STREAM. 
        Store the markers that have higher expression on the former part of one branch pair.
    de_markers_less: `dict` (`adata.uns['de_markers_less']`)
        DE(differentially expressed) markers for each pair of branches deteced by STREAM. 
        Store the markers that have higher expression on the latter part of one branch pair.
    """

    file_path = os.path.join(adata.uns['workdir'],'de_markers')
    if(not os.path.exists(file_path)):
        os.makedirs(file_path)    

    if(marker_list is None):
        print('Scanning all features ...')
        marker_list = adata.var_names.tolist()
    else:
        print('Scanning the specified marker list ...')
        ###remove duplicate keys
        marker_list = list(dict.fromkeys(marker_list)) 
        for marker in marker_list:
            if(marker not in adata.var_names):
                raise ValueError("could not find '%s' in `adata.var_names`"  % (marker))

    flat_tree = adata.uns['flat_tree']
    dict_node_state = nx.get_node_attributes(flat_tree,'label')
    df_marker_detection = adata.obs.copy()
    df_marker_detection.rename(columns={"label": "CELL_LABEL", "branch_lam": "lam"},inplace = True)

    if(use_precomputed and ('scaled_marker_expr' in adata.uns_keys())):
        print('Importing precomputed scaled marker expression matrix ...')
        results = adata.uns['scaled_marker_expr']  
        df_results = pd.DataFrame(results).T 
        if(all(np.isin(marker_list,df_results.columns.tolist()))):
            input_markers_expressed = marker_list
            df_scaled_marker_expr = df_results[input_markers_expressed]
        else:
            raise ValueError("Not all markers in `marker_list` can be found in precomputed scaled marker expression matrix. Please set `use_precomputed=False`")   
    else:
        input_markers = marker_list
        df_sc = pd.DataFrame(index= adata.obs_names.tolist(),
                             data = adata[:,input_markers].X,
                             columns=input_markers)
        #exclude markers that are expressed in fewer than min_num_cells cells
        #min_num_cells = max(5,int(round(df_marker_detection.shape[0]*0.001)))
        print("Filtering out markers that are expressed in less than " + str(min_num_cells) + " cells ...")
        input_markers_expressed = np.array(input_markers)[np.where((df_sc[input_markers]>0).sum(axis=0)>min_num_cells)[0]].tolist()
        df_marker_detection[input_markers_expressed] = df_sc[input_markers_expressed].copy()

        print(str(n_jobs)+' cpus are being used ...')
        params = [(df_marker_detection,x,percentile_expr) for x in input_markers_expressed]
        pool = multiprocessing.Pool(processes=n_jobs)
        results = pool.map(scale_marker_expr,params)
        pool.close()
        adata.uns['scaled_marker_expr'] = results
        df_scaled_marker_expr = pd.DataFrame(results).T

    print(str(len(input_markers_expressed)) + ' markers are being scanned ...')
    df_marker_detection[input_markers_expressed] = df_scaled_marker_expr 

    #### DE (Differentially expressed markers) between sub-branches
    dict_de_greater = dict()
    dict_de_less = dict()
    dict_label_node = {value: key for key,value in nx.get_node_attributes(flat_tree,'label').items()}
    if(preference!=None):
        preference_nodes = [dict_label_node[x] for x in preference]
    else:
        preference_nodes = None
    root_node = dict_label_node[root]
    bfs_edges = bfs_edges_modified(flat_tree,root_node,preference=preference_nodes)
    pairs_branches = list()
#     all_branches = np.unique(df_marker_detection['branch_id']).tolist()
    for node_i in dict_node_state.keys():
        neighbor_branches = [x for x in bfs_edges if node_i in x]
        if(len(neighbor_branches)>1):
            pairs_branches += list(itertools.combinations(neighbor_branches,r=2))
    for pair_i in pairs_branches:
        pair_i_alias = ((dict_node_state[pair_i[0][0]],dict_node_state[pair_i[0][1]]),(dict_node_state[pair_i[1][0]],dict_node_state[pair_i[1][1]]))
        if(pair_i[0] in nx.get_edge_attributes(flat_tree,'id').values()):
            df_cells_sub1 = df_marker_detection[df_marker_detection.branch_id==pair_i[0]]
        else:
            df_cells_sub1 = df_marker_detection[df_marker_detection.branch_id==(pair_i[0][1],pair_i[0][0])]
        if(pair_i[1] in nx.get_edge_attributes(flat_tree,'id').values()):
            df_cells_sub2 = df_marker_detection[df_marker_detection.branch_id==pair_i[1]]  
        else:
            df_cells_sub2 = df_marker_detection[df_marker_detection.branch_id==(pair_i[1][1],pair_i[1][0])]
        #only use Mann-Whitney U test when the number of observation in each sample is > 20
        if(df_cells_sub1.shape[0]>20 and df_cells_sub2.shape[0]>20):
            df_de_pval_qval = pd.DataFrame(columns = ['z_score','U','logfc','mean_up','mean_down','pval','qval'],dtype=float)
            for markername in input_markers_expressed:
                sub1_values = df_cells_sub1.loc[:,markername].tolist()
                sub2_values = df_cells_sub2.loc[:,markername].tolist()
                diff_mean = abs(np.mean(sub1_values) - np.mean(sub2_values))
                if(diff_mean>0):
                    logfc = np.log2(max(np.mean(sub1_values),np.mean(sub2_values))/(min(np.mean(sub1_values),np.mean(sub2_values))+diff_mean/1000.0))
                else:
                    logfc = 0
                if(logfc>cutoff_logfc):
                    df_de_pval_qval.loc[markername] = np.nan
                    df_de_pval_qval.loc[markername,['U','pval']] = mannwhitneyu(sub1_values,sub2_values,alternative='two-sided')
                    df_de_pval_qval.loc[markername,'logfc'] = logfc
                    df_de_pval_qval.loc[markername,'mean_up'] = np.mean(sub1_values)
                    df_de_pval_qval.loc[markername,'mean_down'] = np.mean(sub2_values)
                    sub1_sub2_values = sub1_values + sub2_values
                    len_sub1 = len(sub1_values)
                    len_sub2 = len(sub2_values)
                    len_sub = len_sub1+len_sub2
                    ranks_values = stats.rankdata(sub1_sub2_values)
                    ranks_id, ranks_count = np.unique(ranks_values,return_counts=True)
                    sum_k = 0
                    for i,t in enumerate(ranks_count):
                        sum_k = sum_k + (t**3-t)/float(len_sub*(len_sub-1))
                    mu_U = (len_sub1*len_sub2)/2.0
                    sigma_U = np.sqrt(len_sub1*len_sub2*(len_sub+1-sum_k)/12.0)
                    df_de_pval_qval.loc[markername,'z_score'] = (df_de_pval_qval.loc[markername,'U']-mu_U)/sigma_U
            if(df_de_pval_qval.shape[0]==0):
                print('No DE markers are detected between branches ' + pair_i_alias[0][0]+'_'+pair_i_alias[0][1]+\
                      ' and '+pair_i_alias[1][0]+'_'+pair_i_alias[1][1])
            else:
                p_values = df_de_pval_qval['pval']
                q_values = multipletests(p_values, method='fdr_bh')[1]
                df_de_pval_qval['qval'] = q_values
                dict_de_greater[pair_i_alias] = df_de_pval_qval[(abs(df_de_pval_qval['z_score'])>cutoff_zscore)&
                                                          (df_de_pval_qval['z_score']>0)].sort_values(['z_score'],ascending=False)
                dict_de_greater[pair_i_alias].to_csv(os.path.join(file_path,'de_markers_greater_'+pair_i_alias[0][0]+'_'+pair_i_alias[0][1] + ' and '\
                                        + pair_i_alias[1][0]+'_'+pair_i_alias[1][1] + '.tsv'),sep = '\t',index = True)
                dict_de_less[pair_i_alias] = df_de_pval_qval[(abs(df_de_pval_qval['z_score'])>cutoff_zscore)&
                                                       (df_de_pval_qval['z_score']<0)].sort_values(['z_score'])
                dict_de_less[pair_i_alias].to_csv(os.path.join(file_path,'de_markers_less_'+pair_i_alias[0][0]+'_'+pair_i_alias[0][1] + ' and '\
                                     + pair_i_alias[1][0]+'_'+pair_i_alias[1][1] + '.tsv'),sep = '\t',index = True)   
        else:
            print('There are not sufficient cells (should be greater than 20) between branches '+\
                  pair_i_alias[0][0]+'_'+pair_i_alias[0][1] +' and '+\
                  pair_i_alias[1][0]+'_'+pair_i_alias[1][1]+ '. fold_change is calculated')
            df_de_pval_qval = pd.DataFrame(columns = ['logfc','mean_up','mean_down'],dtype=float)
            for markername in input_markers_expressed:
                sub1_values = df_cells_sub1.loc[:,markername].tolist()
                sub2_values = df_cells_sub2.loc[:,markername].tolist()
                diff_mean = abs(np.mean(sub1_values) - np.mean(sub2_values))
                if(diff_mean>0):
                    logfc = np.log2(max(np.mean(sub1_values),np.mean(sub2_values))/(min(np.mean(sub1_values),np.mean(sub2_values))+diff_mean/1000.0))
                else:
                    logfc = 0
                if(logfc>cutoff_logfc):
                    df_de_pval_qval.loc[markername] = np.nan
                    # #make sure the largest fold change is 5
                    # df_de_pval_qval.loc[markername,'fold_change'] = np.log2((np.mean(sub1_values)+1/24.0)/(np.mean(sub2_values)+1/24.0))
                    df_de_pval_qval.loc[markername,'logfc'] = logfc
                    df_de_pval_qval.loc[markername,'mean_up'] = np.mean(sub1_values)
                    df_de_pval_qval.loc[markername,'mean_down'] = np.mean(sub2_values)
            if(df_de_pval_qval.shape[0]==0):
                print('No DE markers are detected between branches ' + pair_i_alias[0][0]+'_'+pair_i_alias[0][1]+\
                      ' and '+pair_i_alias[1][0]+'_'+pair_i_alias[1][1])
            else:
                dict_de_greater[pair_i_alias] = df_de_pval_qval[(abs(df_de_pval_qval['logfc'])>cutoff_logfc)&
                                                          (df_de_pval_qval['logfc']>0)].sort_values(['logfc'],ascending=False)
                dict_de_greater[pair_i_alias].to_csv(os.path.join(file_path,'de_markers_greater_'+pair_i_alias[0][0]+'_'+pair_i_alias[0][1] + ' and '\
                                        + pair_i_alias[1][0]+'_'+pair_i_alias[1][1] + '.tsv'),sep = '\t',index = True)                
                dict_de_less[pair_i_alias] = df_de_pval_qval[(abs(df_de_pval_qval['logfc'])>cutoff_logfc)&
                                                       (df_de_pval_qval['logfc']<0)].sort_values(['logfc'])
                dict_de_less[pair_i_alias].to_csv(os.path.join(file_path,'de_markers_less_'+pair_i_alias[0][0]+'_'+pair_i_alias[0][1] + ' and '\
                                     + pair_i_alias[1][0]+'_'+pair_i_alias[1][1] + '.tsv'),sep = '\t',index = True)   
    adata.uns['de_markers_greater'] = dict_de_greater
    adata.uns['de_markers_less'] = dict_de_less


def plot_de_markers(adata,num_markers = 15,cutoff_zscore=1,cutoff_logfc = 0.25,
                  save_fig=False,fig_path=None,fig_size=(12,8)):
    if(fig_path is None):
        fig_path = os.path.join(adata.uns['workdir'],'de_markers')
    if(not os.path.exists(fig_path)):
        os.makedirs(fig_path)  

    dict_de_greater = adata.uns['de_markers_greater']
    dict_de_less = adata.uns['de_markers_less']
    flat_tree = adata.uns['flat_tree']
    dict_node_state = nx.get_node_attributes(flat_tree,'label')
    colors = sns.color_palette("Set1", n_colors=8, desat=0.8)
    for sub_edges_i in dict_de_greater.keys():
        fig = plt.figure(figsize=fig_size)
        gs = gridspec.GridSpec(2,1)
        ax = fig.add_subplot(gs[0],adjustable='box')
        if('z_score' in dict_de_greater[sub_edges_i].columns):
            if(not dict_de_greater[sub_edges_i].empty):
                val_greater = dict_de_greater[sub_edges_i].iloc[:num_markers,:]['z_score'].values  # the bar lengths
                pos_greater = np.arange(dict_de_greater[sub_edges_i].iloc[:num_markers,:].shape[0])-1    # the bar centers on the y axis
            else:
                val_greater = np.repeat(0,num_markers)
                pos_greater = np.arange(num_markers)-1
            ax.bar(pos_greater,val_greater, align='center',color = colors[0])
            ax.plot([pos_greater[0]-1,pos_greater[-1]+1], [cutoff_zscore, cutoff_zscore], "k--")
            q_vals = dict_de_greater[sub_edges_i].iloc[:num_markers,:]['qval'].values
            for i, q in enumerate(q_vals):
                alignment = {'horizontalalignment': 'center', 'verticalalignment': 'bottom'}
                ax.text(pos_greater[i], val_greater[i]+.1, \
                        "{:.2E}".format(Decimal(str(q))),color='black',size=0.6*mpl.rcParams['font.size'],
                        **alignment)
            plt.xticks(pos_greater,dict_de_greater[sub_edges_i].index,rotation=-45)
            ax.set_ylabel('z score')
            ax.set_title('DE markers between branches ' + sub_edges_i[0][0]+'_'+sub_edges_i[0][1] + ' and ' + \
                         sub_edges_i[1][0]+'_'+sub_edges_i[1][1]) 
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False) 
            ax.locator_params(axis='y',nbins=5)
        
            ax1 = fig.add_subplot(gs[1], adjustable='box')
            if(not dict_de_less[sub_edges_i].empty):
                val_less = dict_de_less[sub_edges_i].iloc[:num_markers,:]['z_score'].values  # the bar lengths
                pos_less = np.arange(dict_de_less[sub_edges_i].iloc[:num_markers,:].shape[0])-1    # the bar centers on the y axis
            else:
                val_less = np.repeat(0,num_markers)
                pos_less = np.arange(num_markers)-1
            ax1.bar(pos_less,val_less, align='center',color = colors[1])
            ax1.plot([pos_less[0]-1,pos_less[-1]+1], [-cutoff_zscore, -cutoff_zscore], "k--")
            q_vals = dict_de_less[sub_edges_i].iloc[:num_markers,:]['qval'].values
            for i, q in enumerate(q_vals):
                alignment = {'horizontalalignment': 'center', 'verticalalignment': 'top'}
                ax1.text(pos_less[i], val_less[i]-.1, \
                        "{:.2E}".format(Decimal(str(q))),color='black',size=0.6*mpl.rcParams['font.size'],
                         **alignment)
            plt.xticks(pos_less,dict_de_less[sub_edges_i].index)
            ax1.set_xticklabels(dict_de_less[sub_edges_i].index,rotation=-45)
            ax1.set_ylabel('z score')

            ax1.spines['right'].set_visible(False)
            ax1.spines['bottom'].set_visible(False) 
            ax1.xaxis.tick_top()
            ax1.locator_params(axis='y',nbins=5)
            
            plt.tight_layout()
            if(save_fig):
                plt.savefig(os.path.join(fig_path,'de_markers_'+sub_edges_i[0][0]+'_'+sub_edges_i[0][1] + ' and '\
                            + sub_edges_i[1][0]+'_'+sub_edges_i[1][1]+'.pdf'),pad_inches=1,bbox_inches='tight')
                plt.close(fig)
        else:
            if(not dict_de_greater[sub_edges_i].empty):
                val_greater = dict_de_greater[sub_edges_i].iloc[:num_markers,:]['logfc'].values  # the bar lengths
                pos_greater = np.arange(dict_DE_greater[sub_edges_i].iloc[:num_markers,:].shape[0])-1    # the bar centers on the y axis
            else:
                val_greater = np.repeat(0,num_markers)
                pos_greater = np.arange(num_markers)-1
            ax.bar(pos_greater,val_greater, align='center',color = colors[0])
            ax.plot([pos_greater[0]-1,pos_greater[-1]+1], [cutoff_logfc, cutoff_logfc], "k--")
            plt.xticks(pos_greater,dict_de_greater[sub_edges_i].index,rotation=-45)
            ax.set_ylabel('log(fold change)')
            ax.set_title('DE markers between branches ' + sub_edges_i[0][0]+'_'+sub_edges_i[0][1] + ' and ' + \
                         sub_edges_i[1][0]+'_'+sub_edges_i[1][1])
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.locator_params(axis='y',nbins=5)
            
            ax1 = fig.add_subplot(gs[1], adjustable='box')
            if(not dict_de_less[sub_edges_i].empty):
                val_less = dict_de_less[sub_edges_i].iloc[:num_markers,:]['logfc'].values  # the bar lengths
                pos_less = np.arange(dict_de_less[sub_edges_i].iloc[:num_markers,:].shape[0])-1    # the bar centers on the y axis
            else:
                val_less = np.repeat(0,num_markers)
                pos_less = np.arange(num_markers)-1
            ax1.bar(pos_less,val_less, align='center',color = colors[1])
            ax1.plot([pos_less[0]-1,pos_less[-1]+1], [-cutoff_logfc, -cutoff_logfc], "k--")
            plt.xticks(pos_less,dict_de_less[sub_edges_i].index)
            ax1.set_xticklabels(dict_de_less[sub_edges_i].index,rotation=-45)
            ax1.set_ylabel('log(fold change)')

            ax1.spines['right'].set_visible(False)
            ax1.spines['bottom'].set_visible(False) 
            ax1.xaxis.tick_top()
            ax1.locator_params(axis='y',nbins=5)
            
            plt.tight_layout()
            if(save_fig):
                plt.savefig(os.path.join(fig_path,'de_markers_'+sub_edges_i[0][0]+'_'+sub_edges_i[0][1] + ' and '\
                            + sub_edges_i[1][0]+'_'+sub_edges_i[1][1]+'.pdf'),pad_inches=1,bbox_inches='tight')
                plt.close(fig) 


def detect_leaf_markers(adata,marker_list=None,cutoff_zscore=1.,cutoff_pvalue=1e-2,percentile_expr=95,n_jobs = 1,min_num_cells=5,
                        use_precomputed=True, root='S0',preference=None):
    """Detect leaf markers for each branch.
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    marker_list: `list`, optional (default: None): 
        A list of candidate markers to be scanned. Instead of scanning all available genes/peaks/kmers/motifs, this will limit the scanning to a specific list of genes/peaks/kmers/motifs
        If none, all available features (genes/peaks/kmers/motifs) will be scanned.
    cutoff_zscore: `float`, optional (default: 1.5)
        The z-score cutoff used for mean values of all leaf branches.
    cutoff_pvalue: `float`, optional (default: 1e-2)
        The p value cutoff used for Kruskal-Wallis H-test and post-hoc pairwise Conover’s test.
    percentile_expr: `int`, optional (default: 95)
        Between 0 and 100. Between 0 and 100. Specify the percentile of marker expression greater than 0 to filter out some extreme marker expressions. 
    n_jobs: `int`, optional (default: 1)
        The number of parallel jobs to run when scaling the marker expressions .
    min_num_cells: `int`, optional (default: 5)
        The minimum number of cells in which markers are expressed.
    use_precomputed: `bool`, optional (default: True)
        If True, the previously computed scaled marker expression will be used
    root: `str`, optional (default: 'S0'): 
        The starting node
    preference: `list`, optional (default: None): 
        The preference of nodes. The branch with speficied nodes are preferred and will be dealt with first. The higher ranks the node have, The earlier the branch with that node will be analyzed. 
        This will help generate the consistent results as shown in subway map and stream plot.

    Returns
    -------
    updates `adata` with the following fields.
    scaled_marker_expr: `list` (`adata.uns['scaled_marker_expr']`)
        Scaled marker expression for marker marker detection.    
    leaf_markers_all: `pandas.core.frame.DataFrame` (`adata.uns['leaf_markers_all']`)
        All the leaf markers from all leaf branches.
    leaf_markers: `dict` (`adata.uns['leaf_markers']`)
        Leaf markers for each branch.
    """

    file_path = os.path.join(adata.uns['workdir'],'leaf_markers')
    if(not os.path.exists(file_path)):
        os.makedirs(file_path)    

    if(marker_list is None):
        print('Scanning all features ...')
        marker_list = adata.var_names.tolist()
    else:
        print('Scanning the specified marker list ...')
        ###remove duplicate keys
        marker_list = list(dict.fromkeys(marker_list)) 
        for marker in marker_list:
            if(marker not in adata.var_names):
                raise ValueError("could not find '%s' in `adata.var_names`"  % (marker))

    flat_tree = adata.uns['flat_tree']
    dict_node_state = nx.get_node_attributes(flat_tree,'label')
    df_marker_detection = adata.obs.copy()
    df_marker_detection.rename(columns={"label": "CELL_LABEL", "branch_lam": "lam"},inplace = True)

    if(use_precomputed and ('scaled_marker_expr' in adata.uns_keys())):
        print('Importing precomputed scaled marker expression matrix ...')
        results = adata.uns['scaled_marker_expr']          
        df_results = pd.DataFrame(results).T 
        if(all(np.isin(marker_list,df_results.columns.tolist()))):
            input_markers_expressed = marker_list
            df_scaled_marker_expr = df_results[input_markers_expressed]
        else:
            raise ValueError("Not all markers in `marker_list` can be found in precomputed scaled marker expression matrix. Please set `use_precomputed=False`")
    else:
        input_markers = marker_list
        df_sc = pd.DataFrame(index= adata.obs_names.tolist(),
                             data = adata[:,input_markers].X,
                             columns=input_markers)
        #exclude markers that are expressed in fewer than min_num_cells cells
        print("Filtering out markers that are expressed in less than " + str(min_num_cells) + " cells ...")
        input_markers_expressed = np.array(input_markers)[np.where((df_sc[input_markers]>0).sum(axis=0)>min_num_cells)[0]].tolist()
        df_marker_detection[input_markers_expressed] = df_sc[input_markers_expressed].copy()

        print(str(n_jobs)+' cpus are being used ...')
        params = [(df_marker_detection,x,percentile_expr) for x in input_markers_expressed]
        pool = multiprocessing.Pool(processes=n_jobs)
        results = pool.map(scale_marker_expr,params)
        pool.close()
        adata.uns['scaled_marker_expr'] = results
        df_scaled_marker_expr = pd.DataFrame(results).T

    print(str(len(input_markers_expressed)) + ' markers are being scanned ...')
    df_marker_detection[input_markers_expressed] = df_scaled_marker_expr    

    #### find marker markers that are specific to one leaf branch
    dict_label_node = {value: key for key,value in nx.get_node_attributes(flat_tree,'label').items()}
    if(preference!=None):
        preference_nodes = [dict_label_node[x] for x in preference]
    else:
        preference_nodes = None
    root_node = dict_label_node[root]
    bfs_edges = bfs_edges_modified(flat_tree,root_node,preference=preference_nodes)
    leaves = [k for k,v in flat_tree.degree() if v==1]
    leaf_edges = [x for x in bfs_edges if (x[0] in leaves) or (x[1] in leaves)]    

    df_marker_detection['bfs_edges'] = df_marker_detection['branch_id']
    df_marker_detection.astype('object')
    for x in df_marker_detection['branch_id'].unique():
        id_ = df_marker_detection[df_marker_detection['branch_id']==x].index
        if x not in bfs_edges:
            df_marker_detection.loc[id_,'bfs_edges'] =pd.Series(index=id_,data=[(x[1],x[0])]*len(id_))
    
    df_leaf_markers = pd.DataFrame(columns=['zscore','H_statistic','H_pvalue']+leaf_edges)
    for marker in input_markers_expressed:
        meann_values = df_marker_detection[['bfs_edges',marker]].groupby(by = 'bfs_edges')[marker].mean()
        br_values = df_marker_detection[['bfs_edges',marker]].groupby(by = 'bfs_edges')[marker].apply(list)
        leaf_mean_values = meann_values[leaf_edges]
        leaf_mean_values.sort_values(inplace=True)
        leaf_br_values = br_values[leaf_edges]
        if(leaf_mean_values.shape[0]<=2):
            print('There are not enough leaf branches')
        else:
            zscores = stats.zscore(leaf_mean_values)
            if(abs(zscores)[abs(zscores)>cutoff_zscore].shape[0]>=1):
                if(any(zscores>cutoff_zscore)):
                    cand_br = leaf_mean_values.index[-1]
                    cand_zscore = zscores[-1]
                else:
                    cand_br = leaf_mean_values.index[0]
                    cand_zscore = zscores[0]
                list_br_values = [leaf_br_values[x] for x in leaf_edges]
                kurskal_statistic,kurskal_pvalue = stats.kruskal(*list_br_values)
                if(kurskal_pvalue<cutoff_pvalue):  
                    df_conover_pvalues= posthoc_conover(df_marker_detection[[x in leaf_edges for x in df_marker_detection['bfs_edges']]], 
                                                       val_col=marker, group_col='bfs_edges', p_adjust = 'fdr_bh')
                    cand_conover_pvalues = df_conover_pvalues[~df_conover_pvalues.columns.isin([cand_br])][cand_br]
                    if(all(cand_conover_pvalues < cutoff_pvalue)):
                        df_leaf_markers.loc[marker,:] = 1.0
                        df_leaf_markers.loc[marker,['zscore','H_statistic','H_pvalue']] = [cand_zscore,kurskal_statistic,kurskal_pvalue]
                        df_leaf_markers.loc[marker,cand_conover_pvalues.index] = cand_conover_pvalues
    df_leaf_markers.rename(columns={x:dict_node_state[x[0]]+dict_node_state[x[1]]+'_pvalue' for x in leaf_edges},inplace=True)
    df_leaf_markers.sort_values(by=['H_pvalue','zscore'],ascending=[True,False],inplace=True)
    df_leaf_markers.to_csv(os.path.join(file_path,'leaf_markers.tsv'),sep = '\t',index = True)
    dict_leaf_markers = dict()
    for x in leaf_edges:
        x_alias = (dict_node_state[x[0]],dict_node_state[x[1]])
        dict_leaf_markers[x_alias] = df_leaf_markers[df_leaf_markers[x_alias[0]+x_alias[1]+'_pvalue']==1.0]
        dict_leaf_markers[x_alias].to_csv(os.path.join(file_path,'leaf_markers'+x_alias[0]+'_'+x_alias[1] + '.tsv'),sep = '\t',index = True)
    adata.uns['leaf_markers_all'] = df_leaf_markers
    adata.uns['leaf_markers'] = dict_leaf_markers


def detect_markers(adata,marker_list=None,ident='label',cutoff_zscore=1.,cutoff_pvalue=1e-2,percentile_expr=95,n_jobs = 1,min_num_cells=5,
                use_precomputed=True):
    """Detect markers (highly expressed or suppressed) for the specified ident.
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    marker_list: `list`, optional (default: None): 
        A list of candidate markers to be scanned. Instead of scanning all available genes/peaks/kmers/motifs, this will limit the scanning to a specific list of genes/peaks/kmers/motifs
        If none, all available features (genes/peaks/kmers/motifs) will be scanned.
    cutoff_zscore: `float`, optional (default: 1.5)
        The z-score cutoff used for mean values of all leaf branches.
    cutoff_pvalue: `float`, optional (default: 1e-2)
        The p value cutoff used for Kruskal-Wallis H-test and post-hoc pairwise Conover’s test.
    percentile_expr: `int`, optional (default: 95)
        Between 0 and 100. Between 0 and 100. Specify the percentile of marker expression greater than 0 to filter out some extreme marker expressions. 
    n_jobs: `int`, optional (default: all available cpus)
        The number of parallel jobs to run when scaling the marker expressions .
    min_num_cells: `int`, optional (default: 5)
        The minimum number of cells in which markers are expressed.
    use_precomputed: `bool`, optional (default: True)
        If True, the previously computed scaled marker expression will be used

    Returns
    -------
    updates `adata` with the following fields.
    scaled_marker_expr: `list` (`adata.uns['scaled_marker_expr']`)
        Scaled marker expression for marker marker detection.    
    markers_ident_all: `pandas.core.frame.DataFrame` (`adata.uns['markers_ident_all']`)
        All markers for all ident labels.
    markers_ident: `dict` (`adata.uns['markers_']`)
        Markers for each ident label.
    """    

    file_path = os.path.join(adata.uns['workdir'],'markers_found')
    if(not os.path.exists(file_path)):
        os.makedirs(file_path)  

    if(marker_list is None):
        print('Scanning all features ...')
        marker_list = adata.var_names.tolist()
    else:
        print('Scanning the specified marker list ...')
        ###remove duplicate keys
        marker_list = list(dict.fromkeys(marker_list)) 
        for marker in marker_list:
            if(marker not in adata.var_names):
                raise ValueError("could not find '%s' in `adata.var_names`"  % (marker))
        
    if ident not in adata.obs.columns:
        raise ValueError("'%s' does not exist in `adata.obs`" %(ident))

    if(use_precomputed and ('scaled_marker_expr' in adata.uns_keys())):
        print('Importing precomputed scaled marker expression matrix ...')
        results = adata.uns['scaled_marker_expr']
        df_results = pd.DataFrame(results).T 
        if(all(np.isin(marker_list,df_results.columns.tolist()))):
            input_markers_expressed = marker_list
            df_scaled_marker_expr = df_results[input_markers_expressed]
        else:
            raise ValueError("Not all markers in `marker_list` can be found in precomputed scaled marker expression matrix. Please set `use_precomputed=False`")                  
    else:
        input_markers = marker_list
        df_sc = pd.DataFrame(index= adata.obs_names.tolist(),
                            data = adata[:,input_markers].X,
                            columns=input_markers)
        #exclude markers that are expressed in fewer than min_num_cells cells
        print("Filtering out markers that are expressed in less than " + str(min_num_cells) + " cells ...")
        input_markers_expressed = np.array(input_markers)[np.where((df_sc[input_markers]>0).sum(axis=0)>min_num_cells)[0]].tolist()
        df_sc_filtered = df_sc[input_markers_expressed].copy()

        print(str(n_jobs)+' cpus are being used ...')
        params = [(df_sc_filtered,x,percentile_expr) for x in input_markers_expressed]
        pool = multiprocessing.Pool(processes=n_jobs)
        results = pool.map(scale_marker_expr,params)
        pool.close()
        adata.uns['scaled_marker_expr'] = results
        df_scaled_marker_expr = pd.DataFrame(results).T

    print(str(len(input_markers_expressed)) + ' markers are being scanned ...')
    df_input = df_scaled_marker_expr 
    df_input[ident] = adata.obs[ident]
    
    uniq_ident = np.unique(df_input[ident]).tolist()

    df_markers = pd.DataFrame(columns=['zscore','H_statistic','H_pvalue']+uniq_ident)
    for marker in input_markers_expressed:
        mean_values = df_input[[ident,marker]].groupby(by = ident)[marker].mean()
        mean_values.sort_values(inplace=True)
        list_marker_expr = df_input[[ident,marker]].groupby(by = ident)[marker].apply(list)
        if(mean_values.shape[0]<2):
            print('At least two distinct' + ident + 'are required')
        else:
            zscores = stats.zscore(mean_values)
            if(abs(zscores)[abs(zscores)>cutoff_zscore].shape[0]>=1):
                if(any(zscores>cutoff_zscore)):
                    cand_ident = mean_values.index[-1]
                    cand_zscore = zscores[-1]
                else:
                    cand_ident = mean_values.index[0]
                    cand_zscore = zscores[0]
                kurskal_statistic,kurskal_pvalue = stats.kruskal(*list_marker_expr)
                if(kurskal_pvalue<cutoff_pvalue):  
                    df_conover_pvalues= posthoc_conover(df_input, 
                                                       val_col=marker, group_col=ident, p_adjust = 'fdr_bh')
                    cand_conover_pvalues = df_conover_pvalues[cand_ident]
                    if(all(cand_conover_pvalues < cutoff_pvalue)):
    #                     df_markers.loc[marker,:] = 1.0
                        df_markers.loc[marker,['zscore','H_statistic','H_pvalue']] = [cand_zscore,kurskal_statistic,kurskal_pvalue]
                        df_markers.loc[marker,cand_conover_pvalues.index] = cand_conover_pvalues
    df_markers.sort_values(by=['H_pvalue','zscore'],ascending=[True,False],inplace=True)    
    df_markers.to_csv(os.path.join(file_path,'markers_'+ident+'.tsv'),sep = '\t',index = True)
    dict_markers = dict()
    for x in uniq_ident:
        dict_markers[x] = df_markers[df_markers[x]==-1].loc[:,df_markers.columns!=x]
        dict_markers[x].to_csv(os.path.join(file_path,'markers_'+ident+'_'+x+'.tsv'),sep = '\t',index = True)
    adata.uns['markers_'+ident+'_all'] = df_markers
    adata.uns['markers_'+ident] = dict_markers 
def map_new_data(adata_ref,adata_new,use_radius=False):
    """ Map new data to the inferred trajectories
    
    Parameters
    ----------
    adata_ref: reference AnnData
        Annotated data matrix.
    adata_new: new AnnData
        Annotated data matrix for new data (to be mapped).
    use_radius: `bool`, optional (default: False)
        Only valid when `method = 'mlle'`. If True, when searching for the neighbors for each cell in MLLE space, STREAM uses a fixed radius instead of a fixed number of cells.
    first_pc: `bool`, optional (default: False)
        Only valid when `feature='top_pcs'` If True, the first principal component will be included
    top_pcs_feature: `str`, optional (default: None)
        Choose from {{'var_genes'}}
        Only valid when `feature='top_pcs'`.
        Features used for pricipal component analysis
        If None, all the genes will be used.
        IF 'var_genes', the most variable genes obtained from select_variable_genes() will be used.        
    Returns
    -------  
    
    Combined AnnData object. Only the shared obs_keys and var_keys will be preserved
    updates `adata` with the following fields.
    batch: `pandas.Series` (`adata.obs['batch']`,dtype `str`)
        The annotation of each cell. It consists of 'ref' for reference cells and 'new' for the cells to map    

    updates `adata_new` with the following fields.(depending on the `feature` or `method`)
    var_genes: `numpy.ndarray` (`adata_new.obsm['var_genes']`)
        Store #observations × #var_genes data matrix used mapping.
    var_genes: `pandas.core.indexes.base.Index` (`adata_new.uns['var_genes']`)
        The selected variable gene names.
    top_pcs: `numpy.ndarray` (`adata_new.obsm['top_pcs']`)
        Store #observations × n_pc data matrix used for subsequent dimension reduction.
    X_dr : `numpy.ndarray` (`adata_new.obsm['X_dr']`)
        A #observations × n_components data matrix after dimension reduction.
    X_mlle : `numpy.ndarray` (`adata_new.obsm['X_mlle']`)
        Store #observations × n_components data matrix after mlle.
    X_umap : `numpy.ndarray` (`adata_new.obsm['X_umap']`)
        Store #observations × n_components data matrix after umap.
    X_pca : `numpy.ndarray` (`adata_new.obsm['X_pca']`)
        Store #observations × n_components data matrix after pca.
    X_vis_umap: `numpy.ndarray` (`adata.obsm['X_vis_umap']`)
        Store #observations × 2 umap data matrix. 
    X_vis_tsne: `numpy.ndarray` (`adata.obsm['X_vis_tsne']`)
        Store #observations × 2 tsne data matrix.  
    """

    feature = adata_ref.uns['params']['dimension_reduction']['feature']
    method = adata_ref.uns['params']['dimension_reduction']['method']
    assert (method in ['mlle','umap','pca']),"'%s' is not supported. Please use one of 'mlle','umap', and 'pca' for dimension reduction." %method
    if(feature == 'var_genes'):
        print('Top variable genes are being used for mapping ...')
        if('var_genes' not in adata_ref.uns_keys()):
            raise ValueError("variable genes are not selected yet in `adata_ref`")
        for x in adata_ref.uns['var_genes']:
            if(x not in adata_new.var_names):
                raise ValueError("variable gene '%s' does not exist in `adata_new.var_names`"  % (x))
        adata_new.uns['var_genes'] = adata_ref.uns['var_genes'].copy()
        adata_new.obsm['var_genes'] = adata_new[:,adata_new.uns['var_genes']].X.copy()
        input_data = adata_new.obsm['var_genes']
    if(feature == 'all'):
        print('All genes are being used for mapping ...')
        if(not set(adata_ref.var_names) <= set(adata_new.var_names)):
            raise ValueError("`adata_new.var_names` does not contain all the genes in `adata_ref.var_names`")
        input_data = adata_new[:,adata_ref.var.index].X
    if(feature == 'top_pcs'):
        print('Top principal components are being used for mapping ...')
        if('top_pcs' not in adata_ref.uns_keys()):
            raise ValueError("top principal components are not selected yet in `adata_ref`")
        trans = adata_ref.uns['top_pcs']
        top_pcs_feature = adata_ref.uns['params']['select_top_principal_components']['feature']
        if(top_pcs_feature == 'var_genes'):
            if('var_genes' not in adata_ref.uns_keys()):
                raise ValueError("variable genes are not selected yet in `adata_ref`")
            adata_new.uns['var_genes'] = adata_ref.uns['var_genes'].copy()
            adata_new.obsm['var_genes'] = adata_new[:,adata_new.uns['var_genes']].X.copy()
            X_pca = trans.transform(adata_new.obsm['var_genes']) 
        else:
            X_pca = trans.transform(adata_new[:,adata_ref.var.index].X) 
        n_pc = adata_ref.obsm['top_pcs'].shape[1]
        first_pc = adata_ref.uns['params']['select_top_principal_components']['first_pc']
        if(first_pc):
            adata_new.obsm['top_pcs'] = X_pca[:,0:(n_pc)]
        else:
            #discard the first Principal Component
            adata_new.obsm['top_pcs'] = X_pca[:,1:(n_pc+1)]
        input_data = adata_new.obsm['top_pcs']
    adata_new.uns['epg'] = adata_ref.uns['epg'].copy()
    adata_new.uns['flat_tree'] = adata_ref.uns['flat_tree'].copy() 

    print("method '%s' is being used for mapping ..." %method)
    if(method == 'mlle'):
        if('trans_mlle' in adata_ref.uns_keys()):
            trans = adata_ref.uns['trans_mlle']
            if(use_radius):
                dist_nb = trans.nbrs_.kneighbors(input_data, n_neighbors=trans.n_neighbors,return_distance=True)[0]
                ind = trans.nbrs_.radius_neighbors(input_data, radius = dist_nb.max(),return_distance=False)    
                new_X_mlle = np.empty((input_data.shape[0], trans.n_components))
                for i in range(input_data.shape[0]):
                    weights = barycenter_weights_modified(input_data[i], trans.nbrs_._fit_X[ind[i]],reg=trans.reg)
                    new_X_mlle[i] = np.dot(trans.embedding_[ind[i]].T, weights) 
                adata_new.obsm['X_mlle_mapping'] = new_X_mlle              
            else:
                adata_new.obsm['X_mlle_mapping'] = trans.transform(input_data)
            adata_new.obsm['X_dr'] = adata_new.obsm['X_mlle_mapping'].copy()
        else:
            raise Exception("Please run 'st.dimension_reduction()' using 'mlle' first.")  
    if(method == 'umap'):
        if('trans_umap' in adata_ref.uns_keys()):
            trans = adata_ref.uns['trans_umap']
            adata_new.obsm['X_umap_mapping'] = trans.transform(input_data)
            adata_new.obsm['X_dr'] = adata_new.obsm['X_umap_mapping'].copy()
        else:
            raise Exception("Please run 'st.dimension_reduction()' using 'umap' first.")  
    if(method == 'pca'):
        if('trans_pca' in adata_ref.uns_keys()):            
            trans = adata_ref.uns['trans_pca']
            adata_new.obsm['X_pca_mapping'] = trans.transform(input_data)
            adata_new.obsm['X_dr'] = adata_new.obsm['X_pca_mapping'].copy()
        else:
            raise Exception("Please run 'st.dimension_reduction()' using 'pca' first.")  

    if('plot_visualization_2D' in adata_ref.uns['params'].keys()):
        print('Visualizing new cells on 2D plane ...')
        vis_method = adata_ref.uns['params']['plot_visualization_2D']['method']
        vis_trans = 'vis_trans_' + vis_method
        trans = adata_ref.uns[vis_trans]
        adata_new.obsm['X_vis_'+vis_method] = trans.transform(adata_new.obsm['X_dr'])
        if(adata_ref.uns['params']['epg']['use_vis']):
            print("Using the manifold from `plot_visualization_2D()` ")
            adata_new.obsm['X_dr'] = adata_new.obsm['X_vis_'+vis_method]
            
    project_cells_to_epg(adata_new)
    calculate_pseudotime(adata_new)
    adata_combined = adata_ref.concatenate(adata_new,batch_categories=['ref','new'])
    shared_obs_key = [x for x in adata_new.obs_keys() if x in adata_ref.obs_keys()]
    shared_var_key = [x for x in adata_new.var_keys() if x in adata_ref.var_keys()]
    adata_combined.obs = adata_combined.obs[shared_obs_key+['batch']]
    adata_combined.var = adata_combined.var[shared_var_key]
    adata_combined.uns['workdir'] = adata_new.uns['workdir']
    for key in adata_ref.uns_keys():
        if key in ['var_genes', 'epg', 'flat_tree','vis_trans_tsne','vis_trans_umap']:
            adata_combined.uns[key] = adata_ref.uns[key]
        if(key.split('_')[-1]=='color'):
            if(key in adata_new.uns_keys()):
                adata_combined.uns[key] = {**adata_ref.uns[key],**adata_new.uns[key]}
    return(adata_combined)
def save_vr_report(adata,ann_list=None,gene_list=None,file_name='stream_vr_report'):
    """save stream report for single cell VR website http://www.singlecellvr.com/
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    ann_list: `list`, optional (default: None): 
        A list of cell annotation keys. If None, only 'label' will be used.     
    gene_list: `list`, optional (default: None): 
        A list of genes to be displayed
    file_name: `str`, optional (default: 'stream_vr_report')
        Ouput Zip file name.
    
    Returns
    -------
    None

    """ 
    assert (adata.obsm['X_dr'].shape[1]>=3),\
    '''The embedding space should have at least three dimensions. 
    please set `n_component = 3` in `st.dimension_reduction()`'''
    
    if(ann_list is None):
        ann_list = ['label']
    ###remove duplicate keys
    ann_list = list(dict.fromkeys(ann_list)) 
    for ann in ann_list:
        if(ann not in adata.obs.columns):
            raise ValueError("could not find '%s' in `adata.var_names`"  % (ann))
            
    if(gene_list is not None):
        ###remove duplicate keys
        gene_list = list(dict.fromkeys(gene_list)) 
        for gene in gene_list:
            if(gene not in adata.var_names):
                raise ValueError("could not find '%s' in `adata.var_names`"  % (gene))
                
    try:
        file_path = os.path.join(adata.uns['workdir'],file_name)
        if(not os.path.exists(file_path)):
                os.makedirs(file_path)

        flat_tree = adata.uns['flat_tree']
        epg = adata.uns['epg']
        epg_node_pos = nx.get_node_attributes(epg,'pos')
        ft_node_label = nx.get_node_attributes(flat_tree,'label')
        ft_node_pos = nx.get_node_attributes(flat_tree,'pos')
                
        ## output coordinates of stream graph
        list_curves = []
        for edge_i in flat_tree.edges():
            branch_i_pos = np.array([epg_node_pos[i] for i in flat_tree.edges[edge_i]['nodes']])
            df_coord_curve_i = pd.DataFrame(branch_i_pos)
            dict_coord_curves = dict()
            dict_coord_curves['branch_id'] = ft_node_label[edge_i[0]] + '_' + ft_node_label[edge_i[1]]
            dict_coord_curves['xyz'] = [{'x':df_coord_curve_i.iloc[j,0],
                                         'y':df_coord_curve_i.iloc[j,1],
                                         'z':df_coord_curve_i.iloc[j,2]} for j in range(df_coord_curve_i.shape[0])]
            list_curves.append(dict_coord_curves)
        with open(os.path.join(file_path,'graph_paths.json'), 'w') as f:
            json.dump(list_curves, f)     
            
        ## output topology of stream graph
        dict_nodes = dict()
        list_edges = []
        for node_i in flat_tree.nodes():
            dict_nodes_i = dict()
            dict_nodes_i['node_name'] = ft_node_label[node_i]
            dict_nodes_i['xyz'] = {'x':ft_node_pos[node_i][0],
                                   'y':ft_node_pos[node_i][1],
                                   'z':ft_node_pos[node_i][2]}
            dict_nodes[ft_node_label[node_i]] = dict_nodes_i
        for edge_i in flat_tree.edges():
            dict_edges = dict()
            dict_edges['nodes'] = [ft_node_label[edge_i[0]],ft_node_label[edge_i[1]]]
            dict_edges['weight'] = 1
            list_edges.append(dict_edges)
        with open(os.path.join(file_path,'graph_nodes.json'), 'w') as f:
            json.dump(dict_nodes, f)
        with open(os.path.join(file_path,'graph_edges.json'), 'w') as f:
            json.dump(list_edges, f)

        print('STREAM: graph finished!')
   
        ## output coordinates of cells
        list_cells = []
        for i in range(adata.shape[0]):
            dict_coord_cells = dict()
            dict_coord_cells['cell_id'] = adata.obs_names[i]
            dict_coord_cells['x'] = adata.obsm['X_dr'][i,0]
            dict_coord_cells['y'] = adata.obsm['X_dr'][i,1]
            dict_coord_cells['z'] = adata.obsm['X_dr'][i,2]
            list_cells.append(dict_coord_cells)
        with open(os.path.join(file_path,'scatter.json'), 'w') as f:
            json.dump(list_cells, f)    
 
        ## output metadata file of cells
        list_metadata = []
        
        dict_colors = dict()
        for ann in ann_list:
            dict_colors[ann] = get_colors(adata,ann)
        for i in range(adata.shape[0]):
            dict_metadata = dict()
            dict_metadata['cell_id'] = adata.obs_names[i]
            for ann in ann_list:
                dict_metadata[ann] = adata.obs[ann].tolist()[i]
                dict_metadata[ann+'_color'] = dict_colors[ann][i]
            list_metadata.append(dict_metadata)
        with open(os.path.join(file_path,'metadata.json'), 'w') as f:
            json.dump(list_metadata, f)

        ## output gene expression of cells
        if(gene_list is not None):
            print('Generating gene expression of cells ...')
            df_genes = pd.DataFrame(adata.X,index=adata.obs_names,columns=adata.var_names)
            cm = mpl.cm.get_cmap()
            for g in gene_list:
                list_genes = []
                norm = mpl.colors.Normalize(vmin=0, vmax=max(df_genes[g]),clip=True)
                for x in adata.obs_names:
                    dict_genes = dict()
                    dict_genes['cell_id'] = x
                    dict_genes['color'] = mpl.colors.to_hex(cm(norm(df_genes.loc[x,g])))
                    list_genes.append(dict_genes)
                with open(os.path.join(file_path,'gene_'+g+'.json'), 'w') as f:
                    json.dump(list_genes, f) 
        print('STREAM: cells finished!')
        with open(os.path.join(file_path,'index.json'), 'w') as f:
            json.dump({ "tool": 'STREAM' }, f)        
        shutil.make_archive(base_name=os.path.join(adata.uns['workdir'],file_name), format='zip',root_dir=file_path)
        shutil.rmtree(file_path)
    except:
        print('STREAM report failed!')
        raise
    else:
        print(file_name + '.zip is saved at ' + adata.uns['workdir'])
def save_web_report(adata,n_markers=5,file_name='stream_web_report',preference=None,
                    title="experiment name",
                    description="experiment description",
                    starting_node="root node",
                    command_used="command used to generate the analysis",**kwargs):
    """save stream report for the interactive website visualization. A zip file will be generated under the work directory
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    file_name: `str`, optional (default: 'stream_report')
        Ouput Zip file name.
    n_markers: `int`, optional (default: 5)
        Number of top markers selected from each output marker gene file. 
    title: `str`, optional (default: 'experiment name')
        Experiment name displayed on STREAM website.
    description: `str`, optional (default: 'experiment description')
        Experiment description displayed on STREAM website.
    starting_node: `str`, optional (default: 'root node')
        Root node displayed on STREAM website.
    command_used: `str`, optional (default: 'command used to generate the analysis') 
        Command used for the analysis displayed on STREAM website.
    preference: `list`, optional (default: None): 
        The preference of nodes. The branch with speficied nodes are preferred and put on the top part. 
        The higher ranks the node have, the closer to the top the branch with that node is.
    **kwargs: additional arguments to plot_stream()
    
    Returns
    -------
    None
    """    
    if('flat_tree' not in adata.uns_keys()):
        raise ValueError('Please run st.plot_flat_tree(adata) before saving web report')
    if('leaf_markers' not in adata.uns_keys()):
        raise ValueError('Please run st.detect_leaf_markers(adata) before saving web report')
    if('transition_markers' not in adata.uns_keys()):
        raise ValueError('Please run st.detect_transistion_markers(adata) before saving web report')
    if(('de_markers_greater' not in adata.uns_keys()) or ('de_markers_less' not in adata.uns_keys())):
        raise ValueError('Please run st.detect_de_markers(adata) before saving web report')
    
    rootdir = os.path.join(adata.uns['workdir'],file_name)
    reportdir = os.path.join(rootdir,file_name)
    if(not os.path.exists(reportdir)):
        os.makedirs(reportdir)
    epg = adata.uns['epg']
    flat_tree = adata.uns['flat_tree']
    dict_nodes_pos = nx.get_node_attributes(epg,'pos')
    dict_nodes_label = nx.get_node_attributes(flat_tree,'label')
    dict_label_node = {value: key for key,value in nx.get_node_attributes(flat_tree,'label').items()}

    adata.obs[['label','label_color']].to_csv(os.path.join(reportdir,'metadata.tsv'),sep = '\t',index = True)

    print('Generating coordinates of cells in 3D plots...')
    #coordinates of cells in 3D plots
    df_sample = adata.obs[['label','label_color']].copy()
    df_coord = pd.DataFrame(adata.obsm['X_dr'],index=adata.obs_names)
    color = df_sample.sample(frac=1,random_state=100)['label_color'] 
    coord = df_coord.sample(frac=1,random_state=100)
    df_coord_cells = pd.concat([color, coord], axis=1)
    df_coord_cells = df_coord_cells.iloc[:,:4]
    df_coord_cells.columns = ['color','D0','D1','D2']
    df_coord_cells.to_csv(os.path.join(reportdir,'coord_cells.csv'),sep='\t')

    print('Generating coordinates of curves in 3D plots...')
    #coordinates of curves in 3D plots
    for edge_i in flat_tree.edges():
        branch_i_nodes = flat_tree.edges[edge_i]['nodes']
        branch_i_pos = np.array([dict_nodes_pos[i] for i in branch_i_nodes])
        df_coord_curve_i = pd.DataFrame(branch_i_pos)
        df_coord_curve_i.to_csv(os.path.join(reportdir, 'coord_curve_'+dict_nodes_label[edge_i[0]] + '_' + dict_nodes_label[edge_i[1]]+'.csv'),sep='\t',index=False)

    print('Generating coordinates of states/nodes in 3D plots...')
    #coordinates of states(nodes) in 3D plots
    df_coord_states = pd.DataFrame(columns=[0,1,2])
    for x in dict_nodes_label.keys():
        df_coord_states.loc[dict_nodes_label[x]] = dict_nodes_pos[x][:3]
    df_coord_states.sort_index(inplace=True)
    df_coord_states.to_csv(os.path.join(reportdir,'coord_states.csv'),sep='\t')

    print('Generating coordinates of cells in flat tree...')
    #coordinates of cells in flat tree
    df_sample = adata.obs[['label','label_color']].copy()
    df_coord = pd.DataFrame(adata.obsm['X_spring'],index=adata.obs_names)
    color = df_sample.sample(frac=1,random_state=100)['label_color'] 
    coord = df_coord.sample(frac=1,random_state=100)
    df_flat_tree_coord_cells = pd.concat([color, coord], axis=1)
    df_flat_tree_coord_cells.columns = ['color','D0','D1']
    df_flat_tree_coord_cells.to_csv(os.path.join(reportdir, 'flat_tree_coord_cells.csv'),sep='\t')

    print('Generating coordinates of flat tree...')
    #nodes in flat tree
    df_nodes_flat_tree = pd.DataFrame(columns=['D0','D1'])
    for x in dict_nodes_label.keys():
        df_nodes_flat_tree.loc[dict_nodes_label[x]] = flat_tree.node[x]['pos_spring']
    df_nodes_flat_tree.sort_index(inplace=True)
    df_nodes_flat_tree.to_csv(os.path.join(reportdir,'nodes.csv'),sep='\t')

    #edges in flat tree
    df_edges_flat_tree = pd.DataFrame(columns=[0,1])
    for i,edge_i in enumerate(flat_tree.edges()):
        df_edges_flat_tree.loc[i] = [dict_nodes_label[edge_i[0]],dict_nodes_label[edge_i[1]]]
        df_edges_flat_tree.to_csv(os.path.join(reportdir,'edges.tsv'),sep = '\t',index = False,header=False)

    #mapping between original gene names and slugified gene names
    df_genenames = pd.DataFrame(index=adata.var_names,columns=['ori','converted'])
    for x in df_genenames.index:
        df_genenames.loc[x,['ori','converted']] = [x,slugify(x)]
    df_genenames.to_csv(os.path.join(reportdir,'gene_conversion.tsv'),sep = '\t',index = False,header=True)

    print('Generating subway map plots and stream plots of cells...')
    #subway map plots and stream plots of cells
    list_node_start = dict_label_node.keys()
    for ns in list_node_start:
        plot_stream_sc(adata,dist_scale=0.5,root=ns,preference=preference,save_fig=False,fig_path=reportdir)
        plot_stream(adata,root=ns,preference=preference,save_fig=True,fig_path=reportdir,fig_format='png',**kwargs) 

    gene_list = []

    print('Generating transition markers...')
    #transition markers
    file_path = os.path.join(reportdir,'transition_markers')
    if(not os.path.exists(file_path)):
        os.makedirs(file_path)   
    for edge_i in adata.uns['transition_markers'].keys():
        df_tg_i = adata.uns['transition_markers'][edge_i]
        df_tg_i.to_csv(os.path.join(file_path,'transition_markers_'+ edge_i[0]+'_'+edge_i[1] + '.tsv'),sep = '\t',index = True)
        gene_list = gene_list + df_tg_i.index[:n_markers].tolist() 

    print('Generating DE markers...')
    #DE markers
    file_path = os.path.join(reportdir,'de_markers')
    if(not os.path.exists(file_path)):
        os.makedirs(file_path)
    for pair_i in adata.uns['de_markers_greater'].keys():  
        df_de_i_greater = adata.uns['de_markers_greater'][pair_i]
        df_de_i_greater.to_csv(os.path.join(file_path,'de_markers_greater_'+pair_i[0][0]+'_'+pair_i[0][1] + ' and '\
                                            + pair_i[1][0]+'_'+pair_i[1][1] + '.tsv'),sep = '\t',index = True) 
        gene_list = gene_list + df_de_i_greater.index[:n_markers].tolist() 
    for pair_i in adata.uns['de_markers_less'].keys():
        df_de_i_less = adata.uns['de_markers_less'][pair_i]
        df_de_i_less.to_csv(os.path.join(file_path,'de_markers_less_'+pair_i[0][0]+'_'+pair_i[0][1] + ' and '\
                                         + pair_i[1][0]+'_'+pair_i[1][1] + '.tsv'),sep = '\t',index = True)
        gene_list = gene_list + df_de_i_less.index[:n_markers].tolist()

    print('Generating leaf markers...')
    #leaf markers
    file_path = os.path.join(reportdir,'leaf_markers')
    if(not os.path.exists(file_path)):
        os.makedirs(file_path)   
    for leaf_i in adata.uns['leaf_markers'].keys():  
        df_lg_i =  adata.uns['leaf_markers'][leaf_i]
        df_lg_i.to_csv(os.path.join(file_path,'leaf_markers'+leaf_i[0]+'_'+leaf_i[1] + '.tsv'),sep = '\t',index = True)
        gene_list = gene_list + df_lg_i.index[:n_markers].tolist() 

    gene_list = np.unique(gene_list)
    print('Visualizing ' + str(gene_list.shape[0]) + ' markers from detected top marker markers...') 

    print('Generating subway map plots and stream plots of markers...')
    #subway map plots and stream plots of genes
    df_gene_expr = pd.DataFrame(index= adata.obs_names.tolist(),
                                data = adata[:,gene_list].X,
                                columns=gene_list)
    for root in list_node_start:
        df_sample = adata.obs[['label','label_color']].copy()
        df_coord = pd.DataFrame(adata.obsm['X_stream_'+root],index=adata.obs_names)
        color = df_sample.sample(frac=1,random_state=100)['label_color'] 
        coord = df_coord.sample(frac=1,random_state=100)
        df_subwaymap_coord_cells = pd.concat([color, coord], axis=1)
        df_subwaymap_coord_cells.columns = ['color','D0','D1']   
        df_subwaymap_coord_cells.to_csv(os.path.join(reportdir,root,'subway_coord_cells.csv'),sep='\t')
        for edge_i in adata.uns['stream_'+root]['edges'].keys():
            df_edge_pos =  pd.DataFrame(adata.uns['stream_'+root]['edges'][edge_i])
            df_edge_pos.to_csv(os.path.join(reportdir, root, 'subway_coord_line_'+dict_nodes_label[edge_i[0]] + '_' + dict_nodes_label[edge_i[1]]+'.csv'),sep='\t',index=False)
        for g in gene_list:
            gene_expr = df_gene_expr[g].copy()
            df_subwaymap_gene_expr = pd.Series(gene_expr).sample(frac=1,random_state=100)
            df_subwaymap_coord_cells_expr = pd.concat([df_subwaymap_coord_cells[['D0','D1']],df_subwaymap_gene_expr], axis=1)
            df_subwaymap_coord_cells_expr.to_csv(os.path.join(reportdir,root, 'subway_coord_' + slugify(g) + '.csv'),sep='\t')
        plot_stream(adata,color=gene_list,root=root,preference=preference,
                       save_fig=True,fig_path=reportdir,fig_format='png',**kwargs) 
    dict_analysis = {}
    dict_analysis['title']=title
    dict_analysis['description'] = description
    dict_analysis['starting_node'] = starting_node
    dict_analysis['command_used'] = command_used     
    with open(os.path.join(rootdir,'analysis_description.json'), 'w') as fp:
        json.dump(dict_analysis, fp,indent=2)
        
    print('Zipping the files...')
    shutil.make_archive(base_name=os.path.join(adata.uns['workdir'],file_name), format='zip',root_dir=rootdir)
    shutil.rmtree(rootdir)
    print('Done!')
def ltr_dimension_reduction(adata, tsneLearningRate, tsnePerplexity, tsneEarlyExaggeration, n_components=3 ):
    """Perform dimension reduction on LTR data.

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    n_components: `int`, optional (default: 3)
        Number of components to keep.

    Returns
    -------
    updates `adata` with the following fields.

    X_dr : `numpy.ndarray` (`adata.obsm['X_dr']`)
        A #observations × n_components data matrix after dimension reduction.
    X_mlle : `numpy.ndarray` (`adata.obsm['X_mlle']`)
        Store #observations × n_components data matrix after mlle.
    X_se : `numpy.ndarray` (`adata.obsm['X_se']`)
        Store #observations × n_components data matrix after spectral embedding.
    X_umap : `numpy.ndarray` (`adata.obsm['X_umap']`)
        Store #observations × n_components data matrix after umap.
    X_pca : `numpy.ndarray` (`adata.obsm['X_pca']`)
        Store #observations × n_components data matrix after pca.
    trans_mlle : `sklearn.manifold.locally_linear.LocallyLinearEmbedding` (`adata.uns['trans_mlle']`)
        Store mlle object
    trans_se : `sklearn.manifold.spectral_embedding_.SpectralEmbedding` (`adata.uns['trans_se']`)
        Store se object
    trans_umap : `umap.UMAP` (`adata.uns['trans_umap']`)
        Store umap object
    trans_pca : `sklearn.decomposition.PCA` (`adata.uns['trans_pca']`)
        Store pca object
    """
    input_data = adata.X
    # For pcoa
    # from skbio.stats.ordination import pcoa
    # tmpRel = pcoa(input_data)
    # adata.obsm['X_dr'] = (tmpRel.samples.iloc[:, 0:n_components]).to_numpy()

    # For tsne
    tsne = TSNE(n_components=3, metric='precomputed', learning_rate=tsneLearningRate, n_iter=250, perplexity=tsnePerplexity, min_grad_norm=1e-7, n_iter_without_progress=300, early_exaggeration=tsneEarlyExaggeration, random_state=41)
    adata.obsm['X_dr'] = tsne.fit_transform(input_data)

    # For isomap
    # from sklearn.manifold import Isomap
    # isomap = Isomap(n_neighbors=10, n_components=3, metric='precomputed', n_jobs=48)
    # adata.obsm['X_dr'] = isomap.fit_transform(input_data)

    # For tsne on real data
    # from sklearn.manifold import TSNE
    # from skbio.stats.ordination import pcoa
    # tmpRel = pcoa(input_data)
    # pcoaRel = (tmpRel.samples.iloc[:, 0:50]).to_numpy()
    # tsne = TSNE(n_components=3, learning_rate=5, n_iter=250, init='random')
    # adata.obsm['X_dr'] = tsne.fit_transform(pcoaRel)


    if ('params' not in adata.uns_keys()):
        adata.uns['params'] = dict()
    adata.uns['params']['dimension_reduction'] = {'feature': 'all', 'method': 'pcoa', 'n_components': n_components}
    return None
def calOrthogonalDis(centPos, aPos, bPos):
    zoomP = 8
    oaVec = aPos-centPos
    obVec = bPos-centPos
    oaLen = np.linalg.norm(oaVec)
    obLen = np.linalg.norm(obVec)
    if oaLen<obLen:
        aPos, bPos = bPos, aPos
        oaVec, obVec = obVec, oaVec
        oaLen, obLen = obLen, oaLen
    abAngle = calAngle(None, oaVec, obVec)
    if abAngle>np.pi/2:
        return 1e9
    posPairList = [(aPos, bPos), (bPos, aPos)]
    dis = 0
    for posPair in posPairList:
        aPos = posPair[0]
        bPos = posPair[1]
        abVec = bPos-aPos
        angle = calAngle(None, oaVec, abVec)
        if angle>np.pi/2:
            angle = np.pi - angle
        abLen = np.sqrt(abVec.dot(abVec.T))
        tmpDis = abLen*np.cos(angle)/zoomP + abLen*np.sin(angle)*zoomP
        if tmpDis>dis:
            dis=tmpDis
    return dis
def calDisPenalWithMe(paras):
    posData = paras[0]
    ind = paras[1]
    penalFunc = paras[2]
    centerPos = paras[3]
    disList = []
    for i in range(posData.shape[0]):
        if i == ind:
            disList.append(0)
            continue
        angle = calAngle(centerPos, posData[ind], posData[i])
        disDir = posData[ind] - posData[i]
        # l = np.sqrt(disDir.dot(disDir.T))
        l = calOrthogonalDis(centerPos, posData[ind], posData[i])
        rel = l + l * penalFunc(angle)
        disList.append(rel)
    return disList
def z2p3DPos(zPos):
    x, y, z = zPos[0], zPos[1], zPos[2]
    phi = np.angle(x + y*1j)
def calAngle(centPos, aPos, bPos):
    x = None
    y = None
    if centPos is not None:
        x = aPos - centPos
        y = bPos - centPos
    else:
        x = aPos
        y = bPos

    l_x = np.sqrt(x.dot(x.T))
    l_y = np.sqrt(y.dot(y.T))
    dot = x.dot(y.T)
    cosValue = dot / (l_x * l_y)
    if cosValue>1:
        cosValue=1
    if cosValue<-1:
        cosValue=-1
    angle = np.arccos(cosValue)
    return angle
def ppSquarePenal(angle):
    x = angle / (np.pi/2)
    return 2 * x * x
def calPoint2LineDis(linePa, linePb, pPos):
    p1, p2, p3 = linePa, linePb, pPos
    return np.linalg.norm(np.cross(p2-p1, p1-p3))/np.linalg.norm(p2-p1)
def ppNoPenal(angle):
    return 0
def plotCluster3D(adata, labelList, toPlotPosDf=None, outGifFileName = 'cluster3D.gif'):
    from meiPlottter import p3DPlotter
    import pandas as pd
    import numpy as np
    from ttStream import loadWorkData
    import ttUtils
    ltrParaFile = adata.uns["ltrParaFile"]
    figureD = f'{adata.uns["workdir"]}/../../figure'
    typList = []
    for l in labelList:
        if l>=0:
            typList.append(1)
        else:
            typList.append(0)
    toPlotDf = None
    if toPlotPosDf is None:
        toPlotDf = pd.DataFrame(adata.obsm['X_dr'])
        toPlotDf = pd.concat([toPlotDf, pd.DataFrame(dict(label=labelList))], axis=1)
        toPlotDf.columns = ['x', 'y', 'z', 'label']
    else:
        toPlotDf = toPlotPosDf.copy(deep=True)
        toPlotDf = toPlotDf[['x', 'y', 'z']]
        toPlotDf = pd.concat([toPlotDf,
                              pd.DataFrame(dict(label=labelList), index=toPlotDf.index)],
                             axis=1)
    sizeList = [20 for i in range(toPlotDf.shape[0])]
    clusterCenter = toPlotDf[toPlotDf['label']>0].groupby('label').mean()
    for i in range(clusterCenter.shape[0]):
        typList.append(2)
        sizeList.append(100)
    toPlotDf = toPlotDf[['x', 'y', 'z']]
    toPlotDf = pd.concat([toPlotDf, clusterCenter], axis=0)
    typId2annot = {1:'Clustered', 2:'Cluster Center'}
    _3dPlotter = p3DPlotter(initSource=None, initTyp=None, colorList=None, ltrParaFile=ltrParaFile)
    _3dPlotter.plot(adata, typList, typId2annot, f'{figureD}/{outGifFileName}', toPlotDf, False, False, sizeList)
def pfCordZeroIsAnsZero(adata, posDf):
    danteD = f"{adata._uns['workdir']}/.."
    modId2modLen = ttUtils.getModId2modLen(f'{danteD}/toCalRest.modSeq2modId.tab')
    modLenList = []
    posDf = posDf[['x', 'y', 'z']]
    disSer = posDf.apply(np.linalg.norm, axis=1)
    modLenSer = posDf.apply(lambda x: modId2modLen[x.name], axis=1)
    tmpDf = pd.concat([disSer, modLenSer], axis=1)
    tmpDf.columns = ['centOulaDis', 'modLen']
    minModLen = tmpDf['modLen'].min()
    tmpDf = tmpDf.sort_values(by='centOulaDis')
    topMinModLen = tmpDf.iloc[:int(tmpDf.shape[0]/5),:]['modLen'].min()
    posDf = posDf[tmpDf['modLen']>=topMinModLen]
    return posDf, topMinModLen, tmpDf['modLen'].min()
    # return topMinModLen == minModLen, topMinModLen
def getTopShortCenter(adata, posDf, centPerc):
    # posDf should be a dataframe with index that indicate modId
    danteD = f"{adata._uns['workdir']}/.."
    modId2modLen = ttUtils.getModId2modLen(f'{danteD}/toCalRest.modSeq2modId.tab')
    modLenList = []
    totModNum = adata.obsm['X_dr'].shape[0]
    for modId in range(totModNum):
        modLenList.append(modId2modLen[modId])
    tmpDf = pd.DataFrame(dict(modLen=modLenList))
    tmpDf = tmpDf.loc[posDf.index]
    if centPerc is None:
        tmpDf = tmpDf[tmpDf['modLen'] == tmpDf['modLen'].min()]
    else:
        tmpDf = tmpDf.sort_values(by='modLen')
        tmpDf = tmpDf.iloc[:max(2,int(tmpDf.shape[0] * centPerc)), ]
    centerPos = np.array([posDf[['x', 'y', 'z']].loc[list(tmpDf.index)].mean(axis=0)])
    return centerPos[0]
def caliPosDfWithNewCenter(adata, posDf, centPerc):
    # posDf sould be with and only with x, y, z three columns
    centerPos = getTopShortCenter(adata, posDf, centPerc)
    posDf = posDf - centerPos
    posDf = posDf[~((posDf['x']==0)&(posDf['y']==0)&(posDf['z']==0))]
    return posDf
def partitionIntoBlocks(posDf, adata):
    sc = AgglomerativeClustering(n_clusters=6, affinity='cosine', linkage='average').fit(posDf)
    return sc.labels_
def testToCluNum(cosineDisMat, prefix=None, figureD=None):
    silhouetteScoreList = []
    for i in range(2,10):
        sc = AgglomerativeClustering(n_clusters=i, affinity='precomputed', linkage='average').fit(cosineDisMat)
        silhouetteScoreList.append(silhouette_score(cosineDisMat, sc.labels_, metric='precomputed'))
    return silhouetteScoreList.index(max(silhouetteScoreList))+2
def getValidBranchCluIdSet(posDf, minNumOfBranch, adata, numPerc=0.9, piPerc=0.1):
    def removeCentRelaPoint(posDf):
        tmpPosDf = posDf.copy(deep=True)
        tmpPosDf['centOulaDis'] = posDf.apply(np.linalg.norm, axis=1)
        posDf = posDf[ tmpPosDf['centOulaDis'] >= (tmpPosDf['centOulaDis'].max()/10) ]
        return posDf
    def getCluIdPiPerc(posDf, numPerc, minNumOfBranch, adata):
        if posDf.shape[0] < minNumOfBranch:
            return 1
        posDf = posDf[['x', 'y', 'z']]
        # posDf = removeCentRelaPoint(posDf)
        centPos = posDf.mean()
        posDf = posDf[~((posDf['x']==centPos['x'])&(posDf['y']==centPos['y'])&(posDf['z']==centPos['z']))]
        # cosineDisAr = cosine_distances(np.array([centPos]), posDf).astype('float64')
        cosineDisAr = cosine_distances(np.array([centPos]), posDf)
        return np.percentile(cosineDisAr, numPerc*100)
    piPercSer = posDf.groupby('cluId').apply(getCluIdPiPerc, numPerc=numPerc, minNumOfBranch=minNumOfBranch, adata=adata)
    numSer = posDf['cluId'].value_counts()
    toShowDf = pd.concat([piPercSer, numSer],  axis=1)
    print(toShowDf.head(10))
    branchCluIdSet = set()
    for ind, pi in piPercSer.items():
        if pi < piPerc:
            branchCluIdSet.add(ind)
    return branchCluIdSet, piPercSer
def preZoomInValidBranchCluId(adata, prefix, oriDisMat, tarIndList, tsneLearningRate, tsnePerplexity, tsneEarlyExaggeration, minNumOfBranch, cluCentCut, zoomInLevel, maxZoomInLevel):
    print('===')
    print('Pre:', prefix)
    tarDisMat = oriDisMat[tarIndList, :][:, tarIndList]

    tsne = TSNE(n_components=3, metric='precomputed', learning_rate=tsneLearningRate, n_iter=250,
                perplexity=tsnePerplexity, min_grad_norm=1e-7, n_iter_without_progress=300,
                early_exaggeration=tsneEarlyExaggeration, random_state=41)
    posDf = pd.DataFrame(tsne.fit_transform(tarDisMat))
    posDf.index = tarIndList
    posDf.columns = ['x', 'y', 'z']

    tsne = TSNE(n_components=2, metric='precomputed', learning_rate=tsneLearningRate, n_iter=250,
                perplexity=tsnePerplexity, min_grad_norm=1e-7, n_iter_without_progress=300,
                early_exaggeration=tsneEarlyExaggeration, random_state=41)
    posDf_2d = pd.DataFrame(tsne.fit_transform(tarDisMat))
    posDf_2d.index = tarIndList
    posDf_2d.columns = ['x', 'y']

    # if not pfCordZeroIsAnsZero(adata, posDf):
    #     posDf = caliPosDfWithNewCenter(adata, posDf, 0.01)
    posDf, topMinModLen, globMinModLen = pfCordZeroIsAnsZero(adata, posDf)
    posDf_2d = posDf_2d.loc[posDf.index]
    if (maxZoomInLevel>0 and zoomInLevel>maxZoomInLevel) \
            or topMinModLen > globMinModLen+1:
        print('===')
        return [], set(), posDf, posDf_2d, pd.Series([])

    # posDf = caliPosDfWithNewCenter(adata, posDf, 0.01)
    cosineDisMat = cosine_distances(posDf).astype('float64')
    toCluNum = testToCluNum(cosineDisMat)
    cluIdList = AgglomerativeClustering(n_clusters=toCluNum, affinity='precomputed', linkage='average').fit(cosineDisMat).labels_
    posDf = pd.concat([posDf,
                       pd.DataFrame(dict(cluId=cluIdList), index=posDf.index)],
                      axis=1)
    validBranchCluIdSet, centPercScore = getValidBranchCluIdSet(posDf, minNumOfBranch, adata, piPerc=cluCentCut)
    print('===')
    return cluIdList, validBranchCluIdSet, posDf, posDf_2d, centPercScore
def plotAllZoomIn(ii, zoomInAdataDict, prefix, sim=False):
    def plotEach(infoDf, sim=False):
        infoDf = infoDf.drop_duplicates(['modId'])
        infoDf.index = list(infoDf.modId)
        if sim:
            infoDf['Class'] = infoDf.oriId.apply(lambda x:'Simulate_' + chr(int(x.split('_')[0])+ord('A')))
        posDf = infoDf[['x', 'y', 'z', 'Class']]
        posDf_2d = infoDf[['xx', 'yy', 'Class']]
        posDf.columns = ['x', 'y', 'z', 'cluId']
        posDf_2d.columns = ['x', 'y', 'cluId']
        prefix = infoDf.zoomIn.iloc[0]
        plotZoomIn(zoomInAdataDict[infoDf.zoomIn.iloc[0]], posDf, posDf_2d, prefix, treeAnnot=False)
    def transLabel(branchClass):
        alphaList = branchClass.split('_')[1:]
        numList = [str(ord(al)-ord('A')+1) for al in alphaList]
        return f"zoomIn {'.'.join(numList)}"
    def alterClassLabel(infoDf, zoomInAdataDict, prefix):
        leafClassList = infoDf.Class[infoDf.isTerminalCluster==True].drop_duplicates().sort_values().tolist()
        leafClass2label = {lc:chr(ord('A')+ind) for ind, lc in enumerate(leafClassList)}
        branchClassList = infoDf.Class[infoDf.isTerminalCluster==False].drop_duplicates().sort_values().tolist()
        branchClass2label = {bc:transLabel(bc) for bc in branchClassList}
        Class2label = {**leafClass2label, **branchClass2label}
        Class2label[prefix] = 'All'
        infoDf.Class = infoDf.Class.map(Class2label)
        infoDf.zoomIn = infoDf.zoomIn.map(Class2label)
        infoDf.zoomIn = infoDf.zoomIn.fillna('All')
        print(list(zoomInAdataDict))
        print(list(Class2label))
        newDict = {Class2label[zoomIn]:adata for zoomIn, adata in zoomInAdataDict.items()}
        return infoDf, newDict

    infoDf = ii.infoDf.copy(deep=True)
    infoDf, zoomInAdataDict = alterClassLabel(infoDf, zoomInAdataDict, prefix)
    infoDf.groupby('zoomIn', as_index=False).apply(plotEach, sim=sim)
def plotZoomIn(adata, posDf, posDf_2d, prefix, treeAnnot):
    # posDf should contain column 'cluId' that was used as pltId.
    # posDf_2d should only contain columns 'x' and 'y' and corresponded index.
    # posDf_2d's cluId is updated by that of posDf
    from meiPlottter import p3DPlotter
    # from ltrShower import ltrShower
    import pandas as pd

    # Basic parameters are defined by adata
    ltrParaFile = adata.uns["ltrParaFile"]
    figureD = f'{adata.uns["workdir"]}/../../figure'

    if posDf is None:
        posDf = pd.DataFrame(adata.obsm['X_dr'], index=list(adata.obs['modId']))
        posDf = pd.concat([posDf,
                           pd.DataFrame(dict(cluId=list(adata.obs['cluId'])), index=list(adata.obs['modId']))],
                          axis=1)
        posDf = pd.concat([posDf, adata.uns['oToPltPointPosDf']], axis=0)
    # 2dPlot
    if posDf_2d is None:
        posDf_2d = adata.uns['posDf_2d'][['x', 'y']]
    if 'cluId' not in posDf_2d.columns:
        posDf_2d = pd.concat([posDf_2d,
                                   posDf['cluId']],
                                 axis=1)
    '''
    ls = ltrShower(ltrParaFile)
    ls.plot2DWithCluId(posDf=posDf_2d, outPrefix=prefix)
    '''
    outGifFileName = f'{prefix}.gif'
    _3dPlotter = p3DPlotter(initSource=None, initTyp=None, colorList=None, ltrParaFile=ltrParaFile)
    toPlotDf = posDf[['x', 'y', 'z']]

    dedupCluList = posDf.cluId.drop_duplicates().sort_values().reset_index(drop=True).tolist()
    pltId2annot = {ind+1:clu for ind, clu in enumerate(dedupCluList)}
    annot2pltId = {b:a for a,b in pltId2annot.items()}
    typList = [annot2pltId[annot] for ind, annot in posDf.cluId.items()]
    _3dPlotter.plot(adata, typList, pltId2annot, f'{figureD}/{outGifFileName}', toPlotDf, False, treeAnnot, None)
def transBranch22D():
    pass
def plotZoomInBranchTime(posDf, prefix, modId2insTimeSet, adata, cluId2centScoreSer):
    from meiPlottter import zoomInBranch2dPlot
    centDisList = []
    insTimeList = []
    for ind, row in posDf.iterrows():
        insTimeSet = modId2insTimeSet[ind]
        for insTime in insTimeSet:
            insTimeList.append(insTime)
            centDisList.append(np.linalg.norm(np.array(row)))
    toPlotDf = pd.DataFrame(dict(centDis=centDisList, insTime=insTimeList))
    figureD = f'{adata.uns["workdir"]}/../../figure'
    outPdf = f'{figureD}/{prefix}.insTime.pdf'
    zoomInBranch2dPlot(toPlotDf, outPdf)

    # pcaRel = sklearnPCA(n_components=2).fit_transform(posDf)
    # insTimeList =
    # toPlotDf = pd.concat([pd.DataFrame])
    #
    #
    #
    # print(type(pcaRel))
    # print(pcaRel.shape)
def extractBranchFasta(adata, modIdList, prefix, ltrFastaFile):
    import sys
    danteD = f'{adata.uns["workdir"]}/../'
    fastaD = f'{danteD}../fasta'
    if not os.path.exists(fastaD):
        os.system(f'mkdir -p {fastaD}')
    modId2oriIdSet = ttUtils.getModId2OriIdSet(f'{danteD}/finalModSeq.tab', f'{danteD}/toCalRest.modSeq2modId.tab')
    totOriIdSet = set()
    for modId in modIdList:
        for oriId in modId2oriIdSet[modId]:
            totOriIdSet.add(oriId)
    with pysam.FastaFile(ltrFastaFile) as ltrFa, open(f'{fastaD}/{prefix}.fa', 'w') as outFa:
        for oriId in totOriIdSet:
            print(f'>{oriId}', file=outFa)
            print(ltrFa.fetch(oriId), file=outFa)
def selectIndList(posDf, key, idSet):
    indList, idList = [], []
    for id, (ind, row) in enumerate(posDf.iterrows()):
        if row[key] in idSet:
            idList.append(id)
            indList.append(ind)
    return indList, idList
def transCluId2pltId(posDf, cluId2pltId):
    pltIdSer = posDf['cluId'].apply(lambda x: cluId2pltId[x] if x in cluId2pltId else 0)
    posDf['cluId'] = pltIdSer
    return posDf
def ltr_adjCluIdWithBranchId(adata):
    def getRootNode(tr):
        posSer = pd.Series([tr.nodes[n]['pos'] for n in tr.nodes()], index=[n for n in tr.nodes()])
        return posSer.apply(np.linalg.norm).idxmin()
    def dfs_getNodeId2rootId(nowNode, tr, rootNode, visedNode, branchId2rootId):
        visedNode.add(nowNode)
        for n in tr.neighbors(nowNode):
            branchId2rootId[tr[nowNode][n]['id']] = rootNode
        for n in tr.neighbors(nowNode):
            if n not in visedNode:
                dfs_getNodeId2rootId(n, tr, rootNode, visedNode, branchId2rootId)
    def getBranchId2rootId(rootNode, tr):
        branchId2rootId = {}
        for n in tr.neighbors(rootNode):
            visedNode = set([rootNode])
            dfs_getNodeId2rootId(n, tr, n, visedNode, branchId2rootId)
        return branchId2rootId

    rootNode = getRootNode(adata.uns['flat_tree'])
    branchId2rootId = getBranchId2rootId(rootNode, adata.uns['flat_tree'])
    newCluIdList = [branchId2rootId[branchId] for branchId in list(adata.obs['branch_id'])]
    cluIdSet = set(list(dict.fromkeys(newCluIdList)))
    dedupCluIdList = sorted(list(cluIdSet))
    # 0 is for noise, which would be merged later (outside of this function)
    newCluIdList = [dedupCluIdList.index(i)+1 for i in newCluIdList]
    rootId2cluId = {rootId:(i+1) for i,rootId in enumerate(dedupCluIdList)}
    adata.uns['branchId2cluId'] = {branchId:rootId2cluId[rootId] for branchId,rootId in branchId2rootId.items()}
    adata.obs['cluId'] = newCluIdList
def reclusterWithTraj(prefix, adata, posDf, posDf_2d, epgLambda, epgAlpha, epgMu):
    # Only noises were prevent from trajectory reconstruction.
    modId2cluId = {modId:cluId for modId, cluId in zip(list(posDf.index), list(posDf.cluId))}
    validCluIdSet = set(list(np.array(range(posDf.cluId.max()))+1))
    toTrajIndList, toTrajZeroIndList = selectIndList(posDf, 'cluId', validCluIdSet)
    subAdata = ad.AnnData(X=adata.X[toTrajIndList, :][:, toTrajIndList])
    subAdata.obs['modId'] = toTrajIndList
    subAdata.obs['cluId'] = list(posDf.loc[toTrajIndList]['cluId'])
    subAdata.obsm['X_dr'] = np.array(posDf.loc[toTrajIndList][['x', 'y', 'z']])
    subAdata.uns['workdir'] = adata.uns['workdir']
    subAdata.uns['params'] = adata.uns['params']
    subAdata.uns['prefix'] = prefix
    subAdata.uns['oToPltPointPosDf'] = posDf[posDf.cluId==0]
    subAdata.uns['posDf_2d'] = posDf_2d
    branchNodeNum = 4
    ltr_seed_elastic_principal_graph(subAdata)
    elastic_principal_graph(subAdata, epg_lambda=epgLambda, epg_alpha=epgAlpha, epg_mu=epgMu,
                            epg_n_nodes=len(subAdata.obs['cluId'].drop_duplicates()) * branchNodeNum + 1)
    ltr_adjCluIdWithBranchId(subAdata)
    for modId,cluId in zip(list(subAdata.obs['modId']), list(subAdata.obs['cluId'])):
        modId2cluId[modId] = cluId
    newCluIdSer = pd.Series([modId2cluId[modId] for modId in list(posDf.index)], index=posDf.index)
    posDf['cluId'] = newCluIdSer
    return subAdata, posDf
def reclusterNoise(prefix, nowAdata, posDf, posDf_2d):
    def getBranchId2cutDist(adata, q):
        tmpDf = pd.DataFrame(dict(branchId=list(adata.obs['branch_id']),
                                  dist=list(adata.obs['branch_dist'])))
        return tmpDf.groupby('branchId')['dist'].apply(lambda x:x.quantile(q))
    def getClosestBranchIdAndDist(ser, adata, list_br_id, flat_tree, epg, dict_nodes_pos, branchId2cutDist):
        xp = np.array(ser[['x', 'y', 'z']])
        dict_br_matrix = dict()
        for br_id in list_br_id:
            dict_br_matrix[br_id] = np.array([dict_nodes_pos[i] for i in flat_tree.edges[br_id]['nodes']])
        dict_results = dict()
        list_dist_xp = list()
        for br_id in list_br_id:
            dict_results[br_id] = project_point_to_line_segment_matrix(dict_br_matrix[br_id], xp)
            list_dist_xp.append(dict_results[br_id][2])
        br_id = list_br_id[np.argmin(list_dist_xp)]
        results = dict_results[br_id]
        dist = results[2]
        if dist<=branchId2cutDist[br_id]:
            return adata.uns['branchId2cluId'][br_id]
        else:
            return 0

    branchId2cutDist = getBranchId2cutDist(nowAdata, q=.8)

    list_br_id = [br_id for br_id in list(pd.Series(nowAdata.obs['branch_id']).drop_duplicates())]
    flat_tree = nowAdata.uns['flat_tree']
    epg = nowAdata.uns['epg']
    dict_nodes_pos = nx.get_node_attributes(epg, 'pos')
    noiseCluIdSer = posDf[posDf['cluId']==0].apply(getClosestBranchIdAndDist, adata=nowAdata, list_br_id=list_br_id,
                                                   flat_tree=flat_tree, epg=epg, dict_nodes_pos=dict_nodes_pos,
                                                   branchId2cutDist=branchId2cutDist, axis=1)
    posDf.loc[noiseCluIdSer.index, 'cluId'] = noiseCluIdSer
    return posDf, posDf_2d
def removeAdjZeroInd(posDf, indList):
    cordDisSer = posDf.loc[indList, ['x', 'y', 'z']].apply(np.linalg.norm, axis=1)
    cordDisCut = 0.1 * np.percentile(cordDisSer, 99.5)
    return list(cordDisSer[cordDisSer > cordDisCut].index)
def zoomIn(adata, oriDisMat, tarIndList, tsneLearningRate, tsnePerplexity, tsneEarlyExaggeration, epgLambda, epgAlpha, epgMu, prefix, oriId2insTimeSet, ltrFastaFile, zoomInAdataDict, zoomInLevel, maxZoomInLevel, minNumOfBranch=50, cluIdList=None, validCluIdSet=None, posDf=None, posDf_2d=None, max_tsneEarlyExaggeration=None, cluCentCut=0.1):
    def calValidP(validCluIdSet, posDf):
        if len(validCluIdSet)==0:
            return 0
        return np.sum(posDf.cluId.isin(validCluIdSet))/posDf.shape[0]
    if max_tsneEarlyExaggeration is None:
        max_tsneEarlyExaggeration = 1.5*tsneEarlyExaggeration
    validPCutOff = 0.9
    if cluIdList is None:
        validP = 0
        while validP<validPCutOff and tsneEarlyExaggeration<=max_tsneEarlyExaggeration:
            print(tsneEarlyExaggeration, max_tsneEarlyExaggeration)
            cluIdList, validCluIdSet, posDf, posDf_2d, cluId2centScoreSer = preZoomInValidBranchCluId(adata, prefix, oriDisMat, tarIndList, tsneLearningRate, tsnePerplexity, tsneEarlyExaggeration, minNumOfBranch, cluCentCut=cluCentCut, zoomInLevel=zoomInLevel, maxZoomInLevel=maxZoomInLevel)
            validP = calValidP(validCluIdSet, posDf)
            tsneEarlyExaggeration += 0.5
        if validP < validPCutOff:
            return infoIntegrator.init_Null(adata.uns["ltrParaFile"])
    upperCase = 'A'
    cluId2pltId = {}
    notZoomInPltIdSet = set()
    inte = infoIntegrator.init_Null(adata.uns["ltrParaFile"])
    for cluId in sorted(list(validCluIdSet)):
        upperCase = chr(ord(upperCase) + 1)
        cluId2pltId[cluId] = ord(upperCase) - ord('A')
    totPltIdSet = set([cluId2pltId[i] for i in cluId2pltId])
    totPltIdSet.add(0)
    # The zeros in posDf.cluId should not be modified
    posDf = transCluId2pltId(posDf, cluId2pltId)
    nowAdata, posDf = reclusterWithTraj(prefix, adata, posDf, posDf_2d, epgLambda, epgAlpha, epgMu)
    if len(list(posDf.cluId.drop_duplicates()))<2:
        return inte
    posDf, posDf_2d = reclusterNoise(prefix, nowAdata, posDf, posDf_2d)
    zoomInAdataDict[prefix] = nowAdata
    # plotZoomIn(nowAdata, posDf, posDf_2d, prefix, treeAnnot=False)
    cluId2isTer = {}
    for cluId in range(1, max(list(posDf['cluId']))+1):
        isTer = True
        tmpIndList = list(posDf[posDf['cluId'] == cluId].index)
        tmpTsneEarlyExaggeration = max(2, tsneEarlyExaggeration / 3)
        tmpTsnePerplexity = tsnePerplexity
        tmpPrefix = f"{prefix}_{chr(ord('A')+cluId-1)}"
        validP = 0
        preZoomInCluIdList, preZoomInValidCluIdSet, tmpPosDf, tmpPosDf_2d, cluId2centScoreSer = None, None, None, None, None
        while validP < validPCutOff and tmpTsneEarlyExaggeration <= tsneEarlyExaggeration:
            preZoomInCluIdList, preZoomInValidCluIdSet, tmpPosDf, tmpPosDf_2d, cluId2centScoreSer = preZoomInValidBranchCluId(
                adata, tmpPrefix, oriDisMat, tmpIndList, tsneLearningRate, tmpTsnePerplexity, tmpTsneEarlyExaggeration,
                minNumOfBranch, cluCentCut=cluCentCut, zoomInLevel=zoomInLevel+1, maxZoomInLevel=maxZoomInLevel)
            validP = calValidP(preZoomInValidCluIdSet, tmpPosDf)
            tmpTsneEarlyExaggeration += 0.5
        if len(preZoomInValidCluIdSet)>1 and validP>=validPCutOff:
            subInte = zoomIn(adata, oriDisMat, tmpIndList, tsneLearningRate, tsnePerplexity, tmpTsneEarlyExaggeration,
                             epgLambda, epgAlpha, epgMu, tmpPrefix, oriId2insTimeSet, ltrFastaFile, zoomInAdataDict,
                             zoomInLevel+1, maxZoomInLevel,
                             minNumOfBranch, preZoomInCluIdList, preZoomInValidCluIdSet, tmpPosDf, tmpPosDf_2d)
            inte = infoIntegrator.merge(inte, subInte)
            isTer=False
        else:
            notZoomInPltIdSet.add(cluId)
        cluId2isTer[cluId] = isTer
    myInte = infoIntegrator.init_from_posDf(posDf, posDf_2d, prefix, adata.uns["ltrParaFile"], cluId2isTer)
    return infoIntegrator.merge(inte, myInte)
def isComplexScenario(posDf, adata):
    # posDf can be a np.array, with each position of point each row
    # index of posDf is modId, can not be changed
    def hdbscanPf(posDf, adata, minValidNum):
        cosDisMat = cosine_distances(posDf).astype('float64')
        cluRel = hdbscan.HDBSCAN(min_cluster_size=minValidNum,
                                 metric='precomputed').fit(cosDisMat)
        cluRel = cluRel.labels_
        return cluRel

    posDf = caliPosDfWithNewCenter(adata, posDf, 0.01)
    pfMethod = hdbscanPf
    cluRel = pfMethod(posDf, adata, max(50, int(posDf.shape[0]/30)))
    cluNum = np.max(cluRel)+1
    return cluNum>4
def ltr_seed_elastic_principal_graph(adata, threads=48):
    def idxquantile(s, q, *args, **kwargs):
        qv = s.quantile(q, *args, **kwargs)
        return (s.sort_values(ascending=False)<=qv).idxmax()
    def calCluIdSegPointPos(posDf, perc=0.95):
        cordDisSer = posDf[['x', 'y', 'z']].apply(np.linalg.norm, axis=1)
        return idxquantile(cordDisSer, 1-perc), idxquantile(cordDisSer,perc)
    def constructInitTree(terPosSer, posDf, branchNodeNum=4):
        ttg = nx.Graph()
        ttg.add_node(0, pos=np.array([0, 0, 0]))
        nodeId = 1
        for ind, terPos in terPosSer.items():
            posA = np.array(posDf[['x', 'y', 'z']].loc[terPos[0]])
            posB = np.array(posDf[['x', 'y', 'z']].loc[terPos[1]])
            ttg.add_node(nodeId, pos=posA)
            ttg.add_edge(0, nodeId)
            nodeId += 1
            for i in range(branchNodeNum-1):
                tmpPos = posA + (posA+posB)/branchNodeNum*(i+1)
                ttg.add_node(nodeId, pos=tmpPos)
                ttg.add_edge(nodeId-1, nodeId)
                nodeId += 1
        return ttg
    def pfValidInitTreeCluId(posDf, totTerminalDis):
        myTerminalDis = posDf[['x', 'y', 'z']].apply(np.linalg.norm, axis=1).quantile(q=0.99)
        return myTerminalDis>totTerminalDis/3
    posDf = pd.DataFrame(adata.obsm['X_dr'])
    posDf.columns = ['x', 'y', 'z']
    posDf.index = list(adata.obs['modId'])
    posDf = pd.concat([posDf,
                       pd.DataFrame(dict(cluId=list(adata.obs['cluId'])), index=posDf.index)],
                      axis=1)

    # To remove noised valid cluId (those total class very closed to center)
    totTerminalDis = posDf[['x', 'y', 'z']].apply(np.linalg.norm, axis=1).quantile(q=0.99)
    validInitTreeCluId = posDf.groupby('cluId').apply(pfValidInitTreeCluId, totTerminalDis=totTerminalDis)
    validInitTreeCluId = pd.Series(list(validInitTreeCluId[validInitTreeCluId==True].index))
    print(validInitTreeCluId)

    terPosSer = posDf[posDf.cluId.isin(validInitTreeCluId)].groupby('cluId').apply(calCluIdSegPointPos)
    initTree = constructInitTree(terPosSer, posDf)
    adata.obs['sc'] = ['cluster ' + str(x) for x in list(adata.obs['cluId'])]
    epg = initTree
    dict_nodes_pos = {n:epg.nodes[n]['pos'] for n in epg.nodes()}
    dict_branches = extract_branches(epg)
    flat_tree = construct_flat_tree(dict_branches)
    nx.set_node_attributes(flat_tree, values={x: dict_nodes_pos[x] for x in flat_tree.nodes()}, name='pos')
    adata.uns['epg'] = deepcopy(epg)
    adata.uns['flat_tree'] = deepcopy(flat_tree)
    adata.uns['seed_epg'] = deepcopy(epg)
    adata.uns['seed_flat_tree'] = deepcopy(flat_tree)
    project_cells_to_epg(adata)
    calculate_pseudotime(adata)

    if ('params' not in adata.uns_keys()):
        adata.uns['params'] = dict()
    adata.uns['params']['epg'] = {'use_vis': False}
    print('Number of initial branches: ' + str(len(dict_branches)))
