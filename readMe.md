# LTR_Stream
LTR_Stream is designed to achieve sub-lineage level LTR-RT clustering in closely related species, 
discovering valuable genetic markers for genome comparison and LTR-RT modular evolution in host genome.
It takes nucleotide sequences of intact LTR-RTs belonging to the same LTR-lineage as input. A mix 
of LTR-RTs from different LTR-lineages is theoretically acceptable but not recommended. LTR_Stream gives
each LTR-RT a cluster label and automatically evaluates reliability of each cluster.

## Graphical Abstract

<div align=center>
    <img src="https://github.com/xjtu-omics/LTR_Stream/blob/main/.readMe_images/GraphAbstract.png" width="700px" height="443px" />
</div>

## Sub-lineage Clustering of Retand LTR-RTs from Three Papaver Species
<div align=center>
    <img src="https://github.com/xjtu-omics/LTR_Stream/blob/main/.readMe_images/Papaver_Retand_all_cut.gif" width="700px" height="377px" /> 
</div>

## Installation
### Requirements
#### 1.Conda 
Conda should be installed with version >=23.1.0.  
Mamba is recommended for speeding up conda.
#### 2.Git
Please install git with version >=2.34.1.  
Please configure the ssh key of git and make sure `git clone` could work.  
### Install LTR_Stream
#### 1.Clone LTR_Stream from github.
```shell
ltrStreamInstallPath=path_you_want_to_install_LTR_Stream
cd ${ltrStreamInstallPath} && git clone git@github.com:xjtu-omics/LTR_Stream.git
```
#### 2.Run install script.
If mamba is not available, please run:
```shell
cd ${ltrStreamInstallPath}/LTR_Stream && bash Init_LTR_Stream_Env.sh
```
For a speeding up installation with mamba, please run:
```shell
cd ${ltrStreamInstallPath}/LTR_Stream && bash Init_LTR_Stream_Env.sh mamba
```

## Quick start
### Sub-lineage level LTR-RT clustering
```shell
conda activate ltrStream
cd ${ltrStreamInstallPath}/LTR_Stream/src
snakemake -s LTR_Stream.smk -f stream --config ltrParaFile=path_of_ltrPara.tsv -j {threadsNumber}
```

### Config file `ltrPara.tsv`
LTR_Stream will automatically run according to parameters set in this TSV 
(Tab-Separated Values) file, so please make sure all the parameters were set here before 
you start LTR_Stream.smk. (You can modify the file name and path according to your preferences. 
In this documentation, we refer to this configuration parameters file as ltrPara.tsv.)
The following is an example of the file. Lines beginning with a # represent comments. Values of optional parameter in this example represent their default values in LTR_Stream. To facilitate 
parameter debugging, the parameters that significantly impact the clustering results will be introduced first. A standard example of this file is under examples/.

```tsv
# An example for ltrPara.tsv
# All tab seperated.

# Mandatory parameters
# workDir: A blank directory for running LTR_Stream
# The outputs of LTR_Stream are in workDir/figure
workDir /xx/xx/xx

# ltrFasta: The nucleotide sequences of the LTR-RT set you want to 
# analyze. Please ensure it is in standard FASTA format. Names of these
# sequences should follow the format like 'chrxx:stPos-endPos(strand)'. It 
# is recommended to use bedtools to extract sequences from the genome.
ltrFasta    /xx/xx/xx.fa


# Optional parameters

# Important parameters
# minOverLapForNovelModule: Control the number and dispersion of module sequences in the 3-D space.
# It is used in disjoint-set data structure to judge if there should be an edge between two alignment
# regions. It could be set at the range from 0 to 1. Greater minOverLapForNovelModule leads to more
# module sequences and more dispersed result. Default is 0.8.
minOverLapForNovelModule 0.8


# topModNum: Control the number and dispersion of module sequences with minOverLapForNovelModule.
# Greater topModNum leads to more module sequences and more dispersed result. LTR_Stream will output 
# a module number versus covered LTR-RTs (named coverLine.pdf under workDir/figure). The topModNum 
# needs to be set large enough to ensure that about 80% of LTR-RTs have 2-3 modules. It is estimated 
# topModNum should be at range 200-800. Larger minOverLapForNovelModule usually corresponds to larger 
# topModNum. You can adjust the two parameters in coordination. Default is 250.
topModNum   250


# tsneEarlyExaggeration: A crucial parameter in t-SNE dimensionality reduction, directly affects the
# results. An excessively large tsneEarlyExaggeration will result in a linear shape in the 
# three-dimensional space, while an excessively small tsneEarlyExaggeration will lead to a dispersed 
# distribution, hindering sub-lineage identification. It is estimated that tsneEarlyExaggeration 
# should be at range 6-9. Default is 6.
tsneEarlyExaggeration   6


# tsnePerplexity: Larger tsnePerplexity will provide more robust results, while a smaller 
# tsnePerplexity will yield more detailed clustering results. Depending on the size of the dataset, 
# it is not recommended to set tsnePerplexity to less than 3% of the module sequence count for larger 
# datasets, or less than 15 for smaller datasets. Default is 100.
tsnePerplexity  100


# cluCentCut: A parameter used to assess the degree of intra-class distribution aggregation in 3D
# space. A larger cluCentCut will result in coarser clustering. If LTR_Stream indicates clustering 
# failure, please increase this parameter within the range of 0-1. Default is 0.1.
cluCentCut  0.1


# maxZoomInLevel: LTR_Stream achieves fine clustering of LTR-RT in complex scenarios through 
# iterative expansion. This parameter controls the maximum depth of iterative expansion. If you find
# that the number of clusters is too large or some categories within subviews are verified as 
# unreliable, you can set a maximum limit. The default value is -1, which means no limit is set. 
maxZoomInLevel  -1


# Other parameters

# tsneLearningRate: For t-SNE dimensionality reduction, LTR_Stream requires a very small learning rate, with a default value of 6. It is not recommended to set this value higher than 8.
tsneLearningRate    6


# blastEvalue: Used for homology searching in BLASTn. Default is 1e-10. If the LTR-RT sequence set to 
# be analyzed has particularly high similarity, you can reduce this parameter accordingly.
blastEvalue 1e-10


# Parameters used in ElPiGraph
epgLambda   0.01
epgMu   0.01
epgAlpha 0.05
```

### Outputs
All outputs will be saved in workDir/figure
#### workDir/figure/\*\*.gif
GIF files showing clustering results in each 3D-subview.

#### workDir/figure/clusterRel.tsv
TSV file recording final cluster results.

#### workDir/figure/classInfo.tsv
TSV file recording details of clustering including coordinate information in each subview.

#### workDir/figure/clusterNuclVali.tsv
TSV file recording foldchange of inter- and intra-distance and corresponding significance for each cluster. Foldchange that signifcantly greater than one means reliable cluster.

#### workDir/figure/coverLine.pdf
Line plot showing module number and corresponding covered LTR-RT percentage. Used for guiding parameter ajustment.

