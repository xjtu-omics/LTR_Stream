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
cd ${ltrStreamInstallPath}/LTR_Stream && bash install_LTR_Stream.sh
```
For a speeding up installation with mamba, please run:
```shell
cd ${ltrStreamInstallPath}/LTR_Stream && bash install_LTR_Stream.sh mamba
```

## Quick start
#### Sub-lineage level LTR-RT clustering
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
The following is an example of the file. Lines beginning with a hash symbol (#) represent 
comments. Values of optional parameter in this example represent their default values in LTR_Stream. To facilitate 
parameter debugging, the parameters that significantly impact the clustering results will be introduced first.

```tsv
# An example for ltrPara.tsv
# All tab seperated.

# Mandatory parameters
# workDir: A blank directory for running LTR_Stream
# The outputs of LTR_Stream are in workDir/figure
workDir /xx/xx/xx

# ltrFasta: The nucleotide sequences of the LTR-RT set you want to analyze. Please ensure it is in standard FASTA format. Names of these sequences should follow the format like 'chrxx:stPos-endPos(strand)'. It is recommended to use bedtools to extract sequences from the genome.
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






# nClustersForInitTree : Control the initial tree for trajectory reconstruction.
# Greater nClustersForInitTree would generate a initial tree with more nodes, 
# leading to a tree with more details but also more noises. Default is 50.
nClustersForInitTree 50


# minNumInClusterForInitTree: When initializing a tree before performing ElPiGraph,
# clusters with less CEPs are usually widely dispersed and mislead the tree.
# This parameter is used for filter these clusters. Greater minNumInClusterForInitTree
# would filter more noises and also may ignore some details. Default is 10.
minNumInClusterForInitTree 10


# kForGeneticMarker: Used for identifying genetic markers, greater k may helpful for
# identifying more specific markers, but also brings instability for those branches 
# with fewer LTR-RTs. Default is 2.
kForGeneticMarker 2


# onlyDetectFinalBranch: Used for identifying genetic markers. If few markers were
# identified, please set it to False to identify more potential markers. Default is
# True.
onlyDetectFinalBranch True

# Parameters used in ElPiGraph
epgLambda 0.2
epgMu 0.02
epgAlpha 0.002
```

### Usage

### Outputs
#### workDir/figure/finalResult.gif
The 3-D reconstructed trajectories of LTR-RTs.
Each dot represent a pattern sequence that extracted from one or severl
LTR-RTs. The color of the dot represents the lineage type which annotated
at the right side. The black lines represent the reconstructed trajectories.
<div align=center>
    <img src="https://github.com/xjtu-omics/LTR_Stream/blob/main/.readMe_images/e4d6ec59.png" width="300px" height="285px" />
</div>

#### workDir/figure/finalInfo.tab
Tab seperated file containing the details of `finalResult.gif`, including
Name of LTR-RT, position in `finalResult.gif`, ID of branch, lineage from
TESorter and insertion time from LTR_Retriever.
<div align=center>
    <img src="https://github.com/xjtu-omics/LTR_Stream/blob/main/.readMe_images/1a11f7df.png" width="460px" height="80px" />
</div>

#### workDir/figure/tesorter.sta.csv
Number of LTR-RTs of each lineage in each species.
<div align=center>
    <img src="https://github.com/xjtu-omics/LTR_Stream/blob/main/.readMe_images/c3602f28.png" width="500px" height="224px" />
</div>

#### workDir/geneticMarkers
##### kmeans_{k}.circos.png
Genomic distribution of identified genetic markers.
<div align=center>
    <img src="https://github.com/xjtu-omics/LTR_Stream/blob/main/.readMe_images/geneticMarker.png" width="400px" height="400px" />
</div>

##### kmeans_{k}.3D.gif
Identified genetic markers on the reconstructed evolutionary trajectories.
<div align=center>
    <img src="https://github.com/xjtu-omics/LTR_Stream/blob/main/.readMe_images/035d4f7e.png" width="300px" height="290px" />
</div>

##### kmeans_{k}.summary.tsv
Summary of the identified genetic markers.
<div align=center>
    <img src="https://github.com/xjtu-omics/LTR_Stream/blob/main/.readMe_images/5cb3625d.png" width="190px" height="74px" />
</div>

##### kmeans_{k}.oriId.tsv
LTR-RT level genetic marker annotation.  
<div align=center>
    <img src="https://github.com/xjtu-omics/LTR_Stream/blob/main/.readMe_images/573a345a.png" width="184px" height="112px" />
</div>

##### kmeans_{k}.enrichRegion.bed
Bed formatted file that annotates genomic regions with enriched genetic markers.
<!--
#### 2. Distribution of subtypes LTR-RT among genome.
```shell
conda activate ltrStream
cd path_you_want_to_install_LTR_Stream/LTR_Stream
python circosPlot.py -f path_of_ltrPara.Tab -b S1_S2,k2 -b S3_S4,l2 -c S1,k2 -c S2,l2 -o {outputPrefix}
# '-b' means plot the positions of LTR-RTs that on one specific branch.
# In case that the branch is too long, we provide an optional parameter
# 'k' to divided LTR-RTs on this branch into k clusters (using k-means).
# '-c' means plot the positions of LTR-RTs that on branches that related
# to one specific node (if the node is a leaf node, then it is equivalent
# to to set this branch in -b). You can also use 'k' parameters to divide
# the LTR-RTs into k clusters.
# Except for the 'k' parameter, we also provide another optional parameter
# 'l'. It helps to filter out LTR-RTs with simple pattern sequences that
# unlike to be genome- or subgenome-specific molecular markers.
```
-->
