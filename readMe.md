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
    <img src="https://github.com/xjtu-omics/LTR_Stream/blob/main/.readMe_images/Papaver_Retand_all.gif" width="700px" height="592px" /> 
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
cd ${ltrStreamInstallPath} && git clone xxx
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
#### 1. LTR-RT evolutionary trajectories reconstruction and classification
```shell
conda activate ltrStream
cd ${ltrStreamInstallPath}/LTR_Stream
snakemake -s LTR_Stream.smk -f stream --config ltrParaFile=path_of_ltrPara.tsv -j {threadsNumber}
```
#### 2. Genetic maker detection
```shell
conda activate ltrStream
cd ${ltrStreamInstallPath}/LTR_Stream
snakemake -s LTR_Stream.smk -f geneticMarkerDetect --config ltrParaFile=path_of_ltrPara.tsv -j {threadsNumber}
```

### Config files
#### 1.`ltrPara.tsv`
Two config files should be prepared for running LTR_Stream.  
The first file is called `ltrPara.tsv`(Of course you can 
rename the name of the file, we use `ltrPara.tsv` here to 
refer to this config file). The second file is called `ref.tsv`. \
For `ltrPara.tsv`, two mandatory parameters are needed:
```tsv
# An example for ltrPara.tsv
# All tab seperated.

# Mandatory parameters
# workDir: A blank directory for running LTR_Stream
# The outputs of LTR_Stream are in workDir/figure
workDir /data/home/xutun/mei/rice/workDir
# refConfig: Path of the second config file where setting the information of genome assemblies.
refConfig /data/home/xutun/project/LTR_Stream/config/ref.tsv



# Optional parameters
# transAllContigs: LTR_Stream usually extracts main chromosomes for LTR-RT  
# identification. If you wish to include those LTR-RTs from contigs, please
# set transAllContigs to True. Default is False.
transAllContigs False


# minOverLapForNovelModule: Control the dispersion of CEPs in the 3-D space.
# It is used in disjoint-set data structure to judge if there should be an
# edge between two alignment regions. It could be set at the range from 0 to 1. 
# Greater minOverLapForNovelModule leads to more dispersed result. Default is 0.8.
minOverLapForNovelModule 0.8


# maxCNSPercentage: Control the dispersion of CEPs with minOverLapForNovelModule.
# Greater maxCNSPercentage leads to more dispersed result. It should be set at 
# a range from 0 to 1. maxCNSPercentage directly control the number of most frequent
# novel modules used for reconstructing evolutionary trajectories. Default is 0.7.
maxCNSPercentage 0.7


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
#### 2.`ref.tsv` is a tsv file recording the information of genome assemblies.
```tsv
# An example for ref.tsv
speId speName chrNum miu source
s001 Oryza sativa_indica 12 1.3e-8  https://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Oryza_sativa/all_assembly_versions/GCA_001889745.1_Rice_IR8_v1.7/GCA_001889745.1_Rice_IR8_v1.7_genomic.fna.gz
s002 Oryza alta 24 1.3e-8 /data/home/xutun/mei/rice/resources/originRef/GWHAZTO00000000.genome.fasta.gz
s003 Oryza punctata 12 1.3e-8 https://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Oryza_punctata/latest_assembly_versions/GCA_000573905.2_OpunRS2/GCA_000573905.2_OpunRS2_genomic.fna.gz
s004 Oryza brachyantha 12 1.3e-8 https://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Oryza_brachyantha/latest_assembly_versions/GCA_000231095.3_ObraRS2/GCA_000231095.3_ObraRS2_genomic.fna.gz
# speId: An ID for species, should be started with an English letter, where only English
# letters and numbers are permitted.
# speName: The Latin name of the species.
# chrNum: Number of chromosomes used in the analysis, LTR_Stream automatically extracts 
# the first chrNum sequences in the fasta file.
# miu: Neutral mutation rate that used to estimate the insertion time of LTR-RTs.
# source: The source of genome assemblies. It can be a download link or a path for local
# fasta file.
# Can be fasta format or a compressed format (.gz).
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
