# Create a conda environment named ltrStream for LTR_Stream
condaPath=`which conda`
condaPath=${condaPath%/*}
envPath=${condaPath}/../envs
conda=conda
if [ $# == 1 ]; then
  conda=$1
fi
$conda env create -f config/ltrStream.yml

# Create another conda environment for ltrHarvest
#ltrHarvestEnvName='ltrHarvestEnv'
#$conda create -n $ltrHarvestEnvName -c bioconda -c conda-forge python=2.7 -y
#$conda install -n $ltrHarvestEnvName -c bioconda -c conda-forge genometools genometools-genometools -y
#ltrHarvestEnvPath=${envPath}/${ltrHarvestEnvName}/bin

# Install LTR_FINDER_parallel
#ltrFinder=ltrFinder
#mkdir -p $ltrFinder
#cd $ltrFinder
#git clone git@github.com:oushujun/LTR_FINDER_parallel
#cd ..
#ltrFinderEnvPath=`pwd`/$ltrFinder/LTR_FINDER_parallel

# Install dante
#dante=dante
#mkdir -p $dante
#cd $dante
#git clone git@github.com:kavonrtep/dante.git
#cd ..
#danteEnvPath=`pwd`/$dante/dante

# Genenrate configYaml
configYaml='config/envConfig.yaml'
echo '---' > $configYaml
echo 'LTR_Stream: ' \"`pwd`/src\" >> $configYaml
#echo 'danteEnv: ' \"${danteEnvPath}\" >> $configYaml
#echo 'ltrHarvestEnv: ' \"${ltrHarvestEnvPath}\" >> $configYaml
#echo 'ltrFinderEnv: ' \"${ltrFinderEnvPath}\" >> $configYaml
