import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
import matplotlib.pyplot as plt
import seaborn as sns

import Bio
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Align import AlignInfo
import pysam
import numpy as np



import sys
import yaml
import os
from ttUtils import getOriId2tesorterClass
from ttUtils import seqName2bainfo

baseD = '/data/home/testXT/workData/simulateMeiSv'


def plotLenDensity(group, speId, lin):
    # Tekay of rice s001 was selected.
    # group: rice
    # speId: s001
    oriId2teClassTsv = f'/data/home/testXT/testLTR_Stream/{group}/workDir/tesorter/{speId}.tesorter.out.cls.tsv'
    # lin = 'Tekay'
    oriId2teClass = getOriId2tesorterClass(oriId2teClassTsv)
    lenList = []
    for oriId in oriId2teClass:
        if oriId2teClass[oriId] == lin:
            bainfo = seqName2bainfo(oriId)
            lenList.append(bainfo.ed-bainfo.st)
    sns.histplot(x=lenList, bins=40)
    plt.xlabel('LTR-RT Length', fontsize=16)
    plt.ylabel('Count', fontsize=16)
    plt.tick_params(axis='both', labelsize=12)
    # plt.show()
    plt.savefig(f'{baseD}/{group}_{speId}_{lin}_lengthDist.pdf')
    plt.close()

def selectLTR_RT(configYaml):
    with open(configYaml, 'r') as ym:
        paraDict = yaml.safe_load(ym)
    group = paraDict['group']
    speId = paraDict['speId']
    lin = paraDict['lin']
    workDir = paraDict['workDir']
    smlRange = paraDict['smlRange']
    midRange = paraDict['midRange']
    lrgRange = paraDict['lrgRange']

    os.system(f'''
        mkdir -p {workDir}
    ''')

    oriId2teClassTsv = f'/data/home/testXT/testLTR_Stream/{group}/workDir/tesorter/{speId}.tesorter.out.cls.tsv'
    oriId2teClass = getOriId2tesorterClass(oriId2teClassTsv)
    refFaFile = f'/data/home/testXT/testLTR_Stream/{group}/workDir/ref/{speId}.ref.fa'
    selectedFastaFile = f'{workDir}/selected.fasta'
    alignedFastaFile = f'{workDir}/aligned.fasta'
    consensusFastaFile = f'{workDir}/consensus.fasta'
    refFa = pysam.FastaFile(refFaFile)

    # """
    # Generate aligned fasta
    tarLenDict = dict(sml=(smlRange[0], smlRange[1]),
                      mid=(midRange[0], midRange[1]),
                      lrg=(lrgRange[0], lrgRange[1]))
    tar2oriId = dict()
    for oriId in oriId2teClass:
        if oriId2teClass[oriId] == lin:
            bainfo = seqName2bainfo(oriId)
            for k in tarLenDict:
                if (tarLenDict[k][0] <= bainfo.ed-bainfo.st <= tarLenDict[k][1]) \
                   and (k not in tar2oriId):
                    tar2oriId[k] = oriId

    tar2seq = dict()
    for k in tar2oriId:
        bainfo = seqName2bainfo(tar2oriId[k])
        seq = Seq(refFa.fetch(bainfo.chr, bainfo.st, bainfo.ed))
        if bainfo.strand=='-':
            seq = seq.reverse_complement()
        tar2seq[k] = seq
    refFa.close()

    with open(selectedFastaFile, 'w') as of:
        for k in tar2seq:
            print(f'>{k}', file=of)
            print(str(tar2seq[k]), file=of)

    os.system(f'''
        PATH=/data/home/testXT/miniconda3/envs/ltrStream/bin:$PATH
        muscle -align {selectedFastaFile} -output {alignedFastaFile}
    ''')
    # """

    alignment = AlignIO.read(alignedFastaFile, 'fasta')
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensusSeq = str(summary_align.dumb_consensus(1)).upper()
    consensusSeq = consensusSeq.split('X')
    consensusSeq = ''.join(consensusSeq)
    with open(consensusFastaFile, 'w') as of:
        print('>consensus', file=of)
        print(consensusSeq, file=of)

def generateSimulatedFasta(configYaml):
    from svSimulator import svSimulator
    svGenTyp = 'sub'
    # selectLTR_RT()
    with open(configYaml, 'r') as ym:
        paraDict = yaml.safe_load(ym)

    workDir = paraDict['workDir']
    minSvLen = paraDict['minSvLen']
    maxSvLen = paraDict['maxSvLen']
    cpuNum = paraDict['cpuNum']
    randomSeed = paraDict['randomSeed']

    maxPopSize=1000
    ranSvP = 1e-9
    yearStep = 3e5
    cpP = 4e-6
    np.random.seed(randomSeed)

    selectedFastaFile = f'{workDir}/selected.fasta'
    consensusFastaFile = f'{workDir}/consensus.fasta'

    consensusFa = pysam.FastaFile(consensusFastaFile)
    selectedFa = pysam.FastxFile(selectedFastaFile)
    consensusSeq = consensusFa.fetch('consensus')
    with open(f'{workDir}/simulatedTe.fa', 'w') as of:
        for i, read in enumerate(selectedFa):
            tarSeq = read.sequence.upper()
            oriSeq = consensusSeq.upper()
            mySimer = svSimulator(maxPopSize=maxPopSize,
                                  oriSeq=oriSeq,
                                  oriSeqName=f'{i}',
                                  tarSeq=tarSeq,
                                  ranSvP=ranSvP,
                                  yearStep=yearStep,
                                  cpP = cpP,
                                  cpuNum=cpuNum,
                                  minSvLen=minSvLen,
                                  maxSvLen=maxSvLen)
            mySimer.simulate(svGenTyp=svGenTyp)
            mySimer.meiPop.printFasta(year=None, of=of)
        # break

    consensusFa.close()
    selectedFa.close()

