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

import sys
import os
from ttUtils import getOriId2tesorterClass
from ttUtils import seqName2bainfo

baseD = '/data/home/testXT/workData/simulateMeiSv'

def plotLenDensity():
    # Tekay of rice s001 was selected.
    oriId2teClassTsv = '/data/home/testXT/testLTR_Stream/rice/workDir/tesorter/s001.tesorter.out.cls.tsv'
    oriId2teClass = getOriId2tesorterClass(oriId2teClassTsv)
    lenList = []
    for oriId in oriId2teClass:
        if oriId2teClass[oriId] == 'Tekay':
            bainfo = seqName2bainfo(oriId)
            lenList.append(bainfo.ed-bainfo.st)
    sns.histplot(x=lenList, bins=20)
    plt.xlabel('LTR-RT Length', fontsize=16)
    plt.ylabel('Count', fontsize=16)
    plt.tick_params(axis='both', labelsize=12)
    # plt.show()
    plt.savefig(f'{baseD}/dist_tekay_length.pdf')
    plt.close()

def selectLTR_RT():
    oriId2teClassTsv = '/data/home/testXT/testLTR_Stream/rice/workDir/tesorter/s001.tesorter.out.cls.tsv'
    oriId2teClass = getOriId2tesorterClass(oriId2teClassTsv)
    refFaFile = '/data/home/testXT/testLTR_Stream/rice/workDir/ref/s001.ref.fa'
    dataD = f'{baseD}/selectLTR_RT'
    selectedFastaFile = f'{dataD}/selected.fasta'
    alignedFastaFile = f'{dataD}/aligned.fasta'
    refFa = pysam.FastaFile(refFaFile)

    """
    tarLenDict = dict(_8k=(8000, 9000),
                      _12k=(11800, 12200),
                      _15k=(14500, 15500))
    tar2oriId = dict()
    for oriId in oriId2teClass:
        if oriId2teClass[oriId] == 'Tekay':
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

    os.system(f'''
        mkdir -p {dataD}
    ''')
    with open(selectedFastaFile, 'w') as of:
        for k in tar2seq:
            print(f'>{k}', file=of)
            print(str(tar2seq[k]), file=of)

    os.system(f'''
        PATH=/data/home/testXT/miniconda3/envs/ltrStream/bin:$PATH
        muscle -align {selectedFastaFile} -output {alignedFastaFile}
    ''')
    """
    alignment = AlignIO.read(alignedFastaFile, 'fasta')
    summary_align = AlignInfo.SummaryInfo(alignment)
    print(summary_align.dumb_consensus(1))









# plotLenDensity()
selectLTR_RT()