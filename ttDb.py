import pandas as pd
from collections import defaultdict
import re
import ttUtils
import pysam

def getSpeIdList(configFile):
    speInfoD = pd.read_table(configFile,sep='\t',comment='#')
    return pd.Series(speInfoD['speId']).tolist()
def getSpeId2dlWeb(configFile):
    relDict = defaultdict(str)
    speInfoD = pd.read_table(configFile,sep='\t',comment='#')
    for ind,row in speInfoD.iterrows():
        relDict[row['speId']] = row['source']
    return relDict
def getSpeId2chrNum(configFile):
    relDict = defaultdict(str)
    speInfoD = pd.read_table(configFile,sep='\t',comment='#')
    for ind,row in speInfoD.iterrows():
        relDict[row['speId']] = row['chrNum']
    return relDict
def getSpeId2miu(configFile):
    relDict = defaultdict(str)
    speInfoD = pd.read_table(configFile,sep='\t',comment='#')
    for ind,row in speInfoD.iterrows():
        relDict[row['speId']] = row['miu']
    return relDict
def getOriId2inRangeId(configFile,ltrReD,ltrFasta=None):
    # None of configFile indicates Fasta model, which does not support inRangeId
    # (no terminal LTR) ranges identification.
    if configFile is None:
        relDict = {}
        with pysam.FastxFile(ltrFasta) as fh:
            for ltr in fh:
                relDict[ltr.name] = ltr.name
        return relDict
    rel = defaultdict(str)
    for spe in getSpeIdList(configFile):
        rel = ttUtils.getOriId2inRangeId(f'{ltrReD}/{spe}/{spe}.ref.fa.pass.list',relDict=rel)
    return rel
def getOriId2tesorterClass(configFile,tesorterD,tesorterOutFile):
    rel = defaultdict(int)
    if configFile is None:
        return ttUtils.getOriId2tesorterClass(tesorterOutFile, rel)

    for spe in getSpeIdList(configFile):
        tesorterOutFile = f'{tesorterD}/{spe}.tesorter.out.cls.tsv'
        rel = ttUtils.getOriId2tesorterClass(tesorterOutFile,rel)
    return rel
def getOriId2insertTime(configFile,ltrRetrieverD, ltrFasta):
    rel = defaultdict(str)
    if configFile is None:
        with pysam.FastxFile(ltrFasta) as inFa:
            for te in inFa:
                rel[te.name] = -1
        return rel

    for spe in getSpeIdList(configFile):
        rel = ttUtils.getOriId2insertTime(f'{ltrRetrieverD}/{spe}/{spe}.ref.fa.pass.list',preDict=rel)
    return rel
def getSpeId2speName(configFile):
    relDict = defaultdict(str)
    speInfoD = pd.read_table(configFile, sep='\t',comment='#')
    for ind, row in speInfoD.iterrows():
        relDict[row['speId']] = row['speName']
    return relDict