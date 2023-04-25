import pandas as pd
from collections import defaultdict
import re
import ttUtils

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
def getOriId2inRangeId(configFile,ltrReD):
    rel = defaultdict(str)
    for spe in getSpeIdList(configFile):
        rel = ttUtils.getOriId2inRangeId(f'{ltrReD}/{spe}/{spe}.ref.fa.pass.list',relDict=rel)
    return rel
def getOriId2tesorterClass(configFile,tesorterD):
    rel = defaultdict(int)
    for spe in getSpeIdList(configFile):
        tesorterOutFile = f'{tesorterD}/{spe}.tesorter.out.cls.tsv'
        rel = ttUtils.getOriId2tesorterClass(tesorterOutFile,rel)
    return rel
def getOriId2insertTime(configFile,ltrRetrieverD):
    rel = defaultdict(str)
    for spe in getSpeIdList(configFile):
        rel = ttUtils.getOriId2insertTime(f'{ltrRetrieverD}/{spe}/{spe}.ref.fa.pass.list',preDict=rel)
    return rel
def getSpeId2speName(configFile):
    relDict = defaultdict(str)
    speInfoD = pd.read_table(configFile, sep='\t',comment='#')
    for ind, row in speInfoD.iterrows():
        relDict[row['speId']] = row['speName']
    return relDict