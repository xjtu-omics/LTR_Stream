import pandas as pd
from collections import defaultdict
def getParaDict(configFile):
    paraDict = defaultdict(int)
    paraD = pd.read_table(configFile,header=None,sep='\t',comment='#')
    for ind,row in paraD.iterrows():
        paraDict[row[0]] = row[1]
    return paraDict