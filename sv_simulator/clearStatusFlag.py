import sys
import os
workDir = sys.argv[1]

toRmFlags = ['selectNovelUnits', 'getNovelModelSeq', 'staNovelSelectedNumVsCovered',
            'mergeModules', 'getTeclassEachModId', 'calModSeqDistancePre',
            'calModSeqDistance', 'pcoaPlot', 'stream']
for flag in toRmFlags:
    os.system(f'''
        rm -rf {workDir}/status/{flag}.ok
    ''')