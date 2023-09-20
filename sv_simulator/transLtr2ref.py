import sys
import pysam
inputFa = sys.argv[1]
outputFa = sys.argv[2]

with pysam.FastxFile(inputFa) as inputFa, open(outputFa, 'w') as of:
    for te in inputFa:
        print('>', te.name.split(':')[0], sep='', file=of)
        print(te.sequence, file=of)