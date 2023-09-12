import argparse
import pysam
from svSimulator import svSimulator

parser = argparse.ArgumentParser(description='SV simulator for LTR_Stream')
parser.add_argument('--maxPopSize', '-mps', type=int)
parser.add_argument('--seqFa', '-sfa', type=str)
args = vars(parser.parse_args())

maxPopSize = 5000 if args['maxPopSize'] is None else args['maxPopSize']
seqFa = ''

svSim = svSimulator(maxPopSize=maxPopSize)
