import sys
from simulateBase import generateSimulatedFasta, selectLTR_RT

selectLTR_RT(sys.argv[1])
generateSimulatedFasta(sys.argv[1])