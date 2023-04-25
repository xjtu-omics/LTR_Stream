import argparse
import sys
import yaml

with open('envConfig.yaml', 'r') as of:
    ttPara = yaml.safe_load(of)
    sys.path.append(ttPara['LTR_Stream'])

from kmeansBranchLtrDistManager import kmeansBranchLtrDistManager


class interface_detectLtrDistribution:
    def __init__(self):
        self.args = self.init_parserPara()
        self.checkPara()
        self.distManager = self.init_manager()
        divideIdList, detectRelDictList, outPre = \
            self.distManager.detectSubLtrDist(onlyDetectFinalBranch=(self.args.onlyDetectFinalBranch=='True'))

        self.distManager.plotResults(divideIdList=divideIdList,
                                     detectRelDictList=detectRelDictList,
                                     outPre=outPre)


    def checkPara(self):
        if self.args.outDir is None:
            raise AssertionError('Please set --outDir (-d)!')
        if self.args.ltrParaFile is None:
            raise AssertionError('Please set --ltrParaFile (-c)!')
        if self.args.k < 2:
            raise AssertionError('Please set k for kmeans under kmeans mode!')
        '''
        if self.args.mode == 'kmeans':
            if self.args.k < 0:
                raise AssertionError('Please set k for kmeans under kmeans mode!')
        elif self.args.mode == 'branch':
            pass
        else:
            raise AssertionError('Please set detection mode (branch or kmeans)!')
        '''

    def init_parserPara(self):
        pyDescription = 'Detect specific genomic distribution of LTR-RTs.'
        parser = argparse.ArgumentParser(description=pyDescription)
        # mode is not provided currently.
        # parser.add_argument('-m', '--mode',
        #                     choices=('branch', 'kmeans'),
        #                     help='Mode for detecting genomic distribution of LTR-RTs.')
        parser.add_argument('-k', type=int,
                            default=-1,
                            help='K for kmeans.')
        parser.add_argument('-d', '--outDir', type=str,
                            help='Directory for storing analyse results.')
        parser.add_argument('-c', '--ltrParaFile', type=str,
                            help='ltrParaFile used in LTR_Stream.')
        parser.add_argument('-f', '--onlyDetectFinalBranch', type=str,
                            choices=['True', 'False'],
                            help='True if only identify with finalBranch. When few markers detected, please set to False to identify more markers.')
        return parser.parse_args()

    def init_manager(self):
        return kmeansBranchLtrDistManager(ltrParaFile=self.args.ltrParaFile,
                                          outDir=self.args.outDir,
                                          k=self.args.k)

        # This mode is not provided currently.
        """
        if self.args.mode == 'kmeans':
            return kmeansLtrDistManager(ltrParaFile=self.args.ltrParaFile,
                                        outDir=self.args.outDir,
                                        k=self.args.k)
        else:
            return branchLtrDistManager(ltrParaFile=self.args.ltrParaFile,
                                     outDir=self.args.outDir)
        """

detectLtrDistribution = interface_detectLtrDistribution()
