import os
import re
import time
def ttSystem(cmd):
    if os.system(cmd)!=0:
        print(f'{cmd} error!')
        exit()
def constructEnv(dirList):
    for ddd in dirList:
        if not os.path.lexists(ddd):
            ttSystem(f'mkdir -p {ddd}')
def passList2bed(passListF,outF):
    with open(outF,'w') as of:
        for line in open(passListF, 'r'):
            line = line.strip()
            if line[0] == '#':
                continue
            strand = line.split('\t')[8]
            if strand=='?':
                strand = '+'
            line = line.split('\t')[0]
            chr = line.split(':')[0]
            line = line.split(':')[1]
            st = int(line.split('..')[0])
            ed = int(line.split('..')[1])
            print(chr, st, ed, '.', '.', strand, sep='\t',file=of)
class bainfo:
    def __init__(self, _chr, _st=-1, _ed=-1, _trans='', _strand='+'):
        if _chr is None:
            self.chr = _chr
        elif ':' not in _chr:
            self.chr = _chr
            self.st = _st
            self.ed = _ed
            self.trans = _trans
            self.strand = _strand
        else:
            self.chr = _chr.split(':')[0]
            self.st = int(_chr.split(':')[1].split('-')[0])
            self.ed = int(_chr.split(':')[1].split('-')[1].split('(')[0])
            self.strand = _chr.split('(')[1][0]
            self.trans = 'nope'

    def __eq__(self, other):
        return self.chr == other.chr and self.st == other.st

    def __lt__(self, other):
        if self.chr != other.chr:
            if len(re.findall('chr([0-9]+)$', self.chr)) > 0 and \
                    len(re.findall('chr([0-9]+)$', other.chr)) > 0:
                myChrId = int(re.findall('chr([0-9]+)$', self.chr)[0])
                otherId = int(re.findall('chr([0-9]+)$', other.chr)[0])
                return myChrId < otherId
            elif len(re.findall('([0-9]+)$', self.chr)) > 0 and \
                    len(re.findall('([0-9]+)$', other.chr)) > 0:
                myChrId = int(re.findall('([0-9]+)$', self.chr)[0])
                otherId = int(re.findall('([0-9]+)$', other.chr)[0])
                return myChrId < otherId
            else:
                return self.chr < other.chr
        if self.st != other.st:
            return self.st < other.st
        else:
            return self.ed < other.ed

    def __cmp__(self, other):
        if self.chr < other.chr:
            return -1
        elif self.chr > other.chr:
            return 1
        elif self.st < other.st:
            return -1
        elif self.st > other.st:
            return 1
        elif self.st == other.st:
            return 0
        else:
            return

    def print(self, of=None):
        if of == None:
            print(self.chr, self.st, self.ed, self.trans, '.', self.strand, sep='\t')
        else:
            print(self.chr, self.st, self.ed, self.trans, '.', self.strand, sep='\t', file=of)

    def getStrBainfo(self):
        return f'{self.chr}:{self.st}-{self.ed}({self.strand})'
def overlapMax(a,b):
    if a.chr != b.chr:
        return 0
    overDis=max(min(a.ed,b.ed)-max(a.st,b.st),0)
    return max(overDis/(a.ed-a.st+1),overDis/(b.ed-b.st+1))
def overlapMin(a,b):
	overDis=max(min(a.ed,b.ed)-max(a.st,b.st),0)
	return min(overDis/(a.ed-a.st+1),overDis/(b.ed-b.st+1))
def overlapA(a,b):
    if a.chr!=b.chr:
        return 0
    overDis=max(min(a.ed,b.ed)-max(a.st,b.st),0)
    return overDis/(a.ed-a.st+1)
def mergeBainfo(a,b,checkStrand=False):
    if a.chr != b.chr or overlapMax(a,b) == 0 or (checkStrand and a.strand != b.strand):
        print(f'Merge false of {a.getStrBainfo()} and {b.getStrBainfo()}')
        exit()
    else:
        nst = min(a.st,b.st)
        ned = max(a.ed,b.ed)
        return bainfo(a.chr,nst,ned,a.trans,a.strand)
def getRestRange(totalRange,selectedList):
    i = 0
    while i < len(selectedList):
        mergeFlag = False
        for j in range(i+1,len(selectedList)):
            if overlapMax(selectedList[i],selectedList[j])>0:
                selectedList[i] = mergeBainfo(selectedList[i],selectedList[j])
                selectedList[j] = selectedList[-1]
                selectedList.pop()
                mergeFlag = True
                break
        if not mergeFlag:
            i = i+1
        else:
            i = 0
    selectedList.sort()
    restList = []
    if selectedList[0].st>totalRange.st:
        restList.append(bainfo(_chr=totalRange.chr,_st=totalRange.st,_ed=selectedList[0].st-1,_strand=totalRange.strand))
    for i in range(len(selectedList)-1):
        nowb = selectedList[i]
        nextb = selectedList[i+1]
        restList.append(bainfo(_chr=nowb.chr,_st=nowb.ed+1,_ed=nextb.st-1,_strand=totalRange.strand))
    if selectedList[-1].ed < totalRange.ed:
        restList.append(bainfo(_chr=totalRange.chr,_st=selectedList[-1].ed+1,_ed=totalRange.ed,_strand=totalRange.strand))
    return restList
def bainfoList2BedFile(infoList,bedFile,writeTyp='w'):
    with open(bedFile,writeTyp) as of:
        for bai in infoList:
            print(bai.chr, bai.st, bai.ed,'.','.',bai.strand ,sep='\t', file=of)
class ttUFS(object):
    def __init__(self, setSize):
        self.setSize = setSize
        self.fa = []
        self.si = []
        for i in range(setSize):
            self.fa.append(i)
            self.si.append(1)
    def find_set(self, x):
        stPoint = x
        while x!=self.fa[x]:
            x = self.fa[x]
        endFa = x
        x = stPoint
        while x!=self.fa[x]:
            tmpX = x
            x =self.fa[x]
            self.fa[tmpX] = endFa
        return endFa

    def find_size(self,x):
        x = self.find_set(x)
        return self.si[x]

    def union(self, a, b):
        a_father = self.find_set(a)
        b_father = self.find_set(b)
        if(a_father != b_father):
            a_size = self.si[a_father]
            b_size = self.si[b_father]
            if(a_size >= b_size):
                self.fa[b_father] = a_father
                self.si[a_father] = a_size + b_size
                self.si[b_father] = 0
            else:
                self.fa[a_father] = b_father
                self.si[b_father] = a_size + b_size
                self.si[a_father] = 0

    def staNum(self,outFile):
        from collections import defaultdict
        ufsId2num = defaultdict(int)
        for i in range(self.setSize):
            ufsId2num[self.find_set(i)] = ufsId2num[self.find_set(i)]+1
        with open(outFile,'w') as of:
            for ufsId in ufsId2num:
                print(ufsId,ufsId2num[ufsId],file=of)

    def getRotIdSet(self):
        from collections import defaultdict
        relSet = defaultdict(int)
        for i in range(self.setSize):
            relSet[self.find_set(i)] = self.find_size(i)
        return relSet
def plotRotate3DScatter(totalData,outFile):
    import matplotlib as mpl
    import pandas as pd
    mpl.use('agg')
    mpl.rcParams['pdf.fonttype'] = 42
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import animation
    fig = plt.figure(1,figsize=(14,12))
    ax = plt.axes(projection='3d')
    x = totalData.iloc[:,0]
    y = totalData.iloc[:,1]
    z = totalData.iloc[:,2]
    scalex = 1.0/(x.max()-x.min())
    scaley = 1.0/(y.max()-y.min())
    scalez = 1.0/(z.max()-z.min())
    def init():
        ax.scatter3D(x*scalex,y*scaley,z*scalez,linewidth=0.1)
        ax.set_xlabel('$TSNE1$', fontdict={'fontsize': 24})
        ax.set_ylabel('$TSNE2$', fontdict={'fontsize': 24})
        ax.set_zlabel('$TSNE3$', fontdict={'fontsize': 24})
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        plt.title('')
        plt.grid()
    def rotate(angle):
        ax.view_init(elev=10,azim=angle)
    rot_animation = animation.FuncAnimation(fig,rotate,init_func=init,frames=np.arange(0,362,2),interval=50)
    rot_animation.save(outFile,dpi=80,writer='imagemagick')
def plotRotate3DScatterWithColor(totalData,outFile,colorTyp='num',colorDict = None):
    import matplotlib as mpl
    import pandas as pd
    mpl.use('agg')
    mpl.rcParams['pdf.fonttype'] = 42
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import animation
    fiveColors = ['#cdcdcd', '#66CC96', '#8AA0CF', '#D1B36D', '#BB65E0', '#F570A3', '#000000']
    if colorTyp=='num':
        from matplotlib import animation
        fig = plt.figure(1,figsize=(14,12))
        ax = plt.axes(projection='3d')
        x = totalData.iloc[:,0]
        y = totalData.iloc[:,1]
        z = totalData.iloc[:,2]
        colorId = totalData.iloc[:,3]
        colorId = colorId % 7
        totColors = np.array(fiveColors)
        colorId = totColors[colorId]
        def init():
            ax.scatter3D(x,y,z,c=colorId,linewidth=0.1)
            ax.set_xlabel('$Dim1$')
            ax.set_ylabel('$Dim2$')
            ax.set_zlabel('$Dim3$')
            plt.title('Compositional Plot in 3 Dimensions')
            plt.grid()
        def rotate(angle):
            ax.view_init(elev=10,azim=angle)
        rot_animation = animation.FuncAnimation(fig,rotate,init_func=init,frames=np.arange(0,362,2),interval=50)
        rot_animation.save(outFile,dpi=80,writer='imagemagick')
    elif colorTyp=='str':
        from matplotlib import animation
        from collections import defaultdict
        fig = plt.figure(1, figsize=(14, 12))
        ax = plt.axes(projection='3d')
        x = totalData.iloc[:, 0]
        y = totalData.iloc[:, 1]
        z = totalData.iloc[:, 2]
        oriStr2colorId = defaultdict(int)
        idUsed = 0
        colorId = []
        if colorDict is not None:
            for ind,row in totalData.iterrows():
                colorId.append(colorDict[row[3]])
            colorId = np.array(colorId)
        else:
            for ind,row in totalData.iterrows():
                if row[3] not in oriStr2colorId:
                    oriStr2colorId[row[3]] = idUsed
                    idUsed = idUsed+1
                    print(row[3],fiveColors[oriStr2colorId[row[3]]%7])
                colorId.append(oriStr2colorId[row[3]])
            colorId = np.array(colorId)
            colorId = colorId % 7
            totColors = np.array(fiveColors)
            colorId = totColors[colorId]

        def init():
            ax.scatter3D(x, y, z, c=colorId, linewidth=0.1,label=colorId)
            ax.set_xlabel('$PCoA1$',fontdict={'fontsize': 24})
            ax.set_ylabel('$PCoA2$',fontdict={'fontsize': 24})
            ax.set_zlabel('$PCoA3$',fontdict={'fontsize': 24})
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_zticklabels([])
            plt.title('')
            plt.grid()

        def rotate(angle):
            ax.view_init(elev=10, azim=angle)

        rot_animation = animation.FuncAnimation(fig, rotate, init_func=init, frames=np.arange(0, 362, 2), interval=50)
        rot_animation.save(outFile, dpi=80, writer='imagemagick')
def qsub(baseD, commands, taskName, threads, pbsFileName='tmp.pbs', nodes=1):
    pbsF = f'{baseD}/{pbsFileName}'
    with open(pbsF, 'w') as of:
        print(commands, file=of)
    os.system(f'''
        cd {baseD}
        qsub {pbsF} -N {taskName} \
                    -l walltime=9999:00:00 \
                    -l nodes={nodes}:ppn={threads} \
    ''')
    time.sleep(2)

def plotDisVsModNum(tpdf, outPdf):
    from copy import deepcopy
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.stats import spearmanr
    sns.set_style('ticks')
    tpdf = deepcopy(tpdf)

    tpdf = tpdf.dropna()

    tpdf['dis'] = tpdf.apply(
        lambda row: np.linalg.norm(np.array(row[['x', 'y', 'z']])),
        axis=1
    )
    corr, p_value = spearmanr(tpdf.dis, tpdf.num)
    fig, ax = plt.subplots()
    sns.scatterplot(
        tpdf,
        x='dis',
        y='num',
        linewidth=0,
        ax=ax,
        s=8
    )
    ax.text(
        .1, .9,
        f'Spearman: {corr:.2f}',
        fontsize=18,
        transform=plt.gca().transAxes
    )
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    ax.set_xlabel('Distance to center', fontdict={'fontsize': 18})
    ax.set_ylabel('Modular sequence length', fontdict={'fontsize': 18})
    ax.spines[['top', 'right']].set_visible(False)
    plt.savefig(outPdf, bbox_inches='tight')

























