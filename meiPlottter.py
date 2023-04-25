import matplotlib as mpl
mpl.use('agg')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
import sys
import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib import animation
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
from collections import defaultdict
import pycircos
import yaml

with open('envConfig.yaml', 'r') as of:
    ttPara = yaml.safe_load(of)
    sys.path.append(ttPara['LTR_Stream'])

from ttUtils import getChrm2speId
from ttDb import getOriId2insertTime
from ttDb import getSpeId2speName
from fileSysSupplier import fileSysSupplier


class plotterBase:
    def __init__(self, initSource, initTyp, colorList, ltrParaFile, speIdSet=None):
        self.sysSupplier = fileSysSupplier(ltrParaFile)
        self.speIdSet2speName = getSpeId2speName(self.sysSupplier.refParaFile)
        self.colorList = colorList
        self.init_fromSource(initSource, initTyp, speIdSet)


    def init_fromSource(self, initSource, initTyp, speIdSet):
        if initTyp == 'ltrDistManager':
            if speIdSet is None:
                self.speIdSet = {speId: 1 for speId in initSource.distDetector.species2chrNum}
            else:
                self.speIdSet = speIdSet
            self.chrm2len = defaultdict(int)
            self.sortedChrmList = []
            self.genomeLen = 0
            for chrm in initSource.distDetector.chrm2len:
                chrmLen = initSource.distDetector.chrm2len[chrm]
                if getChrm2speId(chrm) in self.speIdSet:
                    self.chrm2len[chrm] = chrmLen
                    self.genomeLen += chrmLen
            for chrm in initSource.distDetector.sortedChrmList:
                if getChrm2speId(chrm) in self.speIdSet:
                    self.sortedChrmList.append(chrm)
            self.seqName2bainfo = initSource.distDetector.seqName2bainfo


    def filterBainList(self, typList, bainList):
        filteredBainList = []
        filteredTypList = []
        for typ, bain in zip(typList, bainList):
            if bain.chr in self.chrm2len:
                filteredTypList.append(typ)
                filteredBainList.append(bain)
        return filteredTypList, filteredBainList

    def filterOriIdList(self, typList, oriIdList):
        filteredOriIdList = []
        filteredBainfoList = []
        filteredTypList = []
        for typ, oriId in zip(typList, oriIdList):
            bain = self.seqName2bainfo(oriId)
            chrm = bain.chr
            if chrm in self.chrm2len:
                filteredOriIdList.append(oriId)
                filteredBainfoList.append(bain)
                filteredTypList.append(typ)
        return filteredTypList, filteredBainfoList, filteredOriIdList

    def checkTypeList(self, lowestId, highestId, typList):
        for typ in typList:
            if typ > highestId or typ < lowestId:
                raise Exception(f'typ in typList used in meiPlotter.circosPlot must be integar from 1-7')

    def getOriId2insertTime(self):
        return getOriId2insertTime(self.sysSupplier.refParaFile, f'{self.sysSupplier.baseD}/ltrRetriever')


class circosPlotter(plotterBase):
    def __init__(self, initSource, initTyp, colorList, ltrParaFile, speIdSet=None):
        super().__init__(initSource, initTyp, colorList, ltrParaFile, speIdSet)
        self.init_circosRela(colorList)

    def init_circosRela(self, colorList):
        if not (colorList is None):
            self.circos_maxPeriBarNum = len(self.colorList)
            # the length of each LTR-RT (Too thin bar is not significant!)
            self.ltrBarLen = self.genomeLen / 2000
            # peri margin for tick labels
            self.circos_rOff = 180
            self.genomeText_rOff = 50
            self.circos_maxTickDis, self.circos_minTickDis = self.initTickDis()
            self.circos_rRangeList = self.initRRangeList()
            self.chrmPeriRange = (950 - self.circos_rOff, 1000 - self.circos_rOff)

    def initTickDis(self):
        longestChrmLen = 0
        for chrm in self.chrm2len:
            longestChrmLen = max(longestChrmLen, self.chrm2len[chrm])
        miniDis = 1e12
        i = 1000000
        while True:
            newMiniDis = abs(longestChrmLen / 4 - i)
            if miniDis > newMiniDis:
                miniDis = newMiniDis
            else:
                # The smallest distance is the previous i
                # int(???/1MB)*1MB is for removing reading difficulty of binary system
                circos_maxTickDis = int(i / 2 / 1000000 / 10) * 1000000 * 10
                circos_minTickDis = circos_maxTickDis / 5
                return int(circos_maxTickDis), int(circos_minTickDis)
            i *= 2

    def initRRangeList(self):
        # A list like [(900,950), (850,900)...]
        circos_rRangeList = []
        nowOff = 50
        for i in range(self.circos_maxPeriBarNum):
            periEd = 1000 - nowOff - self.circos_rOff
            nowOff += 50
            periSt = 1000 - nowOff - self.circos_rOff
            circos_rRangeList.append((periSt, periEd))
        return circos_rRangeList

    def setPeriTicks(self, circle):
        for chrm in self.chrm2len:
            ticklabels = [int(tl / 1000000) for tl in range(0, self.chrm2len[chrm] + 1, self.circos_maxTickDis)]
            circle.tickplot(garc_id=chrm,
                            tickinterval=self.circos_maxTickDis,
                            ticklabels=ticklabels,
                            tickdirection='outer',
                            raxis_range=(1000 - self.circos_rOff + 1, 1000 - self.circos_rOff + 10),
                            ticklabelorientation="horizontal",
                            ticklabelsize=4)
            circle.tickplot(garc_id=chrm,
                            tickinterval=self.circos_minTickDis,
                            tickdirection='outer',
                            raxis_range=(1000 - self.circos_rOff + 1, 1000 - self.circos_rOff + 5))

    def addLegend(self, typList, typId2annot):
        legendElementList = []
        typSet = {typ: 1 for typ in typList}
        for typ in typSet:
            legendElementList.append(Line2D([], [],
                                            marker='|',
                                            color=self.colorList[typ - 1],
                                            label=typId2annot[typ],
                                            markersize=15,
                                            markeredgewidth=1.5,
                                            linestyle='None'))
        plt.legend(handles=legendElementList, loc='center')

    def initFrame(self):
        # Init circos base.
        tmpColorList = ['#6e9ece', '#BDB76B', '#4e9595', '#8d6ab8',
                        '#84574d', '#454545', '#e6928f']
        circle = pycircos.Gcircle()
        tmpSpeIdSet = defaultdict(int)
        for chrm in self.sortedChrmList:
            tmpSpeIdSet[chrm.split('_')[0]] = 1
            # Facecolor should be modified to self.colorList
            arc = pycircos.Garc(arc_id=chrm,
                                size=self.chrm2len[chrm],
                                interspace=1,
                                raxis_range=self.chrmPeriRange,
                                labelposition=85,
                                label_visible=True,
                                label=chrm.split('_')[-1],
                                facecolor=tmpColorList[len(tmpSpeIdSet)-1])
            circle.add_garc(arc)
        circle.set_garcs()
        return circle

    def addGenomeText(self, circle):
        genome2len = defaultdict(int)
        genome2stEdChrList = defaultdict(int)
        totToPlotLen = 0
        for chrm in self.chrm2len:
            genome = chrm.split('_')[0]
            genome2len[genome] += self.chrm2len[chrm]
            totToPlotLen += self.chrm2len[chrm]
        sortedGenomeList = []
        for chrm in self.sortedChrmList:
            genome = chrm.split('_')[0]
            if not (genome in sortedGenomeList):
                sortedGenomeList.append(genome)
            if not (genome in genome2stEdChrList):
                genome2stEdChrList[genome] = [chrm, None]
            else:
                genome2stEdChrList[genome][1] = chrm
        for genome in sortedGenomeList:
            stChrm = genome2stEdChrList[genome][0]
            edChrm = genome2stEdChrList[genome][1]
            pos = circle._garc_dict[stChrm].coordinates[0]
            width = circle._garc_dict[edChrm].coordinates[1] - circle._garc_dict[stChrm].coordinates[0]
            genomeLabel = self.speIdSet2speName[genome]
            rot = (circle._garc_dict[stChrm].coordinates[0] + circle._garc_dict[edChrm].coordinates[1]) / 2
            rot = rot*360/(2*np.pi)
            if 90 < rot < 270:
                rot = 180 - rot
            else:
                rot = -1 * rot
            circle.ax.text(pos + width/2,
                           1000-self.genomeText_rOff,
                           genomeLabel,
                           rotation=rot,
                           ha="center",
                           va="center",
                           fontsize=circle._garc_dict[stChrm].labelsize+5,
                           fontstyle='italic')

    def addPeriSegment(self, circle, toPlotData):
        # toPlotData should be a pandas.dataframe that has four columns
        # chrm, position, width, faceColor, periRange

        # Each peri bar that belongs to each chrm should be plotted seperately.
        for chrm in toPlotData['chrm'].unique():
            for periBarId in toPlotData['faceColor'].unique():
                tmpToPlotData = toPlotData[(toPlotData['chrm'] == chrm) & (toPlotData['faceColor'] == periBarId)]
                if tmpToPlotData.shape[0] == 0:
                    continue
                circle.barplot(chrm,
                               data=[1] * tmpToPlotData.shape[0],
                               positions=tmpToPlotData['position'].to_list(),
                               width=tmpToPlotData['width'].to_list(),
                               facecolor=tmpToPlotData['faceColor'].to_list(),
                               raxis_range=tmpToPlotData['periRange'].to_list()[0])

    def addLtrBar(self, circle, typList, bainfoList):
        # typList is a integer list which decides the periOrder and the color (to be controlled
        # consistent with 3dPlot)
        chrmList = []
        faceColorList = []
        widthList = []
        stPositionList = []
        periRangeList = []


        for typ, bain in zip(typList, bainfoList):
            chrmList.append(bain.chr)
            # faceColorList.append(self.colorList[typ-1])
            faceColorList.append(self.colorList[typ-1])

            periRangeList.append(self.circos_rRangeList[typ-1])
            widthList.append(self.ltrBarLen)
            stPositionList.append(bain.st)
        toPlotData = pd.concat([pd.Series(chrmList),
                                pd.Series(stPositionList),
                                pd.Series(widthList),
                                pd.Series(faceColorList),
                                pd.Series(periRangeList)],
                               axis=1)
        toPlotData.columns = ['chrm', 'position', 'width', 'faceColor', 'periRange']
        self.addPeriSegment(circle, toPlotData)

    def plot(self, typList, oriIdList, outFigFile, typId2annot):
        # typList should be at range [1, circos_maxPeriBarNum]. Very Important.
        # typList decides the periOrder (from outer circos to inner circos).
        # oriIdList is a oriIdList, further filtration would be done here,
        # to prevent unexpected oriId (those mei at scaffolds).
        # If chrDict is None, the default chrDict is the dict from the
        # ltrDistDetector of ltrDistManager.
        self.checkTypeList(lowestId=1,
                           highestId=self.circos_maxPeriBarNum,
                           typList=typList)
        typList, bainfoList, oriIdList = self.filterOriIdList(typList, oriIdList)
        circle = self.initFrame()
        self.setPeriTicks(circle)
        self.addLtrBar(circle, typList, bainfoList)
        self.addLegend(typList, typId2annot)
        self.addGenomeText(circle)
        circle.figure.savefig(outFigFile, bbox_inches='tight')
        outFigFile = outFigFile.split('.pdf')[0] + '.png'
        circle.figure.savefig(outFigFile, bbox_inches='tight', dpi=300)
        plt.close()


class p3DPlotter(plotterBase):
    def __init__(self, initSource, initTyp, colorList, ltrParaFile, speIdSet=None):
        super().__init__(initSource, initTyp, colorList, ltrParaFile, speIdSet)
        self.init_3DRela()

    def init_3DRela(self):
        self._3D_maxTypNum = len(self.colorList)

    def transTypList2plotColorList(self, typList):
        plotColorList = []
        for typ in typList:
            if typ == 0:
                plotColorList.append('#cdcdcd')
            else:
                plotColorList.append(self.colorList[typ - 1])
        return plotColorList

    def getLegendElementList(self, typId2annot):
        legendElementList = []
        for typ in typId2annot:
            if typ == 0:
                continue
            legendElementList.append(Line2D([], [],
                                            marker='.',
                                            color=self.colorList[typ - 1],
                                            label=typId2annot[typ],
                                            markersize=25,
                                            markeredgewidth=0,
                                            linestyle='None'))
        return legendElementList

    def plot(self, adata, typList, typId2annot, outGif, plotText=True, plotTree=True):
        # adata is the loaded adata of LTR_Stream.
        # typList should be at range [0, _3D_maxTypNum]
        # 0 for grey color.
        self.checkTypeList(lowestId=0,
                           highestId=self._3D_maxTypNum,
                           typList=typList)
        plotColorList = self.transTypList2plotColorList(typList)
        legendElementList = self.getLegendElementList(typId2annot)
        comp1, comp2, comp3 = 0, 1, 2
        df_plot = pd.DataFrame(index=adata.obs.index,
                               data=adata.obsm['X_dr'],
                               columns=['Dim' + str(x + 1) for x in range(adata.obsm['X_dr'].shape[1])])
        epg = adata.uns['epg']
        flat_tree = adata.uns['flat_tree']
        ft_node_pos = nx.get_node_attributes(flat_tree, 'pos')
        ft_node_label = nx.get_node_attributes(flat_tree, 'label')
        epg_node_pos = nx.get_node_attributes(epg, 'pos')
        fig = plt.figure(1, figsize=(13, 11))
        ax = plt.axes(projection='3d')

        def init():
            # Plot each frame
            ax.scatter3D(df_plot['Dim' + str(comp1 + 1)],
                         df_plot['Dim' + str(comp2 + 1)],
                         df_plot['Dim' + str(comp3 + 1)],
                         c=plotColorList,
                         linewidth=0.1)
            ax.set_xlabel('PCoA' + str(comp1 + 1), fontdict={'fontsize': 24})
            ax.set_ylabel('PCoA' + str(comp2 + 1), fontdict={'fontsize': 24})
            ax.set_zlabel('PCoA' + str(comp3 + 1), fontdict={'fontsize': 24})
            if plotTree:
                for edge_i in flat_tree.edges():
                    branch_i_pos = np.array([epg_node_pos[i] for i in flat_tree.edges[edge_i]['nodes']])
                    curve_i = pd.DataFrame(branch_i_pos, columns=np.arange(branch_i_pos.shape[1]))
                    ax.plot(curve_i[comp1], curve_i[comp2], curve_i[comp3], c='black', zorder=5)
            if plotText:
                for node_i in flat_tree.nodes():
                    # matplotlib has zcorder issues for 3D plot.To avoid the incorrect order,
                    # here the node point is removed for 3D
                    # ax_i.scatter(ft_node_pos[node_i][comp1],ft_node_pos[node_i][comp2],ft_node_pos[node_i][comp3],
                    #              color='#767070',s=1.5*(mpl.rcParams['lines.markersize']**2),zorder=10)
                    ax.text(ft_node_pos[node_i][comp1], ft_node_pos[node_i][comp2], ft_node_pos[node_i][comp3],
                            ft_node_label[node_i],
                            color='black', fontsize=mpl.rcParams['font.size'],
                            ha='left', va='bottom',
                            zorder=10)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_zticklabels([])
            ax.legend(handles=legendElementList, loc='upper center', ncol=3, fontsize=20)
            plt.grid()

        def rotate(angle):
            ax.view_init(elev=10, azim=angle)
        rot_animation = animation.FuncAnimation(fig, rotate,
                                                init_func=init,
                                                frames=np.arange(0, 362, 2),
                                                interval=50)
        rot_animation.save(outGif, dpi=170, writer='imagemagick')
        plt.close()


class meiPlotter:
    def __init__(self, initSource, initTyp, ltrParaFile):
        if not (initSource is None):

            self.colorList = ['#e6928f', '#6e9ece', '#BDB76B', '#4e9595', '#8d6ab8',
                            '#84574d', '#454545', ]
            self.maxTypNum = len(self.colorList)
            self.initSource = initSource
            self.initTyp = initTyp
            self.myCircosPlotter = circosPlotter(initSource=initSource,
                                                 initTyp=initTyp,
                                                 ltrParaFile=ltrParaFile,
                                                 colorList=self.colorList)
            self.myP3DPlotter = p3DPlotter(initSource=initSource,
                                           initTyp=initTyp,
                                           ltrParaFile=ltrParaFile,
                                           colorList=self.colorList)

    def init_3DRela(self):
        self._3D_maxTypNum = len(self.colorList)

    def p3DPlot(self, adata, typList, typId2annot, outGif, plotText=False, plotTree=True):
        # adata is the loaded adata of LTR_Stream.
        # typList should be at range [0, _3D_maxTypNum]
        # 0 for grey color.
        self.myP3DPlotter.plot(adata, typList, typId2annot, outGif, plotText=plotText, plotTree=plotTree)

    def circosPlot(self, typList, oriIdList, outFigFile, typId2annot):
        # typList should be at range [1, circos_maxPeriBarNum]. Very Important.
        # typList decides the periOrder (from outer circos to inner circos).
        # oriIdList is a oriIdList, further filtration would be done here,
        # to prevent unexpected oriId (those mei at scaffolds).
        # If chrDict is None, the default chrDict is the dict from the
        # ltrDistDetector of ltrDistManager.
        self.myCircosPlotter.plot(typList, oriIdList, outFigFile, typId2annot)
