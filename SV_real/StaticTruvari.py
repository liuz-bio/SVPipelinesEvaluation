import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

class StaticTruvari():
    def __init__(self):
        self.platform = ["Nanopore", "Pacbio"]
        self.depth = ["5", "10", "15", "25", "30"]
        self.svType = ["DEL", "INS", "INV", "DUP"]
        #self.mapTool = ["BWA-MEM", "lordfast", "minialign", "minimap2", "ngmlr", "pbmm2", "winnowmap"]
        self.mapTool = ["lordfast", "ngmlr", "pbmm2", "winnowmap", "minialign", "minimap2"]
        #self.callTool = ["cutesv", "debreak", "delly", "NanoSV", "nanovar", "pbsv", "picky", "sniffles", "svim"]
        self.callTool = ["cutesv", "debreak", "delly", "pbsv", "sniffles", "svim", "picky"]
        self.lines = []
        self.callToolLinesDict = {}
        self.mapToolLinesDict = {}
        self.svTypeLinesDict = {}
        self.depthLinesDict = {}
        self.platformLinesDict = {}
        self.callToolLines()
        self.mapToolLines()
        self.svTypeLines()
        self.depthLines()
        self.platformLines()
        self.ComGroup = [self.callToolLinesDict, self.mapToolLinesDict, self.svTypeLinesDict, self.depthLinesDict, self.platformLinesDict]

    def callToolLines(self):
        for platform in self.platform:
            for depth in self.depth:
                for svType in self.svType:
                    for mapTool in self.mapTool:
                        id = "_".join([platform, mapTool, "callTool", svType, depth])
                        summeryFiles = {}
                        for callTool in self.callTool:
                            summeryFiles[callTool] = os.path.join("truvari", platform, mapTool, callTool, svType, depth, "_".join([platform.lower(), svType.lower() ,depth+'x']), "summary.txt")
                        self.callToolLinesDict[id] = summeryFiles


    def mapToolLines(self):
        for platform in self.platform:
            for depth in self.depth:
                for svType in self.svType:
                    for callTool in self.callTool:
                        id = "_".join([platform, "mapTool", callTool, svType, depth])
                        summeryFiles = {}
                        for mapTool in self.mapTool:
                            summeryFiles[mapTool] = os.path.join("truvari", platform, mapTool, callTool, svType, depth, "_".join([platform.lower(), svType.lower() ,depth+'x']), "summary.txt")
                        self.mapToolLinesDict[id] = summeryFiles

    def svTypeLines(self):
        for platform in self.platform:
            for depth in self.depth:
                for callTool in self.callTool:
                    for mapTool in self.mapTool:
                        id = "_".join([platform, mapTool, callTool, "svType", depth])
                        summeryFiles = {}
                        for svType in self.svType:
                            summeryFiles[svType] = os.path.join("truvari", platform, mapTool, callTool, svType, depth, "_".join([platform.lower(), svType.lower() ,depth+'x']), "summary.txt")
                        self.svTypeLinesDict[id] = summeryFiles

    def depthLines(self):
        for platform in self.platform:
            for callTool in self.callTool:
                for svType in self.svType:
                    for mapTool in self.mapTool:
                        id = "_".join([platform, mapTool, callTool, svType, "depth"])
                        summeryFiles = {}
                        for depth in self.depth:
                            summeryFiles[depth] = os.path.join("truvari", platform, mapTool, callTool, svType, depth, "_".join([platform.lower(), svType.lower() ,depth+'x']), "summary.txt")
                        self.depthLinesDict[id] = summeryFiles

    def platformLines(self):
        for callTool in self.callTool:
            for depth in self.depth:
                for svType in self.svType:
                    for mapTool in self.mapTool:
                        id = "_".join(["platform", mapTool, callTool, svType, depth])
                        summeryFiles = {}
                        for platform in self.platform:
                            summeryFiles[platform] = os.path.join("truvari", platform, mapTool, callTool, svType, depth, "_".join([platform.lower(), svType.lower() ,depth+'x']), "summary.txt")
                        self.platformLinesDict[id] = summeryFiles

    def summaryRead(self,summary):
        oneSumm = eval(open(summary,'r').read())
        line = []
        for i in ['TP-base', 'TP-call', 'FP', 'FN', 'precision', 'recall', 'f1', 'base cnt', 'call cnt', 'TP-call_TP-gt', 'TP-call_FP-gt', 'TP-base_TP-gt', 'TP-base_FP-gt', 'gt_precision', 'gt_recall', 'gt_f1']:
            line.append(oneSumm[i])
        return line

    def plotbar(self,df, figName):
        df = df.fillna(0, inplace=False)
        if "_delly_" in figName and "_svType_" in figName:
            try:
                df =df.loc[["DEL","INS"],:]
            except KeyError:
                print(figName)
                print(df)
        labels = [i for i in df.index.tolist()]
        groups = [df[i].tolist() for i in df.columns.tolist()]
        groupsLabe = [i for i in df.columns.tolist()]
        x = np.arange(len(labels))  # the label locations

        width = 0.9/ (len(labels) + len(groups))  # the width of the bars
        span = 0.3 / (len(labels) + len(groups) - 1)

        fig, ax = plt.subplots(figsize=(7, 4.326))
        axx = len(groups) * (width) + (len(groups) - 1) * span / 2
        a = x + width / 2 - (len(groups) * (width) + (len(groups) - 1) * span / 2)
        xPos = [axx / 2 - span / 2 + a + (width + span) * i for i in range(0, len(groups) + 0)]

        for i in range(len(groups)):
            if groupsLabe[i] == "TP-call":
                ax2 = ax.twinx()
                rects2 = ax2.bar(xPos[i], groups[i], width, color="red", label=groupsLabe[i])
                ax2.bar_label(rects2, padding=3)
            else:
                rects1 = ax.bar(xPos[i], groups[i], width, label=groupsLabe[i])
                ax.bar_label(rects1, fmt="%.2f", padding=3)
        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_ylabel('Percentage')
        ax.set_title(figName.split('/')[-1].split('.')[0])
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylim(0, 1.15)
        ax.legend(loc=8, bbox_to_anchor=(0.6, -0.17), ncol=4, frameon=False)
        ax2.legend(loc=8, bbox_to_anchor=(0.2, -0.17), ncol=4, frameon=False)
        ax2.set_ylabel("Number of SV")
        ax2.set_ylim(0, max(df.loc[:,'TP-call'].tolist())+500)


        #fig.tight_layout()
        plt.savefig(figName, format='png')
        plt.close()

    def run(self):
        for groups in self.ComGroup:
            indexs = []
            for groupK, groupV in groups.items():
                groupSumm = []
                saveCsv = True
                for k, summary in groupV.items():
                    if os.path.exists(summary):
                        groupSumm.append(self.summaryRead(summary))
                        #print(summary)
                    else:
                        #if ("cutesv" in summary) and ("ngmlr" in summary) and ("25" in summary):
                        print(summary)
                        saveCsv = False
                        break;
                    if k not  in indexs:
                        indexs.append(k)
                #print(groupSumm)
                #print(groupV)
                if saveCsv:
                    cols = ['TP-base', 'TP-call', 'FP', 'FN', 'precision', 'recall', 'f1', 'base-cnt', 'call-cnt', 'TP-call_TP-gt', 'TP-call_FP-gt', 'TP-base_TP-gt', 'TP-base_FP-gt', 'gt_precision', 'gt_recall', 'gt_f1']
                    try:
                        df = pd.DataFrame(groupSumm, columns=cols, index=indexs, dtype='float32')
                    except:
                        print(groupSumm)
                        print(len(groupSumm))
                        print(len(cols))
                        print(indexs)
                    os.makedirs("summary", exist_ok=True)
                    df.to_csv("summary/%s.csv"%groupK, encoding="utf-8", header=True, index=True)
                    figName = "summary/%s.png"%groupK
                    self.plotbar(df[['TP-call','precision', 'recall', 'f1']], figName)



