import os
import sys


class overlapVcf():
    def __init__(self):
        self.platform = ["Nanopore", "Pacbio"]
        self.depth = ["5", "10", "15", "25", "30"]
        self.svType = ["DEL", "INS", "INV", "DUP"]
        #self.mapTool = ["BWA-MEM", "lordfast", "minialign", "minimap2", "ngmlr", "pbmm2", "winnowmap"]
        self.mapTool = ["lordfast", "ngmlr", "pbmm2", "winnowmap", "minialign", "minimap2"]
        #self.callTool = ["cutesv", "debreak", "delly", "NanoSV", "nanovar", "pbsv", "picky", "sniffles", "svim"]
        self.callTool = ["cutesv", "debreak", "delly", "pbsv", "sniffles", "svim"]
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
                            summeryFiles[callTool] = os.path.join("CallOutPath", platform, mapTool, callTool, svType, depth, "_".join([platform.lower(), svType.lower() ,depth+'x.vcf']))
                        self.callToolLinesDict[id] = summeryFiles


    def mapToolLines(self):
        for platform in self.platform:
            for depth in self.depth:
                for svType in self.svType:
                    for callTool in self.callTool:
                        id = "_".join([platform, "mapTool", callTool, svType, depth])
                        summeryFiles = {}
                        for mapTool in self.mapTool:
                            summeryFiles[mapTool] = os.path.join("truvari", platform, mapTool, callTool, svType, depth, "_".join([platform.lower(), svType.lower() ,depth+'x.vcf']))
                        self.mapToolLinesDict[id] = summeryFiles

    def svTypeLines(self):
        for platform in self.platform:
            for depth in self.depth:
                for callTool in self.callTool:
                    for mapTool in self.mapTool:
                        id = "_".join([platform, mapTool, callTool, "svType", depth])
                        summeryFiles = {}
                        for svType in self.svType:
                            summeryFiles[svType] = os.path.join("truvari", platform, mapTool, callTool, svType, depth, "_".join([platform.lower(), svType.lower() ,depth+'x.vcf']))
                        self.svTypeLinesDict[id] = summeryFiles

    def depthLines(self):
        for platform in self.platform:
            for callTool in self.callTool:
                for svType in self.svType:
                    for mapTool in self.mapTool:
                        id = "_".join([platform, mapTool, callTool, svType, "depth"])
                        summeryFiles = {}
                        for depth in self.depth:
                            summeryFiles[depth] = os.path.join("truvari", platform, mapTool, callTool, svType, depth, "_".join([platform.lower(), svType.lower() ,depth+'x.vcf']))
                        self.depthLinesDict[id] = summeryFiles

    def platformLines(self):
        for callTool in self.callTool:
            for depth in self.depth:
                for svType in self.svType:
                    for mapTool in self.mapTool:
                        id = "_".join(["platform", mapTool, callTool, svType, depth])
                        summeryFiles = {}
                        for platform in self.platform:
                            summeryFiles[platform] = os.path.join("truvari", platform, mapTool, callTool, svType, depth, "_".join([platform.lower(), svType.lower() ,depth+'x.vcf']))
                        self.platformLinesDict[id] = summeryFiles

    def run(self):
        for groups in self.ComGroup:
            indexs = []
            for groupK, groupV in groups.items():
                groupSumm = []
                saveCsv = True
                with open('tmp.info','w') as f:
                    f.write('\n'.join(groupV.values())+'\n')

                os.system("SURVIVOR merge tmp.info 300 2 1 1 0 30 tmp.info.vcf")

                if os.path.exists("tmp.info.vcf"):
                    groupName = list(groupV.keys())
                    groupDict = dict(zip(groupName, [[] for i in range(len(groupName))]))
                    for line in open("tmp.info.vcf",'r'):
                        if "#" in line:
                            continue
                        else:
                            line = line.strip().split('\t')
                            svID = line[2]
                            samples = [ix.split(':')[7].split('.')[0].split('-') for ix in line[-len(groupName):]]
                            print(samples)
                            print(svID)
                            print(groupName)
                            for sample in groupName:
                                for sampleList in samples:
                                    if sample in sampleList:
                                        groupDict[sample].append(svID)
                    os.makedirs("overlap", exist_ok=True)
                    for sample in groupName:
                        os.makedirs("overlap/%s"%groupK, exist_ok=True)
                        with open("overlap/%s/%s.txt"%(groupK, sample), 'w') as f:
                            f.write('\n'.join(groupDict[sample])+'\n')
                    print("overlap/%s/%s.txt"%(groupK, sample))
                    sys.exit()
                else:
                    print('############################################')
                    continue
