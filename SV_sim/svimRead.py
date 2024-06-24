from VCFRead import VCFRead
class svimRead(VCFRead):
    def __init__(self, platform=None, depth=None, svType=None, vcfFile=None,outFile=None,
                 ref="/home/lz/Data_sequence/Ref/hs37d5.fa",mapTool=None, svTool=None,RE_nu=None, nux=None):
        VCFRead.__init__(self, platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,outFile=outFile,
                 ref=ref,mapTool=mapTool,svTool=svTool,RE_nu=RE_nu, nux=nux)

    def chSVTYPE(self,line):
        if line.info['SVTYPE']=='TRA':
            return 'BND'
        elif line.info['SVTYPE'] == "DUP:TANDEM" and self.svTool.lower() == "svim":
            print("-----------------------------------",line.info['SVTYPE'],self.svTool.lower(),self.svType)
            return "DUP"
        else:
            return str(line.info['SVTYPE'])

    def getInfo(self,line):
        svTYPE = 'SVTYPE=' + self.svType
        svEND = 'END=' + str(line.stop)
        try:
            try:
                svLEN = 'SVLEN=' + str(line.info['SVLEN'][0])
                if str(line.info['SVTYPE']) in ["BND", "TRA"]:
                    svLEN = 'SVLEN=' + "NULL"
                    CHR2Pos = self.getCHR2Pos(line)
                    svTYPE = "SVTYPE=" + "BND;CHR2="+CHR2Pos[0]
                    svEND = 'END=' + CHR2Pos[1]
            except TypeError:
                svL = str(line.info['SVLEN'])
                svLEN = 'SVLEN=' + svL
                if str(line.info['SVTYPE']) in ["BND", "TRA"]:
                    svLEN = 'SVLEN=' + "NULL"
                    CHR2Pos = self.getCHR2Pos(line)
                    svTYPE = "SVTYPE=" + "BND;CHR2="+CHR2Pos[0]
                    svEND = 'END=' + CHR2Pos[1]
            infos = ';'.join(['PRECISE', svEND, svTYPE, svLEN])
        except KeyError:
            if str(line.info['SVTYPE']) == "INS":
                 svLEN = 'SVLEN=' + str(len(line.alts[0]))
            elif self.svType in ["DEL","INV","DUP"]:
                svLEN = 'SVLEN=' + str(int(self.getSVEND(line).split('=')[-1]) - int(self.getPos(line)))
            elif self.svType in ["BND", "TRA"]:
                svLEN = 'SVLEN=' + "NULL"
                CHR2Pos = self.getCHR2Pos(line)
                svTYPE = "SVTYPE=" + "BND;CHR2="+CHR2Pos[0]
                svEND = 'END=' + CHR2Pos[1]
            infos = ';'.join(['PRECISE', svEND, svTYPE, svLEN])
        return infos

    def getSample(self,line):
        GT = '/'.join([str(ix) for ix in line.samples.values()[0]['GT']])
        #AD = line.samples.values()[0]['AD']
        DR = '5' #str(AD[0])
        DV = '5' #str(AD[1])
        if DR == 'None' or DV == 'None':
            DR = '5'
            DV = '5'
        sample = ':'.join([GT, DR, DV])
        return sample, GT

    def getRE(self,line):
        tmp = [i for i in str(line).split('\t')[7].split(';') if 'SUPPORT=' in i]
        return int(tmp[0].split('=')[1])

