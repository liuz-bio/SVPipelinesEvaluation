from VCFRead import VCFRead
class nanosvRead(VCFRead):
    def __init__(self, platform=None, depth=None, svType=None, vcfFile=None,outFile=None,
                 ref="/home/lz/Data_sequence/Ref/hs37d5.fa",mapTool=None, svTool=None, RE_nu=None, nux=None):
        VCFRead.__init__(self, platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,outFile=outFile,
                 ref=ref,mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)

    def getSample(self,line):
        GT = '/'.join([str(ix) for ix in line.samples.values()[0]['GT']])
        DV = line.samples.values()[0]['DV']
        DR = line.samples.values()[0]['DR']
        DR = str(DR[0])
        DV = str(DV[0])
        sample = ':'.join([GT, DR, DV])
        return sample, GT

    def getRE(self,line):
        RE = line.samples.values()[0]['DV'][0]
        return int(RE)


    def get_passinfo(self,line):
        if 'LowQual' not in str(line.filter.keys()[0]):
            passinfo = "PASS"
        else:
            passinfo = "LowQual"

        return passinfo
