from VCFRead import VCFRead
class svimRead(VCFRead):
    def __init__(self, platform=None, depth=None, svType=None, vcfFile=None,outFile=None,
                 ref="/home/lz/Data_sequence/Ref/hs37d5.fa",mapTool=None, svTool=None,RE_nu=None, nux=None):
        VCFRead.__init__(self, platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,outFile=outFile,
                 ref=ref,mapTool=mapTool,svTool=svTool,RE_nu=RE_nu, nux=nux)

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


