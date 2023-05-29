from VCFRead import VCFRead
class dellyRead(VCFRead):
    def __init__(self, platform=None, depth=None, svType=None, vcfFile=None,outFile=None,
                 ref="/home/lz/Data_sequence/Ref/hs37d5.fa",mapTool=None, svTool=None,RE_nu=None, nux=None):
        VCFRead.__init__(self, platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,outFile=outFile,
                 ref=ref,mapTool=mapTool,svTool=svTool,RE_nu=RE_nu, nux=nux)

    def getSample(self,line):
        GT = '/'.join([str(ix) for ix in line.samples.values()[0]['GT']])
        DR = '5' #str(line.samples.values()[0]['RR'])
        DV = '5' #str(line.samples.values()[0]['RV'])
        sample = ':'.join([GT, DR, DV])
        return sample,GT

    def getRE(self,line):
        tmp = [i for i in str(line).split('\t')[7].split(';') if 'SR=' in i]
        return int(tmp[0].split('=')[1])
