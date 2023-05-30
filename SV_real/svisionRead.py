from VCFRead import VCFRead
class svisionRead(VCFRead):
    def __init__(self, platform=None, depth=None, svType=None, vcfFile=None, outFile=None,
                 ref="/home/lz/Data_sequence/Ref/hs37d5.fa",mapTool=None, svTool=None,RE_nu=None, nux=None):
        VCFRead.__init__(self, platform=platform, depth=depth, svType=svType, vcfFile=vcfFile, outFile=outFile,
                         ref=ref,mapTool=mapTool,svTool=svTool,RE_nu=RE_nu, nux=nux)

    def getSample(self,line):
        GT = '/'.join([str(ix) for ix in line.samples.values()[0]['GT']])
        DR = '5' #str(line.samples.values()[0]['DR'])
        DV = '5' #str(line.samples.values()[0]['DV'])
        sample = ':'.join([GT, DR, DV])
        return sample, GT

    def getRE(self,line):
        tmp = [i for i in str(line).split('\t')[7].split(';') if 'SUPPORT=' in i]
        return int(tmp[0].split('=')[1])
   
    def get_passinfo(self,line):
        if 'Covered' in str(line.filter.keys()[0]):
            passinfo = 'PASS'
        else:
            passinfo = 'Covered'

        return passinfo

    def getSVEND(self,line):
        if self.svType == 'INS':
            return 'END=' + str(line.pos+1)
        else:
            return 'END=' + str(line.stop)

