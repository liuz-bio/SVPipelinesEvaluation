from VCFRead import VCFRead


class pbsvRead(VCFRead):
    def __init__(self, platform=None, depth=None, svType=None, vcfFile=None,outFile=None,
                 ref="/home/lz/Data_sequence/Ref/hs37d5.fa",mapTool=None, svTool=None,RE_nu=None, nux=None):
        VCFRead.__init__(self, platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,outFile=outFile,
                 ref=ref,mapTool=mapTool,svTool=svTool,RE_nu=RE_nu, nux=nux)

    def getSample(self, line):
        GT = '/'.join([str(ix) for ix in line.samples.values()[0]['GT']])
        #AD = line.samples.values()[0]['AD']
        DR = '5' #str(AD[0])
        DV = '5' #str(AD[1])
        sample = ':'.join([GT, DR, DV])
        return sample, GT

    def getRE(self,line):
        AD = [int(i) for i in line.samples.values()[0]['AD']]
        GT = '/'.join([str(ix) for ix in line.samples.values()[0]['GT']])
        if GT == '1/1':
            return max(AD)
        elif GT == '0/1':
            return abs(AD[0]-AD[1])
        return 1

    def writetmp(self):
        with open(self.tmpFileName,'w') as out:
            for hl in open("Head.pbsv.info",'r'):
                out.write(hl)
            with open(self.vcfFile, 'r') as inp:
                for li in inp:
                    if "#" in li:
                        continue
                    elif '=>' in li:
                        li = li.replace("=>",'=')
                        out.write(li)
                    else:
                        out.write(li)
