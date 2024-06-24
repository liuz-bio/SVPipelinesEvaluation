from VCFRead import VCFRead
class pickyRead(VCFRead):
    def __init__(self, platform=None, depth=None, svType=None, vcfFile=None,outFile=None,
                 ref="/home/lz/Data_sequence/Ref/hs37d5.fa",mapTool=None, svTool=None,RE_nu=None, nux=None):
        VCFRead.__init__(self, platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,outFile=outFile,
                 ref=ref,mapTool=mapTool,svTool=svTool,RE_nu=RE_nu, nux=nux)

    def getSample(self,line):
        #GT = '/'.join([str(ix) for ix in line.samples.values()[0]['GT']])
        #AD = line.samples.values()[0]['AD']
        #DR = str(AD[0])
        #DV = str(AD[1])
        GT = '0/1'
        DR = '12'
        DV = '12'
        sample = ':'.join([GT, DR, DV])
        return sample, GT

if __name__=="__main__":
    SvToolVcf = pickyRead(platform='CCS', depth='25x',svType='DEL', vcfFile='tmp.9.3940.CCS.25x.winnowmap.picky.vcf',outFile='w/CHM13.25x.vcf', ref='/home/lz/Data_sequence/2023_11_14/Genome/hg38/genome.fa', mapTool='winnowmap', svTool='picky', RE_nu=2, nux=0)

    SvToolVcf.run()
