from VCFRead import VCFRead
class sniffles2Read(VCFRead):
    def __init__(self, platform=None, depth=None, svType=None, vcfFile=None,outFile=None,
                 ref="/home/lz/Data_sequence/Ref/hs37d5.fa",mapTool=None, svTool=None,RE_nu=None, nux=None):
        VCFRead.__init__(self, platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,outFile=outFile,
                         ref=ref,mapTool=mapTool,svTool=svTool,RE_nu=RE_nu, nux=nux)

    def getRE(self,line):
        tmp = [i for i in str(line).split('\t')[7].split(';') if 'SUPPORT=' in i]
        return int(tmp[0].split('=')[1])
    
if __name__=="__main__":
    SvToolVcf = snifflesRead(platform='CCS', depth='25x',svType='BND', vcfFile='/home/lz/Data_sequence/2023_11_14/working/Pipelines/CCS/CHM13/lra/vcf/sniffles2/CHM13.vcf',outFile='CallOutPath/CCS/lra/sniffles2/BND/25x/2/CHM13.25x.vcf', ref='/home/lz/Data_sequence/2023_11_14/Genome/hg38/genome.fa', mapTool='lra', svTool='sniffles2', RE_nu=2, nux=0)

    SvToolVcf.run()
