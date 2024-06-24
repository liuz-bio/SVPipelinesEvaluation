from VCFRead import VCFRead
class snifflesRead(VCFRead):
    def __init__(self, platform=None, depth=None, svType=None, vcfFile=None,outFile=None,
                 ref="/home/lz/Data_sequence/Ref/hs37d5.fa",mapTool=None, svTool=None,RE_nu=None, nux=None):
        VCFRead.__init__(self, platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,outFile=outFile,
                         ref=ref,mapTool=mapTool,svTool=svTool,RE_nu=RE_nu, nux=nux)

#tmp.9.6144.CCS.25x.minimap2.sniffles.vcf
if __name__=="__main__":
    SvToolVcf = snifflesRead(platform='CCS', depth='25x',svType='BND', vcfFile='tmp.9.6144.CCS.25x.minimap2.sniffles.vcf',outFile='CallOutPath/CCS/minimap2/sniffles/BND/25x/2/CHM13.25x.vcf', ref='/home/lz/Data_sequence/2023_11_14/Genome/hg38/genome.fa', mapTool='minimap2', svTool='sniffles', RE_nu=2, nux=0)
#    SvToolVcf = snifflesRead(platform='CCS', depth='25x',svType='BND', vcfFile='/home/lz/Data_sequence/2023_11_14/working/Pipelines/CCS/CHM13/minimap2/vcf/sniffles/CHM13.vcf',outFile='CallOutPath/CCS/minimap2/sniffles/BND/25x/2/CHM13.25x.vcf', ref='/home/lz/Data_sequence/2023_11_14/Genome/hg38/genome.fa', mapTool='minimap2', svTool='sniffles', RE_nu=2, nux=0)

    SvToolVcf.run()
