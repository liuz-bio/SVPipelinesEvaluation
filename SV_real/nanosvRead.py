from VCFRead import VCFRead
class nanosvRead(VCFRead):
    def __init__(self, platform=None, depth=None, svType=None, vcfFile=None,outFile=None,
                 ref="/home/lz/Data_sequence/Ref/hs37d5.fa",mapTool=None, svTool=None, RE_nu=None, nux=None):
        VCFRead.__init__(self, platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,outFile=outFile,
                 ref=ref,mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)

    def getSample(self,line):
        try:
            GT = '/'.join([str(ix) for ix in line.samples.values()[0]['GT']])
        except:
            return ':'.join(['0/0', '5,' '5']), False
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

#tmp.9.85412.ONT.25x.ngmlr.nanosv.vcf
if __name__=="__main__":
    #SvToolVcf = nanosvRead(platform='ONT', depth='25x',svType='BND', vcfFile='/home/lz/Data_sequence/2023_11_14/working/Pipelines/ONT/CHM13/ngmlr/vcf/NanoSV/CHM13.vcf',outFile='CallOutPath/ONT/ngmlr/nanosv/BND/25x/2/CHM13.25x.vcf', ref='/home/lz/Data_sequence/2023_11_14/Genome/hg38/genome.fa', mapTool='ngmlr', svTool='nanosv', RE_nu=2, nux=0)
    SvToolVcf = nanosvRead(platform='ONT', depth='25x',svType='BND', vcfFile='tmp.13.104891.ONT.25x.pbmm2.nanosv.vcf',outFile='CallOutPath/ONT/ngmlr/nanosv/BND/25x/2/CHM13.25x.vcf', ref='/home/lz/Data_sequence/2023_11_14/Genome/hg38/genome.fa', mapTool='ngmlr', svTool='nanosv', RE_nu=2, nux=0)
    SvToolVcf.run()

