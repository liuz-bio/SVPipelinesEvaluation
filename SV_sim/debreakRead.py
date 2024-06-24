from VCFRead import VCFRead
import os


class debreakRead(VCFRead):
    def __init__(self, platform=None, depth=None, svType=None, vcfFile=None,outFile=None,
                 ref="/home/lz/Data_sequence/Ref/hs37d5.fa",mapTool=None, svTool=None,RE_nu=None, nux=None):
        VCFRead.__init__(self, platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,outFile=outFile,
                 ref=ref,mapTool=mapTool,svTool=svTool,RE_nu=RE_nu, nux=nux)

    def getSample(self,line):
        GT = '/'.join([str(ix) for ix in line.samples.values()[0]['GT']])
        DR = '5' #str(int(int(self.depth)/2))
        DV = '5' #str(int(int(self.depth)/2))
        sample = ':'.join([GT, DR, DV])
        return sample, GT

    def newHead(self):
        self.vcfHead.insert(-1,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
        config = {"1": "249250621", "2": "243199373", "3": "198022430", \
         "4": "191154276", "5": "180915260", "6": "171115067", \
         "7": "159138663", "8": "146364022", "9": "141213431", \
         "10": "135534747", "11": "135006516", "12": "133851895", \
         "13": "115169878", "14": "107349540", "15": "102531392", \
         "16": "90354753", "17": "81195210", "18": "78077248", \
         "19": "59128983", "20": "63025520", "21": "48129895", \
         "22": "51304566", "X": "155270560", "Y": "59373566", \
         "MT": "16569"}
        nu=4
        for k,v in config.items():
            self.vcfHead.insert(nu,"##contig=<ID=%s,length=%s>"%(k,v))
            nu+=1

    def OutPut(self):
        os.makedirs('/'.join(self.outFile.split('/')[:-1]), exist_ok=True)
        out = open(self.outFile, "w")
        #self.newHead()
        for line in self.vcfHead + self.SvLines:
            out.write(line + '\n')
        out.close()
        #os.system("cat %s|vcf-sort|bgzip -c >%s" % (self.outFile, self.outFile + '.gz'))
        #os.system("tabix %s" % self.outFile + '.gz')

    def getRE(self,line):
        tmp = [i for i in str(line).split('\t')[7].split(';') if 'SUPPREAD=' in i]
        return int(tmp[0].split('=')[1])

    def getCHR2Pos(self,line):
        tmp = [ i for i in str(line).split('\t')[7].split(';') if 'CHR2=' in i]
        CHR2 = tmp[0].split('=')[1]
        return [CHR2, str(line.stop)]
