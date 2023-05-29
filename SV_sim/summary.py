import os
import sys

class summary():
    def __int__(self,vcfFile=None):
        self.vcfFile = vcfFile
        self.F1_GT = [] 
        self.F1_no_GT = []
        self.precision_GT = []
        self.precision_no_GT = []
        self.recall_GT = []
        self.recall_no_GT = []
        self.lsv_base_GT = []
        self.lsv_comm_GT = []
        self.lsv_base_no_GT = []
        self.lsv_comm_no_GT = []
        self.bc1_GT = []
        self.bc2_GT = []
        self.lsv_bias_GT = []
        self.bc1_no_GT = []
        self.bc2_no_GT = []
        self.lsv_bias_no_GT = []
        self.bc1_bc2_no_GT = []
        self.bc1_bc2_GT = []

    def outWrite(self,outFile,outL):
        with open(outFile,'w') as out:
            for line in outL:
                out.write(line+'\n')

    def summarySV(InfoName):
        for line in open(InfoName,'r'):
            lines = line.strip().split('\t')
            vcfFileName = lines[-1]
            depth = lines[1]
            platform = lines[0]
            svType = lines[2].upper()
            if svType == "TRA":
                svType = "BND"
            mapTool = lines[3]
            svTool = lines[4]
            EvalOutFile = os.path.join("EvalOutFile",platform, mapTool, svTool, svType, depth, '.'.join(vcfFileName.split('.')[:-1]))

-rw-r--r-- 1 lz lab1 8.6K Oct  1 20:57 bias.txt
-rw-r--r-- 1 lz lab1 337K Oct  1 20:57 FN_GT.vcf
-rw-r--r-- 1 lz lab1 322K Oct  1 20:57 FN_no_GT.vcf
-rw-r--r-- 1 lz lab1  57K Oct  1 20:57 FP_GT.vcf
-rw-r--r-- 1 lz lab1  29K Oct  1 20:57 FP_no_GT.vcf
-rw-r--r-- 1 lz lab1  193 Oct  1 20:57 lsv.txt
-rw-r--r-- 1 lz lab1  186 Oct  1 20:57 summary.txt
-rw-r--r-- 1 lz lab1  30K Oct  1 20:57 TP_GT.vcf
-rw-r--r-- 1 lz lab1  59K Oct  1 20:57 TP_no_GT.vcf 
if __name__=='__main__':
   sumy = summary('aa.tab')
