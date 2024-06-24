import pysam
import os


class createSimulateSV():
    def __init__(self,SimulateBedFile=None,svType=None,outFile=None,callSVFile=None,ref=None,depth=None):
        self.SimulateBedFile = SimulateBedFile
        self.svType = svType
        self.outFile = outFile
        self.depth = depth
        print('self.svType',self.svType)
        print('SimulateBedFile',SimulateBedFile)
        #self.callSVFile = callSVFile
        #self.vcfHead = str(pysam.VariantFile(callSVFile).header).strip().split('\n')
        self.vcfHead = open("Head.info", "r").read().strip().split("\n")
        self.svLines = []
        self.genome = pysam.FastaFile(ref)
        self.nus = 1

    def creatDEL(self,line):
        svChr = str(line[0])
        svPos = str(line[1])
        svID = "Sim."+str(self.nus)
        self.nus +=1
        svEnd = 'END=' + str(line[2])
        svType = "SVTYPE=" + str(self.svType)
        aSeq = self.genome.fetch(svChr, int(svPos)+1, int(line[2]) + 1)
        bSeq = self.genome.fetch(svChr, int(svPos), int(svPos) + 1)
        svQual = '.'
        svLEN = 'SVLEN=-' + str(len(aSeq))
        svFilter = 'PASS'
        #GT = "0/1"
        GT = line[-1].split(':')[0]
        DR = str(int(int(self.depth)/2))
        DV = str(int(int(self.depth)/2))
        svinfo = ';'.join(['PRECISE', svEnd, svType,svLEN])
        sample = ':'.join([GT, DR, DV])
        self.svLines.append('\t'.join([svChr,svPos,svID,aSeq,bSeq,svQual,svFilter,svinfo,"GT:DR:DV",sample]))

    def creatINS(self,line):
        svChr = str(line[0])
        svPos = str(line[1])
        svID = "Sim."+str(self.nus)
        self.nus += 1
        svEnd = 'END=' + str(line[2])
        svType = "SVTYPE=" + str(self.svType)
        aSeq = self.genome.fetch(svChr, int(svPos), int(svPos) + 1)
        bSeq = str(line[4]).upper()
        svLEN = "SVLEN=" + str(len(bSeq))
        svQual = '.'
        svFilter = 'PASS'
        #GT = "0/1"
        GT = line[-1].split(':')[0]
        DR = str(int(int(self.depth) / 2))
        DV = str(int(int(self.depth) / 2))
        svinfo = ';'.join(['PRECISE', svEnd, svType, svLEN])
        sample = ':'.join([GT, DR, DV])
        self.svLines.append('\t'.join([svChr, svPos, svID, aSeq, bSeq, svQual, svFilter, svinfo, "GT:DR:DV", sample]))

    def creatINV(self,line):
        svChr = str(line[0])
        svPos = str(line[1])
        svID = "Sim."+str(self.nus)
        self.nus += 1
        svEnd = 'END=' + str(line[2])
        svLEN = 'SVLEN=' + str(int(line[2]) - int(line[1]))
        svType = "SVTYPE=" + str(self.svType)
        aSeq = self.genome.fetch(svChr, int(svPos), int(svPos) + 1)
        bSeq = "<INV>"
        svQual = '.'
        svFilter = 'PASS'
        #GT = "0/1"
        GT = line[-1].split(':')[0]
        DR = str(int(int(self.depth) / 2))
        DV = str(int(int(self.depth) / 2))
        svinfo = ';'.join(['PRECISE', svEnd, svType, svLEN])
        sample = ':'.join([GT, DR, DV])
        self.svLines.append('\t'.join([svChr, svPos, svID, aSeq, bSeq, svQual, svFilter, svinfo, "GT:DR:DV", sample]))

    def creatDUP(self,line):
        svChr = str(line[0])
        svPos = str(line[1])
        svID = "Sim."+str(self.nus)
        self.nus += 1
        svEnd = 'END=' + str(line[2])
        svLEN = 'SVLEN=' + str(int(line[2]) - int(line[1]))
        svType = "SVTYPE=" + str(self.svType)
        aSeq = self.genome.fetch(svChr, int(svPos), int(svPos) + 1)
        bSeq = "<DUP>"
        svQual = '.'
        svFilter = 'PASS'
        #GT = "0/1"
        GT = line[-1].split(':')[0]
        DR = str(int(int(self.depth) / 2))
        DV = str(int(int(self.depth) / 2))
        svinfo = ';'.join(['PRECISE', svEnd, svType, svLEN])
        sample = ':'.join([GT, DR, DV])
        self.svLines.append('\t'.join([svChr, svPos, svID, aSeq, bSeq, svQual, svFilter, svinfo, "GT:DR:DV", sample]))

    def creatTRA(self,line):
       svChr = str(line[0])
       svPos = str(line[1])
       aSeq = self.genome.fetch(svChr, int(svPos), int(svPos) + 1)
       bSeq = "<TRA>"
       svID = "Sim."+str(self.nus)
       self.nus += 1
       svEnd = 'END=' + str(line[4].split(':')[2])
       svLEN = 'SVLEN=NULL'
       svFilter = 'PASS'
       svQual = '.'
       GT = "0/1"
       DR = str(int(int(self.depth) / 2))
       DV = str(int(int(self.depth) / 2))
       svType = "SVTYPE=" + "TRA;CHR2="+line[4].split(':')[1]
       sample = ':'.join([GT, DR, DV])
       svinfo = ';'.join(['PRECISE', svEnd, svType, svLEN])
       self.svLines.append('\t'.join([svChr, svPos, svID, aSeq, bSeq, svQual, svFilter, svinfo, "GT:DR:DV", sample]))
    
    def OutPut(self):
        os.makedirs('/'.join(self.outFile.split('/')[:-1]), exist_ok=True)
        out = open(self.outFile,"w")
        for line in self.vcfHead+self.svLines:
            out.write(line+'\n')
        out.close()
        os.system("cat %s|vcf-sort|bgzip -c >%s"%(self.outFile,self.outFile+'.gz'))
        os.system("tabix %s"%self.outFile+'.gz')

    def run(self):
        for line in open(self.SimulateBedFile,'r'):
            lines = line.strip().split('\t')
            if self.svType == "DEL":
                self.creatDEL(lines)
            elif self.svType == "INS":
                self.creatINS(lines)
            elif self.svType == "DUP":
                self.creatDUP(lines)
            elif self.svType == "INV":
                self.creatINV(lines)
            elif self.svType == "BND":
                self.creatTRA(lines)

        self.OutPut()

