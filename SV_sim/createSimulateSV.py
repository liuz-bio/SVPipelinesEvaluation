import pysam
import os


class createSimulateSV():
    def __init__(self,SimulateBedFile=None,svType=None,outFile=None,callSVFile=None,ref=None,depth=None):
        self.SimulateBedFile = SimulateBedFile
        self.svType = svType
        self.outFile = outFile
        self.depth = depth
        #self.callSVFile = callSVFile
        #self.vcfHead = str(pysam.VariantFile(callSVFile).header).strip().split('\n')
        self.vcfHead = open("Head.info", "r").read().strip().split("\n")
        self.svLines = []
        self.genome = pysam.FastaFile(ref)
        self.nus = 1
        self.Genotype = {'1':'0/1','2':'1/1','3':'0/1','4':'0/1','5':'1/1','6':'0/1','7':'1/1','8':'0/1','9':'1/1','10':'0/1','11':'1/1','12':'0/1','13':'1/1','14':'0/1','15':'1/1','16':'0/1','17':'1/1','18':'0/1','19':'1/1','20':'0/1','21':'1/1','22':'0/1','X':'1/1','Y':'1/1','MT':'1/1'}

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
        GT = self.Genotype[svChr]
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
        GT = self.Genotype[svChr]
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
        GT = self.Genotype[svChr]
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
        GT = self.Genotype[svChr]
        DR = str(int(int(self.depth) / 2))
        DV = str(int(int(self.depth) / 2))
        svinfo = ';'.join(['PRECISE', svEnd, svType, svLEN])
        sample = ':'.join([GT, DR, DV])
        self.svLines.append('\t'.join([svChr, svPos, svID, aSeq, bSeq, svQual, svFilter, svinfo, "GT:DR:DV", sample]))

    def creatTRA(self,line):
        BND_types = {'forward:forward':{'A':{'A1B2':'N[P[','A2B1':']P]N','A3B4':'N[P[','A4B3':']P]N'},
                                        'B':{'B1A2':'N[P[','B2A1':']P]N','B3A4':'N[P[','B4A3':']P]N'}},
                     'forward:reverse':{'A':{'A1B3':'N]P]','A2B1':']P]N','A3B4':'N[P[','A4B2':'[P[N'},
                                        'B':{'B1A3':'N]P]','B2A1':']P]N','B3A4':'N[P[','B4A2':'[P[N'}},
                     'reverse:forward':{'A':{'A1B2':'N[P[','A2B4':'N[P[','A3B1':']P]N','A4B3':']P]N'},
                                        'B':{'B1A2':'N[P[','B2A4':'N[P[','B3A1':']P]N','B4A3':']P]N'}},
                     'reverse:reverse':{'A':{'A1B3':'N]P]','A2B4':'N[P[','A3B1':']P]N','A4B2':'[P[N'},
                                        'B':{'B1A3':'N]P]','B2A4':'N[P[','B3A1':']P]N','B4A2':'[P[N'}}}
        line4 = line[4].split(':')
        CHR1 = str(line[0])
        CHR2 = str(line4[1])
        CHR1A2 = int(line[1]) 
        CHR1A1 = CHR1A2 - 1
        CHR1A3 = int(line[2])
        CHR1A4 = CHR1A3 + 1
        CHR2B2 = int(line4[2])
        CHR2B1 = CHR2B2 - 1
        CHR2B3 = CHR2B2 + (CHR1A3-CHR1A2) 
        CHR2B4 = CHR2B3 + 1
        BND_type = ':'.join(line4[-2:])
        dataInfo = {'A':CHR1,'B':CHR2,'A1':str(CHR1A1),'A2':str(CHR1A2),'A3':str(CHR1A3),'A4':str(CHR1A4),'B1':str(CHR2B1),'B2':str(CHR2B2),'B3':str(CHR2B3),'B4':str(CHR2B4)}
        for k,v in BND_types[BND_type].items():
            for keyID,Value in v.items(): 
                CHR1X = keyID[:2]
                CHR2X = keyID[-2:]
                svChr = dataInfo[CHR1X[0]]
                svChr2 = dataInfo[CHR2X[0]]
                svPos = dataInfo[CHR1X]
                svID = "Sim."+str(self.nus)
                aSeq = self.genome.fetch(svChr, int(svPos), int(svPos) + 1)
                bSeq = Value
                svQual = '.'
                svFilter = 'PASS'
                svEnd = 'END='+dataInfo[CHR2X]
                svType = "SVTYPE=BND;CHR2="+svChr2
                svLEN = 'SVLEN=NULL'
                svinfo = ';'.join(['PRECISE', svEnd, svType, svLEN])
                DR, DV =  '5', '5'
                if self.Genotype[svChr] == self.Genotype[svChr2]:
                    GT = self.Genotype[svChr]
                else:
                    GT = '0/1'
                sample = ':'.join([GT, DR, DV])
                self.svLines.append('\t'.join([svChr, svPos, svID, aSeq, bSeq, svQual, svFilter, svinfo, "GT:DR:DV", sample]))
                self.nus+=1
    
    def OutPut(self):
        os.makedirs('/'.join(self.outFile.split('/')[:-1]), exist_ok=True)
        out = open(self.outFile,"w")
        print(self.outFile)
        for line in self.vcfHead+self.svLines:
            out.write(line+'\n')
        out.close()
        #os.system("cat %s|vcf-sort|bgzip -c >%s"%(self.outFile,self.outFile+'.gz'))
        #os.system("tabix %s"%self.outFile+'.gz')

    def run(self):
        try: 
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
        except BaseException as e:
            print(e)
        try:
            self.OutPut()
        except BaseException as e:
            print(e)

