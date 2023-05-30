import pysam
import os
import re


class createSimulateSV():
    def __init__(self,svType=None,outFile=None,callSVFile=None,ref=None):
        self.svType = svType
        self.outFile = outFile
        self.callSVFile = callSVFile
        self.vcfHead = open("Head.info", "r").read().strip().split("\n")
        self.SvLines = []
        self.SVIDprefix = '-'.join(['Sim',self.svType])
        self.genome = pysam.FastaFile(ref)
        self.nus = 1

    def getChrID(self,line):
        return str(line.chrom)

    def getPos(self,line):
        return str(line.pos)

    def getSVID(self,line):
        return  str(line.id)

    def getRE(self,line):
        tmp = [i for i in str(line).split('\t')[7].split(';') if 'SUPP=' in i]
        return int(tmp[0].split('=')[1])

    def getPASS(self,line):
        line =str(line)
        if 'PASS' in line:
            return False
        else:
            return True

    def getASeq(self,line):
        a_seq = 'N'
        return a_seq

    def getBSeq(self,line):
        b_seq = str(line.alts[0]).upper()
        if ('[' in b_seq):
            b_seq = re.sub('[A,T,C,G,N]','N',b_seq)
            bb = b_seq.split('[')
            bb[1] = 'P'
            return '['.join(bb)
        elif (']' in b_seq):
            b_seq = re.sub('[A,T,C,G,N]','N',b_seq)
            bb = b_seq.split(']')
            bb[1] = 'P'
            return ']'.join(bb)
        else:
            return b_seq

    def getSVEND(self,line):
        return 'END=' + str(line.stop)

    def getSVType(self,line):
        return 'SVTYPE=' + str(line.info['SVTYPE'])

    def getQual(self,line):
        return '.'

    def getFilter(self,line):
        return 'PASS'

    def getCHR2Pos(self,line):
        b_seq = line.alts[0]
        lines = str(line).split('\t')[7]
        if 'CHR2' in lines:
            return [[i for i in lines.split(';') if 'CHR2=' in i ][0].split('=')[1],str(line.stop)]
        elif '[' in b_seq:
            return b_seq.split('[')[1].split(':')
        elif ']':
            return b_seq.split(']')[1].split(':')


    def getInfo(self,line):
        svTYPE = 'SVTYPE=' + str(line.info['SVTYPE'])
        #svEND = 'END=' + str(line.stop)
        svEND = 'END=' + str(line.stop)
        try:
            try:
                svLEN = 'SVLEN=' + str(line.info['SVLEN'][0])
                if str(line.info['SVTYPE']) in ["BND", "TRA"]:
                    svLEN = 'SVLEN=' + "."
                    CHR2Pos = self.getCHR2Pos(line)
                    svTYPE = "SVTYPE=" + "TRA;CHR2="+CHR2Pos[0]
                    svEND = 'END=' + CHR2Pos[1]
            except TypeError:
                svLEN = 'SVLEN=' + str(line.info['SVLEN'])
                if str(line.info['SVTYPE']) in ["BND", "TRA"]:
                    svLEN = 'SVLEN=' + "."
                    CHR2Pos = self.getCHR2Pos(line)
                    svTYPE = "SVTYPE=" + "TRA;CHR2="+CHR2Pos[0]
                    svEND = 'END=' + CHR2Pos[1]
            infos = ';'.join(['PRECISE', svEND, svTYPE, svLEN])
        except KeyError:
            if str(line.info['SVTYPE']) == "INS":
                 svLEN = 'SVLEN=' + str(len(line.alts[0]))
            elif str(line.info['SVTYPE']) in ["DEL","INV","DUP"]:
                svLEN = 'SVLEN=' + str(int(self.getSVEND(line).split('=')[-1]) - int(self.getPos(line)))
            elif str(line.info['SVTYPE']) in  ["BND","TRA"]:
                svLEN = 'SVLEN=' + "."
                CHR2Pos = self.getCHR2Pos(line)
                print(CHR2Pos)
                svTYPE = "SVTYPE=" + "TRA;CHR2="+CHR2Pos[0]
                svEND = 'END=' + CHR2Pos[1]
            infos = ';'.join(['PRECISE', svEND, svTYPE, svLEN])
        return infos

    def getSample(self,line):
        lines = str(line).strip().split('\t')
        if str(line.info['SVTYPE']) in  ["BND","TRA"]:
            GT = '/'.join([str(ix) for ix in line.samples.values()[0]['GT']])
            if "./." == GT:
                GT = '0/1'
        else:
            GTs = [lines[-2].split(':')[0], lines[-1].split(':')[0]]
            if "./." in GTs:
                GTs.remove('./.')
            if len(GTs)==1:
                GT = '0/1'
            else:
                GT = GTs[1]
        DR = '5'
        DV = '5'
        sample = ':'.join([GT, DR, DV])
        #print('GT:',GT)
        return sample,GT

    def OutPut(self):
        os.makedirs('/'.join(self.outFile.split('/')[:-1]), exist_ok=True)
        print('/'.join(self.outFile.split('/')[:-1]))
        out = open(self.outFile,"w")
        for line in self.vcfHead+self.SvLines:
            out.write(line+'\n')
        out.close()
        #os.system("cat %s|vcf-sort|bgzip -c >%s"%(self.outFile,self.outFile+'.gz'))
        #os.system("tabix %s"%self.outFile+'.gz')

    def run(self):
        IDnu = 1
        chrsID = [str(i) for i in range(1, 23, 1)] + ['X', 'Y', 'MT']
        for line in pysam.VariantFile(self.callSVFile):
            if ('#' in str(line)) or (str(line.chrom) not in chrsID): \
               #(str(line.info['SVTYPE']) != self.svType) or \
#               (str(line.chrom) not in chrsID): #or \
               #(str(line.filter.keys()[0]) not in  ["PASS"]):
               #E print(line)
                continue
            else:
                tmpSample, tmpGT = self.getSample(line)
                lines = str(line).strip().split('\t')
                self.SvLines.append('\t'.join([self.getChrID(line), \
                                             self.getPos(line), \
                                             self.SVIDprefix+'.'+str(IDnu), \
                                             self.getASeq(line),
                                             self.getBSeq(line),
                                             self.getQual(line),
                                             self.getFilter(line),
                                             self.getInfo(line),
                                             'GT:DR:DV',
                                             tmpSample]))
            IDnu=IDnu+1
        self.OutPut()

