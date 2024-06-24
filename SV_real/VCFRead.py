import pysam
import os
import re

class VCFRead():
    def __init__(self, platform=None, depth=None, svType=None, vcfFile=None,outFile=None,
                 ref="data/hs37d5.fa",mapTool=None, svTool=None, RE_nu=2, nux=None):
        self.platform = platform
        self.depth = depth
        self.svType = svType
        self.vcfFile = vcfFile
        self.genome = pysam.FastaFile(ref)
        self.nux = nux
        self.RE_nu = RE_nu
        self.outFile = outFile
        self.mapTool = mapTool
        self.svTool = svTool

        self.tmpFileName = 'tmp.%s.%s.%s.%s.%s.%s.vcf'%(str(self.RE_nu),str(self.nux),platform,depth,mapTool,svTool)
        self.writetmp()
        self.vcfIn = pysam.VariantFile(self.tmpFileName)
        #print('tmp.%s.vcf'%str(nux))
        #self.vcfHead = str(self.vcfIn.header).strip().split('\n')
        self.vcfHead = open("Head.info", "r").read().strip().split("\n")
        self.SvLines = []
        print(self.vcfFile)
        if svType=='BND':
            self.SVIDprefix = '-'.join([platform,depth,'TRA',mapTool,svTool])
        else:
            self.SVIDprefix = '-'.join([platform,depth,svType,mapTool,svTool])

    def writetmp(self):
        with open(self.tmpFileName,'w') as out:
            if self.svTool == 'picky':
                out.write(re.sub("INFO\tFORMAT.*","INFO\tFORMAT\tNUL",open('Head.picky.info','r').read()))
                with open(self.vcfFile, 'r') as inp:
                    for li in inp:
                        if "#" in li:
                            continue
                        elif '=>' in li:
                            li = li.replace("=>",'=')
                            out.write(li)
                        else:
                            out.write(li)

            else:
                with open(self.vcfFile, 'r') as inp:
                    for li in inp:
                        if ("#" in li) and ("INFO\tFORMAT" in li):
                            li = re.sub("INFO\tFORMAT.*","INFO\tFORMAT\tNUL",li)
                            out.write(li)
                        elif "#" in li:
                            out.write(li)
                        elif '=>' in li:
                            li = li.replace("=>",'=')
                            out.write(li)
                        else:
                            out.write(li)
    
    def getChrID(self,line):
        return str(line.chrom)

    def getPos(self,line):
        return str(line.pos)

    def getSVID(self,line):
        return  str(line.id)
 
    def getRE(self,line):
        tmp = [i for i in str(line).split('\t')[7].split(';') if 'RE=' in i]
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
        b_seq = str(line.alts[0])#.upper()
        if ('[' in b_seq):
            return b_seq
            #b_seq = re.sub('[A,T,C,G,N]','N',b_seq)
            #bb = b_seq.split('[')
            #bb[1] = 'P'
            #return '['.join(bb)
        elif (']' in b_seq):
            return b_seq
            #b_seq = re.sub('[A,T,C,G,N]','N',b_seq)
            #bb = b_seq.split(']')
            #bb[1] = 'P'
            #return ']'.join(bb)
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
        if '[' in b_seq:
            return b_seq.split('[')[1].split(':')
        elif ']':
            return b_seq.split(']')[1].split(':')
        elif 'CHR2' in lines:
            return [[i for i in lines.split(';') if 'CHR2' in i ][0].split('=')[1],str(line.stop)]

    def getInfo(self,line):
        svTYPE = 'SVTYPE=' + str(self.get_line_SVTYPE(line))
        svEND = self.getSVEND(line)
        try:
            try:
                svLEN = 'SVLEN=' + str(line.info['SVLEN'][0])
                if str(line.info['SVTYPE']) in ["BND", "TRA"]:
                    svLEN = 'SVLEN=' + "NULL"
                    CHR2Pos = self.getCHR2Pos(line)
                    svTYPE = "SVTYPE=" + "BND;CHR2="+CHR2Pos[0]
                    svEND = 'END=' + CHR2Pos[1]
            except TypeError:
                svL = str(line.info['SVLEN'])
                svLEN = 'SVLEN=' + svL
                if str(line.info['SVTYPE']) in ["BND", "TRA"]:
                    svLEN = 'SVLEN=' + "NULL"
                    CHR2Pos = self.getCHR2Pos(line)
                    svTYPE = "SVTYPE=" + "BND;CHR2="+CHR2Pos[0]
                    svEND = 'END=' + CHR2Pos[1]
            infos = ';'.join(['PRECISE', svEND, svTYPE, svLEN])
        except KeyError:
            if str(line.info['SVTYPE']) == "INS":
                 svLEN = 'SVLEN=' + str(len(line.alts[0]))
            elif str(line.info['SVTYPE']) in ["DEL","INV","DUP"]:
                svLEN = 'SVLEN=' + str(int(self.getSVEND(line).split('=')[-1]) - int(self.getPos(line)))
            elif str(line.info['SVTYPE']) in ["BND", "TRA"]:
                svLEN = 'SVLEN=' + "NULL"
                CHR2Pos = self.getCHR2Pos(line)
                svTYPE = "SVTYPE=" + "BND;CHR2="+CHR2Pos[0]
                svEND = 'END=' + CHR2Pos[1]
            infos = ';'.join(['PRECISE', svEND, svTYPE, svLEN])
        return infos

    def getSample(self,line):
        GT = '/'.join([str(ix) for ix in line.samples.values()[0]['GT']])
        DR = '5'#line.samples.values()[0]['DR']
        DV = '5'#line.samples.values()[0]['DV']
        sample = ':'.join([GT, DR, DV])
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

    def chSVTYPE(self,line):
        if line.info['SVTYPE']=='TRA':
            return 'BND'
        else:
            
            return str(line.info['SVTYPE'])

    def get_passinfo(self,line):
        #print(line.filter.keys())
        passinfo = str(line.filter.keys()[0])
        return passinfo

    def get_line_SVTYPE(self,line):
        line_SVTYPE = self.chSVTYPE(line)
        if self.svType in line_SVTYPE:
            line_SVTYPE = self.svType
        return line_SVTYPE
    
    def run(self):
        IDnu = 1
        chrsID = ["chr"+str(i) for i in range(1, 23, 1)] + ['chrX', 'chrY', 'chrMT']
        for line in self.vcfIn:
            passinfo = self.get_passinfo(line)
            line_SVTYPE = self.get_line_SVTYPE(line)
            #print('a',passinfo,line_SVTYPE,self.svType,line.info['SVTYPE'])
            if (line_SVTYPE != self.svType) or \
                    (str(line.chrom) not in chrsID) or \
                    (passinfo not in  ["PASS"]):
                    #(str(line.filter.keys()[0]) not in  ["PASS"]):
                continue
            else:
                if (self.svType in ['TRA','BND']) and (self.getCHR2Pos(line)[0] not in chrsID):
                    continue
                tmpSample, tmpGT = self.getSample(line)
                if "NanoSV" in self.vcfFile and not tmpGT:
                    continue
                if self.svType == 'DUP':
                    if tmpGT == "None/None" :
                        if self.getRE(line)>=self.RE_nu:
                            tmpSample = tmpSample.replace("None/None",'0/1')
                            self.SvLines.append('\t'.join([self.getChrID(line),
                                                           self.getPos(line),
                                                           self.SVIDprefix+'.'+str(IDnu),                       #self.getSVID(line),
                                                           self.getASeq(line),
                                                           self.getBSeq(line),
                                                           self.getQual(line),
                                                           self.getFilter(line),
                                                           self.getInfo(line),
                                                           'GT:DR:DV',
                                                           tmpSample]))
                            IDnu+=1
                    
                    elif tmpGT != "None/None" and self.getRE(line)>=self.RE_nu:
                        self.SvLines.append('\t'.join([self.getChrID(line),
                                                       self.getPos(line),
                                                       self.SVIDprefix+'.'+str(IDnu),                       #self.getSVID(line),
                                                       self.getASeq(line),
                                                       self.getBSeq(line),
                                                       self.getQual(line),
                                                       self.getFilter(line),
                                                       self.getInfo(line),
                                                       'GT:DR:DV',
                                                       tmpSample]))
                        IDnu+=1
                    else:
                        continue
                else:
                    if tmpGT == "None/None" and 'picky' in self.vcfFile :
                        if self.getRE(line)>=self.RE_nu:
                            tmpSample = tmpSample.replace("None/None",'0/1')
                            self.SvLines.append('\t'.join([self.getChrID(line),
                                                           self.getPos(line),
                                                           self.SVIDprefix+'.'+str(IDnu),                       #self.getSVID(line),
                                                           self.getASeq(line),
                                                           self.getBSeq(line),
                                                           self.getQual(line),
                                                           self.getFilter(line),
                                                           self.getInfo(line),
                                                           'GT:DR:DV',
                                                           tmpSample]))
                            IDnu+=1
 
                    elif tmpGT != "None/None" and self.getRE(line)>=self.RE_nu:
                
                        self.SvLines.append('\t'.join([self.getChrID(line),
                                                       self.getPos(line),
                                                       self.SVIDprefix+'.'+str(IDnu),                       #self.getSVID(line),
                                                       self.getASeq(line),
                                                       self.getBSeq(line),
                                                       self.getQual(line),
                                                       self.getFilter(line),
                                                       self.getInfo(line),
                                                       'GT:DR:DV',
                                                       tmpSample]))
                        IDnu+=1                     
        self.OutPut()
        os.remove(self.tmpFileName)
        print('aaa')
