import pysam
import sys
import os

class EvaluationSV():
    def __init__(self, VcfBaseFile=None, VcfCommFile=None, svType=None, outFile=None):
         #self.vcfbase = pysam.VariantFile(VcfBaseFile)
         #self.vcfcomm = pysam.VariantFile(VcfCommFile)
         self.location = 1000
         self.svType = svType
         self.outFile = outFile
         self.vcfhead = "Head.info"
         self.TP_GT_comm = []
         self.TP_no_GT_comm = []
         self.TP_GT_base = []
         self.TP_no_GT_base = []
         self.lsv_GT = []
         self.bc1_GT = []
         self.bc2_GT = []
         self.lsv_no_GT = []
         self.bc1_no_GT = []
         self.bc2_no_GT = []
         self.lsv_base_GT_distr = {'0-50':0,'50-100':0,'100-200':0,'200-300':0,'300-500':0,'500-1000':0,'1000-2500':0,'2500-99999999999':0}
         self.lsv_comm_GT_distr = {'0-50':0,'50-100':0,'100-200':0,'200-300':0,'300-500':0,'500-1000':0,'1000-2500':0,'2500-99999999999':0}
         self.lsv_base_no_GT_distr = {'0-50':0,'50-100':0,'100-200':0,'200-300':0,'300-500':0,'500-1000':0,'1000-2500':0,'2500-99999999999':0}
         self.lsv_comm_no_GT_distr = {'0-50':0,'50-100':0,'100-200':0,'200-300':0,'300-500':0,'500-1000':0,'1000-2500':0,'2500-99999999999':0}
         self.VcfCommFile = VcfCommFile
         self.baseSize, self.vcfbase, self.vcfbaselins = self.readvcf(pysam.VariantFile(VcfBaseFile))
         self.commSize, self.vcfcomm, self.vcfcommlins = self.readvcf(pysam.VariantFile(VcfCommFile))
         
    def getsvLEN(self,line):
        try:
            try:
                svLEN =  abs(int(line.info['SVLEN'][0]))
            except TypeError:
                svLEN = abs(int(line.info['SVLEN']))
        except KeyError:
            if str(line.info['SVTYPE']) == "INS":
                 svLEN = abs(int(len(line.alts[0])))
            elif str(line.info['SVTYPE']) in ["DEL","INV","DUP"]:
                svLEN = abs(int(int(self.getSVEND(line).split('=')[-1]) - int(self.getPos(line))))
        return svLEN

    def readvcf(self, VcfFile):
        vcfDict = {}
        lineSize = 0
        vcfLines = []
        for line in VcfFile.fetch():
            lineSize+=1
            chrID = str(line.chrom)
            start = int(line.pos) 
            end   = int(line.stop)
            geneType =  '/'.join([str(ix) for ix in line.samples.values()[0]['GT']])
            b_seq = str(line.alts[0])
            vcfLines.append(str(line))

            if self.svType in  ['TRA','BND']:
                CHR2 = [i.split('=')[-1] for i in str(line).split('\t')[7].split(';') if 'CHR2=' in i][0]
                svlen = 0
            else:
                CHR2 = 'None'
                svlen = self.getsvLEN(line)
            if chrID in vcfDict:
                vcfDict[chrID].append([chrID, start, end, svlen, CHR2, geneType, b_seq, str(line)])
            else:
                vcfDict[chrID] = [[chrID, start, end, svlen, CHR2, geneType, b_seq, str(line)]]
        return lineSize, vcfDict, vcfLines

    def del_ins_inv_dup(self, vcfbaseline, vcfcommline):
        Result = False
        bc1 = vcfcommline[1] - vcfbaseline[1]
        bc2 = vcfcommline[2] - vcfbaseline[2]
        lsv = vcfcommline[3] - vcfbaseline[3]
        if max([vcfcommline[1]-self.location, vcfbaseline[1]]) <= min([vcfcommline[2]+self.location, vcfbaseline[2]]):
           if min([vcfbaseline[3], vcfcommline[3]])/max([vcfbaseline[3], vcfcommline[3]])>=0.7:
               Result = True
               if abs(bc1)>1000:
                   print(vcfbaseline)
                   print(vcfcommline)
                   print('bc1---------------------')
               elif abs(bc2)>1000:
                   print(vcfbaseline)
                   print(vcfcommline)
                   print('bc2---------------------')
        return Result, bc1, bc2, lsv

    def tra_bnd(self, vcfbaseline, vcfcommline):
        Result = False
        bc1 = vcfcommline[1] - vcfbaseline[1]
        bc2 = vcfcommline[2] - vcfbaseline[2]
        if abs(bc1) <= self.location:
            if abs(bc2) <= self.location:
                if (vcfbaseline[4] == vcfcommline[4]):
                    if 'debreak' in self.VcfCommFile:
                        Result = True
                    elif vcfbaseline[6] == vcfcommline[6]:
                        Result = True
        return Result, bc1, bc2

    def oneselcet(self,number):
        if number>=1:
            return 1
        else:
            return number

    def writeReslut(self, precision_GT,recall_GT,F1_GT,precision_no_GT,recall_no_GT,F1_no_GT):
        os.makedirs(self.outFile, exist_ok=True)
        headlines = open(self.vcfhead,'r').read().strip().split('\n')
        with open('%s/summary.txt'%self.outFile,'w') as out:
            out.write('precision_GT,'+str(precision_GT)+'\n')
            out.write('recall_GT,'+str(recall_GT)+'\n')
            out.write('F1_GT,'+str(F1_GT)+'\n')
            out.write('precision_no_GT,'+str(precision_no_GT)+'\n')
            out.write('recall_no_GT,'+str(recall_no_GT)+'\n')
            out.write('F1_no_GT,'+str(F1_no_GT)+'\n')
        with open('%s/bias.txt'%self.outFile,'w') as out:
            out.write('bc1_GT,'+','.join(list(map(str, self.bc1_GT)))+'\n')
            out.write('bc2_GT,'+','.join(list(map(str, self.bc2_GT)))+'\n')
            out.write('bc1_no_GT,'+','.join(list(map(str, self.bc1_no_GT)))+'\n')
            out.write('bc2_no_GT,'+','.join(list(map(str, self.bc2_no_GT)))+'\n')
            if self.svType in ['DEL','INS','INV','DUP']:
                out.write('lsv_GT,'+','.join(list(map(str, self.lsv_GT)))+'\n')
                out.write('lsv_no_GT,'+','.join(list(map(str, self.lsv_no_GT)))+'\n')

        with open('%s/TP_GT.vcf'%self.outFile,'w') as out:
            out.write('\n'.join(headlines)+'\n')
            for lx in self.TP_GT_comm:
                out.write(lx)
        with open('%s/FP_GT.vcf'%self.outFile,'w') as out:
            out.write('\n'.join(headlines)+'\n')
            for lx in self.vcfcommlins:
                if lx not in self.TP_GT_comm:
                    out.write(lx)
        with open('%s/FN_GT.vcf'%self.outFile,'w') as out:
            out.write('\n'.join(headlines)+'\n')
            for lx in self.vcfbaselins:
                if lx not in self.TP_GT_base:
                    out.write(lx)

        with open('%s/TP_no_GT.vcf'%self.outFile,'w') as out:
            out.write('\n'.join(headlines)+'\n')
            for lx in self.TP_no_GT_comm:
                out.write(lx)
        with open('%s/FP_no_GT.vcf'%self.outFile,'w') as out:
            out.write('\n'.join(headlines)+'\n')
            for lx in self.vcfcommlins:
                if lx not in self.TP_no_GT_comm:
                    out.write(lx)
        with open('%s/FN_no_GT.vcf'%self.outFile,'w') as out:
            out.write('\n'.join(headlines)+'\n')
            for lx in self.vcfbaselins:
                if lx not in self.TP_no_GT_base:
                    out.write(lx)

        with open('%s/lsv.txt'%self.outFile,'w') as out:
            out.write(','.join(['ID','base_GT','comm_GT','base_no_GT','comm_no_GT'])+'\n')
            for k,v in self.lsv_base_GT_distr.items():
                out.write(','.join([k,str(v),str(self.lsv_comm_GT_distr[k]),str(self.lsv_base_no_GT_distr[k]),str(self.lsv_comm_no_GT_distr[k])])+'\n')
                

    def F1_recall_precision(self):
        try:
            precision_GT = self.oneselcet(len(set(self.TP_GT_comm))/float(self.commSize))
            recall_GT =  self.oneselcet(len(set(self.TP_GT_base))/float(self.baseSize))
            F1_GT = 2*precision_GT*recall_GT/(precision_GT+recall_GT)
        except ZeroDivisionError:
            precision_GT,recall_GT,F1_GT = 0,0,0

        try:
            precision_no_GT = self.oneselcet(len(set(self.TP_no_GT_comm))/float(self.commSize))
            recall_no_GT =  self.oneselcet(len(set(self.TP_no_GT_base))/float(self.baseSize))
            F1_no_GT = 2*precision_no_GT*recall_no_GT/(precision_no_GT+recall_no_GT)
        except ZeroDivisionError:
            precision_no_GT,recall_no_GT,F1_no_GT = 0,0,0

        self.writeReslut(precision_GT,recall_GT,F1_GT,precision_no_GT,recall_no_GT,F1_no_GT)
        
        #print(self.bc1_no_GT)
        print(self.outFile)
        #print("precision_GT:",precision_GT)
        #print("recall_GT:",recall_GT)
        #print("F1_GT:",F1_GT)
        #print("precision_no_GT:",precision_no_GT)
        #print("recall_no_GT:",recall_no_GT)
        #print("F1_no_GT:",F1_no_GT)
        #print(self.bc2_no_GT)

    def lsv_Dict(self,linex,lsvType):
        lsv = int(linex[3])
        for k,v in self.lsv_base_GT_distr.items():
            kk = k.split('-')
            if int(kk[0])<int(lsv) <=int(kk[1]):
                if lsvType=='base_no_GT':
                    self.lsv_base_no_GT_distr[k]+=1
                elif lsvType=='comm_no_GT':
                    self.lsv_comm_no_GT_distr[k]+=1
                elif lsvType=='base_GT':
                    self.lsv_base_GT_distr[k]+=1
                elif lsvType=='comm_GT':
                    self.lsv_comm_GT_distr[k]+=1

    def base_comm(self):
        for chrx, linesbase in self.vcfbase.items():
            for  linebase in linesbase:
                if chrx in self.vcfcomm:
                    for linecomm in self.vcfcomm[chrx]:
                        if (self.svType in ['DEL', 'INS', 'DUP', 'INV']):
                            Result, bc1, bc2, lsv = self.del_ins_inv_dup(linebase,linecomm)
                            if Result:
                                self.TP_no_GT_comm.append(linecomm[-1])
                                self.TP_no_GT_base.append(linebase[-1])
                                self.lsv_Dict(linebase,'base_no_GT')
                                self.lsv_Dict(linecomm,'comm_no_GT')
                                self.bc1_no_GT.append(bc1)
                                self.bc2_no_GT.append(bc2)
                                self.lsv_no_GT.append(lsv)           
                                if linebase[5] == linecomm[5]:
                                    self.TP_GT_comm.append(linecomm[-1])
                                    self.TP_GT_base.append(linebase[-1])
                                    self.lsv_Dict(linebase,'base_GT')
                                    self.lsv_Dict(linecomm,'comm_GT')
                                    self.bc1_GT.append(bc1)
                                    self.bc2_GT.append(bc2)
                                    self.lsv_GT.append(lsv)
                        elif (self.svType in ['TRA', 'BND']):
                            Result, bc1, bc2 = self.tra_bnd(linebase,linecomm)
                            if Result:
                                self.TP_no_GT_comm.append(linecomm[-1])
                                self.TP_no_GT_base.append(linebase[-1])
                                self.bc1_no_GT.append(bc1)
                                self.bc2_no_GT.append(bc2)
                                if linebase[5] == linecomm[5]:
                                    self.TP_GT_comm.append(linecomm[-1])
                                    self.TP_GT_base.append(linebase[-1])
                                    self.bc1_GT.append(bc1)
                                    self.bc2_GT.append(bc2)
        self.F1_recall_precision()        
          
    def run(self):
        self.base_comm()
#    def base_comm(self):
#        for chrx, linesbase in self.vcfbase.items():
#            for  linebase in linesbase:           
#                if chrx in self.vcfcomm:
#                    for linecomm in self.vcfcomm[chrx]:
#                        if (self.svType in ['DEL', 'INS', 'DUP', 'INV']) and (self.del_ins_inv_dup(linebase,linecomm)):
#                            self.TP_no_GT_comm.append(linecomm[-1])
#                            self.TP_no_GT_base.append(linebase[-1])
#                            if linebase[5] == linecomm[5]:
#                                self.TP_GT_comm.append(linecomm[-1])
#                                self.TP_GT_base.append(linebase[-1])
#                        elif (self.svType in ['TRA', 'BND']) and (self.tra_bnd(linebase,linecomm)):
#                            self.TP_no_GT_comm.append(linecomm[-1])
#                            self.TP_no_GT_base.append(linebase[-1])
#                            if linebase[5] == linecomm[5]:
#                                self.TP_GT_comm.append(linecomm[-1]) 
#                                self.TP_GT_base.append(linebase[-1])
#        self.F1_recall_precision()

if __name__ == '__main__':
    evalsv = EvaluationSV(VcfBaseFile=sys.argv[1], VcfCommFile=sys.argv[2], svType=sys.argv[3], outFile="tet")
    evalsv.base_comm()
    
