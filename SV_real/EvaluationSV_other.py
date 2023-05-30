import pysam
import sys
import os

class EvaluationSV():
    def __init__(self, VcfBaseFile=None, VcfCommFile=None, svType=None, outFile=None, RE_nu= None, nux = None):
         #self.vcfbase = pysam.VariantFile(VcfBaseFile)
         #self.vcfcomm = pysam.VariantFile(VcfCommFile)
         self.location = 1000
         self.RE_nu = RE_nu
         self.nux = nux
         self.svType = svType
         self.outFile = outFile
         self.vcfhead = "Head.info"
         self.TP_GT_comm = {'50-100':[],'100-500':[],'500-1000':[],'1000-2500':[],'2500-99999999999':[]}
         self.TP_no_GT_comm = {'50-100':[],'100-500':[],'500-1000':[],'1000-2500':[],'2500-99999999999':[]}
         self.TP_GT_base = {'50-100':[],'100-500':[],'500-1000':[],'1000-2500':[],'2500-99999999999':[]}
         self.TP_no_GT_base = {'50-100':[],'100-500':[],'500-1000':[],'1000-2500':[],'2500-99999999999':[]}
         self.lsv_GT = {'50-100':[],'100-500':[],'500-1000':[],'1000-2500':[],'2500-99999999999':[]}
         self.bc1_GT = {'50-100':[],'100-500':[],'500-1000':[],'1000-2500':[],'2500-99999999999':[]}
         self.bc2_GT = {'50-100':[],'100-500':[],'500-1000':[],'1000-2500':[],'2500-99999999999':[]}
         self.lsv_no_GT = {'50-100':[],'100-500':[],'500-1000':[],'1000-2500':[],'2500-99999999999':[]}
         self.bc1_no_GT = {'50-100':[],'100-500':[],'500-1000':[],'1000-2500':[],'2500-99999999999':[]}
         self.bc2_no_GT = {'50-100':[],'100-500':[],'500-1000':[],'1000-2500':[],'2500-99999999999':[]}
         self.VcfCommFile = VcfCommFile
         #self.bedDict = self.bedfile("/home/lz/Data_sequence/ref/annotions/giab/HG002_SVs_Tier1_v0.6.bed")
         bedDicts = {'BND':"/home/lz/Data_sequence/2021_4_2/SV_static/SV_sim/bed/BND.bed",
                     'TRA':"/home/lz/Data_sequence/2021_4_2/SV_static/SV_sim/bed/BND.bed",
                     'DEL':"/home/lz/Data_sequence/work/Test/Truvari/data/HG002_SVs_Tier1_v0.6.bed",
                     'INS':"/home/lz/Data_sequence/work/Test/Truvari/data/HG002_SVs_Tier1_v0.6.bed",
                     'INV':"/home/lz/Data_sequence/2021_4_2/SV_static/SV_HG002/bed/INV.bed",
                     'DUP':"/home/lz/Data_sequence/2021_4_2/SV_static/SV_HG002/bed/DUP.bed"}
         self.bedDict = self.bedfile(bedDicts[svType])
         self.baseSize, self.vcfbase, self.vcfbaselins, self.lenbase, self.nnuubaseDict = self.readvcf(pysam.VariantFile(VcfBaseFile))
         self.commSize, self.vcfcomm, self.vcfcommlins, self.lencomm, self.nnuucommDict = self.readvcf(pysam.VariantFile(VcfCommFile))
         print(self.lenbase)
         print(self.lencomm)
         print(self.baseSize)
         print(self.commSize)

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

    def bedfile(self,bedfile):
        bedDict ={}
        for line in open(bedfile,'r'):
            line = line.split('\t')
            a = int(line[1])
            b = int(line[2])
            if a > b:
                if line[0] in bedDict:
                    bedDict[line[0]].append([b, a])
                else:
                    bedDict[line[0]] = [[b, a]]
            else:
                if line[0] in bedDict:
                    bedDict[line[0]].append([a, b])
                else:
                    bedDict[line[0]] = [[a, b]]
        return bedDict

    def bedPass(self, line ):
        if (line[3]) < 50:
            return False
        chrid = line[0]
        chr2 = line[4]
        try:
            chridList = self.bedDict[chrid]
        except:
            return False

        if chr2 == "None":
            for lx in chridList:
                if (lx[0]-1000 <= int(line[1]) <= lx[1]+1000) or (lx[0]-1000 <= int(line[2]) <= lx[1]+1000):
                    return True
        else:
            try:
                chridList2 = self.bedDict[chr2]
            except:
                return False
            for lx in chridList:
                for lxx in chridList2:
                    if (lx[0]-1000 <= int(line[1]) <= lx[1]+1000) or (lxx[0]-1000 <= int(line[2]) <= lxx +1000):
                        return True
        return False

    def readvcf(self, VcfFile):
        vcfDict = {}
        lineSize = 0
        vcfLines = []
        tmpDict = {'50-100':0,'100-500':0,'500-1000':0,'1000-2500':0,'2500-99999999999':0}
        nnuu = 0
        nnuuDict = {}
        for line in VcfFile.fetch():
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
                if self.bedPass([chrID, start, end, svlen, CHR2, geneType, b_seq,"ID"+str(nnuu), str(line)]):
                    vcfDict[chrID].append([chrID, start, end, svlen, CHR2, geneType, b_seq,"ID"+str(nnuu), str(line)])
                    vcfLines.append(str(line))
                    lineSize+=1
                else:
                    continue
            else:
                if self.bedPass([chrID, start, end, svlen, CHR2, geneType, b_seq,"ID"+str(nnuu), str(line)]):
                    vcfDict[chrID] = [[chrID, start, end, svlen, CHR2, geneType, b_seq,"ID"+str(nnuu), str(line)]]
                    vcfLines.append(str(line))
                    lineSize+=1
                else:
                    continue
            
            for kk in tmpDict.keys():
                k = kk.strip().split('-')
                s, e = int(k[0]), int(k[1])
                if s <= svlen < e:
                    tmpDict[kk] += 1
                    break
                    
            nnuuDict["ID"+str(nnuu)] = str(line)
            nnuu+=1
        return lineSize, vcfDict, vcfLines, tmpDict, nnuuDict

    def del_ins_inv_dup(self, vcfbaseline, vcfcommline):
        Result = False
        tmp = ['50-100','100-500','500-1000','1000-2500','2500-99999999999']
        bc1 = vcfcommline[1] - vcfbaseline[1]
        bc2 = vcfcommline[2] - vcfbaseline[2]
        lsv = vcfcommline[3] - vcfbaseline[3]
        if max([vcfcommline[1]-self.location, vcfbaseline[1]]) <= min([vcfcommline[2]+self.location, vcfbaseline[2]]):
           if min([vcfbaseline[3], vcfcommline[3]])/max([vcfbaseline[3], vcfcommline[3]])>=0.7:
               Result = True
        for ix in tmp:
            i=ix.strip().split('-')
            s, e = int(i[0]), int(i[1])
            if s <= vcfcommline[3] < e:
                rangeId = ix
                break
        return Result, bc1, bc2, lsv, rangeId

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

    def writeReslut(self, F1_recall_precisionGTDict,F1_recall_precisionNGTDict):
        os.makedirs(self.outFile, exist_ok=True)
        headlines = open(self.vcfhead,'r').read().strip().split('\n')
        with open('%s/summary.txt'%self.outFile,'w') as out:
            out.write('precision_GT,'+str(F1_recall_precisionGTDict["all"][0])+'\n')
            out.write('recall_GT,'+str(F1_recall_precisionGTDict["all"][1])+'\n')
            out.write('F1_GT,'+str(F1_recall_precisionGTDict["all"][2])+'\n')
            out.write('precision_no_GT,'+str(F1_recall_precisionNGTDict["all"][0])+'\n')
            out.write('recall_no_GT,'+str(F1_recall_precisionNGTDict["all"][1])+'\n')
            out.write('F1_no_GT,'+str(F1_recall_precisionNGTDict["all"][2])+'\n')
            out.write('TP_GT_comm,'+str(sum([len(set(i)) for i in self.TP_GT_comm.values()]))+'\n')
            out.write('TP_GT_base,'+str(sum([len(set(i)) for i in self.TP_GT_base.values()]))+'\n')
            out.write('TP_no_GT_comm,'+str(sum([len(set(i)) for i in self.TP_no_GT_comm.values()]))+'\n')
            out.write('TP_no_GT_base,'+str(sum([len(set(i)) for i in self.TP_no_GT_base.values()]))+'\n')
            out.write('commSize,'+str(sum(self.lencomm.values()))+'\n')
            out.write('baseSize,'+str(sum(self.lenbase.values()))+'\n')
        with open('%s/bias.txt'%self.outFile,'w') as out:
            out.write('bc1,'+','.join(list(map(str, [ii for i in self.bc1_no_GT.values() for ii in i])))+'\n')
            out.write('bc2,'+','.join(list(map(str, [ii for i in self.bc2_no_GT.values() for ii in i])))+'\n')
            if self.svType in ['DEL','INS','INV','DUP']:
                out.write('lsv,'+','.join(list(map(str, [ii for i in self.lsv_no_GT.values() for ii in i])))+'\n')

        with open('%s/lsv_F1.txt'%self.outFile,'w') as out:
            #out.write(','.join(['ID','base_GT','comm_GT','base_no_GT','comm_no_GT'])+'\n')
            out.write(','.join(['ID','baseT','commT','baseAll','commAll','precision','recall','F1'])+'\n')
            for k in self.bc1_no_GT.keys():
                out.write(','.join([k,str(len(set(self.TP_no_GT_base[k]))),str(len(set(self.TP_no_GT_comm[k]))),
                                      str(self.lenbase[k]),str(self.lencomm[k]),
                                      str(F1_recall_precisionNGTDict["precision"][k]),str(F1_recall_precisionNGTDict["recall"][k]),str(F1_recall_precisionNGTDict['F1'][k])])+'\n')

        with open('%s/lsv_bias.txt'%self.outFile,'w') as out:
            for k in self.bc1_no_GT.keys():
                out.write(','.join([k,'bc1',','.join(list(map(str,self.bc1_no_GT[k])))])+'\n')
                out.write(','.join([k,'bc2',','.join(list(map(str,self.bc2_no_GT[k])))])+'\n')
                out.write(','.join([k,'lsv',','.join(list(map(str,self.lsv_no_GT[k])))])+'\n')
        
        with open("%s/tp-call.vcf"%self.outFile,'w') as out:
            out.write(open(self.vcfhead,'r').read())
            for lx in self.TP_no_GT_comm.values():
                out.write('\n'.join([self.nnuucommDict[lxx] for lxx in lx]))
        
        with open("%s/tp-base.vcf"%self.outFile,'w') as out:
            out.write(open(self.vcfhead,'r').read())
            for lx in self.TP_no_GT_base.values():
                out.write('\n'.join([self.nnuubaseDict[lxx] for lxx in lx]))
                
        with open("%s/fn.vcf"%self.outFile,'w') as out:
            out.write(open(self.vcfhead,'r').read())
            for lx in (set(self.nnuubaseDict.keys()) - set([ii for i in self.TP_no_GT_base.values() for ii in i])):
                out.write(self.nnuubaseDict[lx])
                    
        with open("%s/fp.vcf"%self.outFile,'w') as out:
            out.write(open(self.vcfhead,'r').read())
            for lx in (set(self.nnuucommDict.keys()) - set([ii for i in self.TP_no_GT_comm.values() for ii in i])):
                out.write(self.nnuucommDict[lx])

    #计算不同长度分布SV的precision recall F1
    def lenDist(self, lenDictCommTP, lenDictCommAll, lenDictBaseTP, lenDictBaseAll):
        tmp_F1 = self.tmpDict.copy()
        tmp_recall = tmp_F1.copy()
        tmp_precision = tmp_F1.copy()
        for k in tmp_F1.keys():
            try:
                tmp_precision[k] = int(lenDictCommTP[k])/float(lenDictCommAll[k])
                tmp_recall[k] = int(lenDictBaseTP[k])/float(lenDictBaseAll[k])
                tmp_F1[k] = 2*tmp_precision[k]*tmp_recall[k]/(tmp_precision[k] + tmp_recall[k])
            except ZeroDivisionError:
                tmp_precision[k], tmp_recall[k], tmp_F1[k] = 0, 0, 0
        return tmp_precision, tmp_recall, tmp_F1

    def lsv_bc(self,gs):
        bc1Dict = dict(zip(self.tmpDict.keys(),[[] for i in self.tmpDict.keys()]))
        bc2Dict = bc1Dict.copy()
        lsvDict = bc1Dict.copy()
        bc1s = []
        bc2s = []
        lsvs = []
        for line in gs:
            line = line.strip().split(':')
            GT_base = line[0]
            SVID_base = line[1]
            GT_comm = line[2]
            SVID_comm = line[3]
            try:
                baseL = self.vcfbase[SVID_base] #[chrID, start, end, svlen, CHR2]
                commL = self.vcfcomm[SVID_comm]
            except:
                print("erro!!", SVID_base, SVID_comm)
            if (baseL[0] == commL[0]) and (baseL[-1] == commL[-1]):
                for k in  bc1Dict.keys():
                    kk = [int(i) for i in k.split('-')]
                    if (kk[0] < commL[3] <= kk[1]) or (kk[0] < baseL[3] <= kk[1]):
                        bc1Dict[k].append(int(commL[1]) - int(baseL[1]))
                        bc2Dict[k].append(int(commL[2]) - int(baseL[2]))
                        lsvDict[k].append(int(commL[3]) - int(baseL[3]))
                        continue
        bc1s = [item for subl in list(bc1Dict.values()) for item in subl]
        bc2s = [item for subl in list(bc2Dict.values()) for item in subl]
        lsvs = [item for subl in list(lsvDict.values()) for item in subl]
        return bc1s,bc2s,lsvs,bc1Dict,bc2Dict,lsvDict

    def F1_recall_precision(self):
        F1_recall_precisionGTDict ={"F1":{},"recall":{},"precision":{},"all":[]}
        F1_recall_precisionNGTDict ={"F1":{},"recall":{},"precision":{},"all":[]}
        for k in self.TP_no_GT_comm.keys():
            baseGT = len(set(self.TP_GT_base[k]))
            commGT = len(set(self.TP_GT_comm[k]))
            baseNGT = len(set(self.TP_no_GT_base[k]))
            commNGT = len(set(self.TP_no_GT_comm[k]))
            try:
                precision_GT = self.oneselcet(int((commGT)/float(self.lencomm[k])))
                recall_GT =  self.oneselcet(self.oneselcet(int((baseGT)/float(self.lenbase[k]))))
                F1_GT = 2*precision_GT*recall_GT/(precision_GT+recall_GT)
            except ZeroDivisionError:
                precision_GT,recall_GT,F1_GT = 0,0,0
            F1_recall_precisionGTDict["F1"][k] = F1_GT
            F1_recall_precisionGTDict["recall"][k] = recall_GT
            F1_recall_precisionGTDict["precision"][k] = precision_GT

            try:
                precision_no_GT = self.oneselcet(int((commNGT)/float(self.lencomm[k])))
                recall_no_GT =  self.oneselcet(int((baseNGT)/float(self.lenbase[k])))
                F1_no_GT = 2*precision_no_GT*recall_no_GT/(precision_no_GT+recall_no_GT)
            except ZeroDivisionError:
                precision_no_GT,recall_no_GT,F1_no_GT = 0,0,0
            F1_recall_precisionNGTDict["F1"][k] = F1_no_GT
            F1_recall_precisionNGTDict["recall"][k] = recall_no_GT
            F1_recall_precisionNGTDict["precision"][k] = precision_no_GT
        
        try:
            precisionAllGT = self.oneselcet(int(sum([len(set(i)) for i in self.TP_GT_comm.values()]))/float(sum(self.lencomm.values())))
            precisionAllNGT = self.oneselcet(int(sum([len(set(i)) for i in self.TP_no_GT_comm.values()]))/float(sum(self.lencomm.values())))
            recallAllGT = self.oneselcet(int(sum([len(set(i)) for i in self.TP_GT_base.values()]))/float(sum(self.lenbase.values())))
            recallAllNGT = self.oneselcet(int(sum([len(set(i)) for i in self.TP_no_GT_base.values()]))/float(sum(self.lenbase.values())))
            F1AllGT = 2*precisionAllGT*recallAllGT/(precisionAllGT+recallAllGT)                              
            F1AllNGT = 2*precisionAllNGT*recallAllNGT/(precisionAllNGT+recallAllNGT)
        except ZeroDivisionError:
            precisionAllGT,precisionAllNGT,recallAllGT,recallAllNGT,F1AllGT,F1AllNGT = 0,0,0,0,0,0
        F1_recall_precisionGTDict["all"] = [precisionAllGT,recallAllGT,F1AllGT]
        F1_recall_precisionNGTDict["all"] = [precisionAllNGT,recallAllNGT,F1AllNGT]
        self.writeReslut(F1_recall_precisionGTDict,F1_recall_precisionNGTDict)


    def base_comm(self):
        for chrx, linesbase in self.vcfbase.items():
            for  linebase in linesbase:
                if chrx in self.vcfcomm:
                    for linecomm in self.vcfcomm[chrx]:
                        if (self.svType in ['DEL', 'INS', 'DUP', 'INV']):
                            Result, bc1, bc2, lsv, rangeId = self.del_ins_inv_dup(linebase,linecomm)
                            if Result:
                                self.TP_no_GT_comm[rangeId].append(linecomm[-2])
                                self.TP_no_GT_base[rangeId].append(linebase[-2])
                                self.bc1_no_GT[rangeId].append(bc1)
                                self.bc2_no_GT[rangeId].append(bc2)
                                self.lsv_no_GT[rangeId].append(lsv)
                                if linebase[5] == linecomm[5]:
                                    self.TP_GT_comm[rangeId].append(linecomm[-2])
                                    self.TP_GT_base[rangeId].append(linebase[-2])
                                    self.bc1_GT[rangeId].append(bc1)
                                    self.bc2_GT[rangeId].append(bc2)
                                    self.lsv_GT[rangeId].append(lsv)
                        elif (self.svType in ['TRA', 'BND']):
                            Result, bc1, bc2 = self.tra_bnd(linebase,linecomm)
                            rangeId = "50-100"
                            if Result:
                                self.TP_no_GT_comm[rangeId].append(linecomm[-2])
                                self.TP_no_GT_base[rangeId].append(linebase[-2])
                                self.bc1_no_GT[rangeId].append(bc1)
                                self.bc2_no_GT[rangeId].append(bc2)
                                if linebase[5] == linecomm[5]:
                                    self.TP_GT_comm[rangeId].append(linecomm[-2])
                                    self.TP_GT_base[rangeId].append(linebase[-2])
                                    self.bc1_GT[rangeId].append(bc1)
                                    self.bc2_GT[rangeId].append(bc2)
        self.F1_recall_precision()

    def run(self):
        self.base_comm()
        #bc1s,bc2s,lsvs,bc1Dict,bc2Dict,lsvDict = self.lsv_bc(gs)#计算偏差
        #precisionDict, recallDict, F1Dict = self.lenDist(lenDictCommTP, lenDictCommAll, lenDictBaseTP, lenDictBaseAll)#计算不同长度SV的F1 precision recall
        #return bc1s,bc2s,lsvs,lenDictCommTP,lenDictCommAll,lenDictBaseTP,lenDictBaseAll,precisionDict,recallDict,F1Dict,bc1Dict,bc2Dict,lsvDict

if __name__ == '__main__':
    evalsv = EvaluationSV(VcfBaseFile=sys.argv[1], VcfCommFile=sys.argv[2], svType=sys.argv[3], outFile="tet", RE_nu= 1, nux = 1)
    evalsv.base_comm()
