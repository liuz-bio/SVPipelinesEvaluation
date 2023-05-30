import pysam
import sys
import os
import json
import shutil
import gzip

class EvaluationSV():
    def __init__(self, VcfBaseFile=None, VcfCommFile=None, svType=None, outFile=None, sample=None, pwa_path=None, RE_nu= None, nux = None):
         #self.vcfbase = pysam.VariantFile(VcfBaseFile)
         #self.vcfcomm = pysam.VariantFile(VcfCommFile)
         self.location = 1000
         self.sample = sample
         self.nux = nux
         self.RE_nu = RE_nu
         self.svType = svType
         self.pwa_path = pwa_path
         self.outFile = outFile
         self.vcfhead = "/home/lz/Data_sequence/2020_5_11/SV_work/three_static/Head.info"
         self.SummaryDict = {}
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
         bedDicts = {"DEL.HG003":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/DEL.HG003.bed",
                     "DEL.HG004":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/DEL.HG004.bed",
                     "DEL.HG005":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/DEL.HG005.bed",
                     "DEL.HG006":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/DEL.HG006.bed",
                     "DEL.HG007":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/DEL.HG007.bed",
                     "DEL.NA12778":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/DEL.NA12778.bed",
                     "DUP.HG003":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/DUP.HG003.bed",
                     "DUP.HG004":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/DUP.HG004.bed",
                     "DUP.HG005":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/DUP.HG005.bed",
                     "DUP.HG006":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/DUP.HG006.bed",
                     "DUP.HG007":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/DUP.HG007.bed",
                     "DUP.NA12778":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/DUP.NA12778.bed",
                     "INS.HG003":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/INS.HG003.bed",
                     "INS.HG004":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/INS.HG004.bed",
                     "INS.HG005":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/INS.HG005.bed",
                     "INS.HG006":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/INS.HG006.bed",
                     "INS.HG007":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/INS.HG007.bed",
                     "INS.NA12778":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/INS.NA12778.bed",
                     "INV.HG003":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/INV.HG003.bed",
                     "INV.HG004":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/INV.HG004.bed",
                     "INV.HG005":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/INV.HG005.bed",
                     "INV.HG006":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/INV.HG006.bed",
                     "INV.HG007":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/INV.HG007.bed",
                     "INV.NA12778":"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/bed/INV.NA12778.bed"}
                     #'DUP':"/home/lz/Data_sequence/2021_4_2/SV_static/SV_HG002/bed/DUP.bed"}
         self.bed = bedDicts[svType+'.'+sample]
         self.lsv_base_GT_distr = {'50-100':0,'100-500':0,'500-1000':0,'1000-2500':0,'2500-99999999999':0}
         self.lsv_comm_GT_distr = {'50-100':0,'100-500':0,'500-1000':0,'1000-2500':0,'2500-99999999999':0}
         self.lsv_base_no_GT_distr = {'50-100':0,'100-500':0,'500-1000':0,'1000-2500':0,'2500-99999999999':0}
         self.lsv_comm_no_GT_distr = {'50-100':0,'100-500':0,'500-1000':0,'1000-2500':0,'2500-99999999999':0}
         self.tmpDict =  {'50-100':0,'100-500':0,'500-1000':0,'1000-2500':0,'2500-99999999999':0}
         self.VcfCommFile = VcfCommFile
         self.VcfBaseFile = VcfBaseFile
         self.vcfbase = self.readvcf(pysam.VariantFile(VcfBaseFile))
         self.vcfcomm = self.readvcf(pysam.VariantFile(VcfCommFile))
         self.vcfcomm_T = {}
         self.vcfbase_T = {}
    
    def getBedFile(self, VcfCommFile, outBed):
        out = open(outBed,'w')
        for i in open(VcfCommFile):
            if "#" in i:
                continue
            else:
                i = i.strip().split('\t')
                chrId = i[0]
                start = int(i[1])-1000
                end = str(int([ii for ii in i[7].split(';') if "END=" in ii][0].split('=')[-1])+1000)
                print(chrId,start,end)
            if "MT" in chrId:
                continue
            else:
                if start < 0:
                    start = 0
                out.write(chrId+'\t'+str(start)+'\t'+end+'\n')
        out.close()

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
        for line in VcfFile.fetch():
            chrID = str(line.chrom)
            lineID = str(line.id)#.strip().split('\t')[2]
            start = int(line.pos) 
            end   = int(line.stop)

            if self.svType in  ['TRA','BND']:
                CHR2 = [i.split('=')[-1] for i in str(line).split('\t')[7].split(';') if 'CHR2=' in i][0]
                svlen = 0
            else:
                CHR2 = 'None'
                svlen = self.getsvLEN(line)

            vcfDict[lineID] = [chrID, start, end, svlen, CHR2]
            
        return vcfDict

    def writeReslut(self, precision_GT,recall_GT,F1_GT,precision_no_GT,recall_no_GT,F1_no_GT,
                         TP_GT_comm,TP_GT_base,TP_no_GT_comm,TP_no_GT_base, commSize, baseSize,
                         bc1s,bc2s,lsvs,lenDictCommTP,lenDictCommAll,lenDictBaseTP,lenDictBaseAll,
                         precisionDict,recallDict,F1Dict,bc1Dict,bc2Dict,lsvDict):
        os.makedirs(self.outFile, exist_ok=True)
        headlines = open(self.vcfhead,'r').read().strip().split('\n')
        with open('%s/summary.txt'%self.outFile,'w') as out:
            out.write('precision_GT,'+str(precision_GT)+'\n')
            out.write('recall_GT,'+str(recall_GT)+'\n')
            out.write('F1_GT,'+str(F1_GT)+'\n')
            out.write('precision_no_GT,'+str(precision_no_GT)+'\n')
            out.write('recall_no_GT,'+str(recall_no_GT)+'\n')
            out.write('F1_no_GT,'+str(F1_no_GT)+'\n')
            out.write('TP_GT_comm,'+str(TP_GT_comm)+'\n')
            out.write('TP_GT_base,'+str(TP_GT_base)+'\n')
            out.write('TP_no_GT_comm,'+str(TP_no_GT_comm)+'\n')
            out.write('TP_no_GT_base,'+str(TP_no_GT_base)+'\n')
            out.write('commSize,'+str(commSize)+'\n')
            out.write('baseSize,'+str(baseSize)+'\n')
        with open('%s/bias.txt'%self.outFile,'w') as out:
            out.write('bc1,'+','.join(list(map(str, bc1s)))+'\n')
            out.write('bc2,'+','.join(list(map(str, bc2s)))+'\n')
            if self.svType in ['DEL','INS','INV','DUP']:
                out.write('lsv,'+','.join(list(map(str, lsvs)))+'\n')

        with open('%s/lsv_F1.txt'%self.outFile,'w') as out:
            #out.write(','.join(['ID','base_GT','comm_GT','base_no_GT','comm_no_GT'])+'\n')
            out.write(','.join(['ID','baseT','commT','baseAll','commAll','precision','recall','F1'])+'\n')
            tmp = self.tmpDict.copy()
            for k in tmp.keys():
                out.write(','.join([k,str(lenDictBaseTP[k]),str(lenDictCommTP[k]),
                                      str(lenDictBaseAll[k]),str(lenDictCommAll[k]),
                                      str(precisionDict[k]),str(recallDict[k]),str(F1Dict[k])])+'\n')

        with open('%s/lsv_bias.txt'%self.outFile,'w') as out:
            tmp = self.tmpDict.copy()
            for k in tmp.keys():
                out.write(','.join([k,'bc1',','.join(list(map(str,bc1Dict[k])))])+'\n')
                out.write(','.join([k,'bc2',','.join(list(map(str,bc2Dict[k])))])+'\n')
                out.write(','.join([k,'lsv',','.join(list(map(str,lsvDict[k])))])+'\n')

                
    #计算正确SV和所有的长度分布
    def lsv_Dict(self,vcfDict):
        tmp_dict =  self.tmpDict.copy()
        for k1, v1 in vcfDict.items():
            svLen = abs(int(v1[3]))# SV length
            for k,v in tmp_dict.items():
                kk = k.split('-')
                if int(kk[0])< svLen <=int(kk[1]):
                    tmp_dict[k]+=1
                    continue
        return tmp_dict

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

    def sort_tabix(self, vcfFile, outFileGz):
        os.system("cat %s |vcf-sort|bgzip > %s;tabix %s"%(vcfFile,outFileGz,outFileGz))

    def lsv_bc(self,gs):
        bc1Dict = dict(zip(self.tmpDict.keys(),[[] for i in self.tmpDict.keys()]))
        bc2Dict = dict(zip(self.tmpDict.keys(),[[] for i in self.tmpDict.keys()]))
        lsvDict = dict(zip(self.tmpDict.keys(),[[] for i in self.tmpDict.keys()]))
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

    def truvari(self, baseFile, commFile, outFile):
        bc1s = []
        bc2s = []
        if os.path.exists(outFile):
            shutil.rmtree(outFile)
        try:
            os.system("/home/lz/miniconda3/bin/truvari bench -b %s -c %s -o %s -r 1000 -p 0.00 --pctsim 0 --passonly --sizemax 500000 --reference /home/lz/Data_sequence/2021_4_2/visortest/simulation_pipline/Nanopore/vcf/TS/SV/data/hs37d5.fa --includebed  %s "%(baseFile,commFile,outFile,self.bed))
        except BaseException as e:
            print('truvari bench',e)
        try:
            self.vcfcomm_T = self.readvcf(pysam.VariantFile("%s/tp-call.vcf"%outFile))
            self.vcfbase_T = self.readvcf(pysam.VariantFile("%s/tp-base.vcf"%outFile))
            lenDictCommTP = self.lsv_Dict(self.vcfcomm_T)
            lenDictBaseTP = self.lsv_Dict(self.vcfbase_T)
            lenDictCommAll = self.lsv_Dict(self.vcfcomm)
            lenDictBaseAll = self.lsv_Dict(self.vcfbase)
        except BaseException as e:
            print('ccccccccccccccccccccc',e)
        try:
            tabfile = "%s/a.%s.%s.tab"%(self.pwa_path,str(self.nux),str(self.RE_nu))
            with open(tabfile, 'w') as tabout:
                tabout.write("%s/tp-base.vcf\n%s/tp-call.vcf\n"%(self.outFile,self.outFile))
        except BaseException as e:
            print(e)
#        os.system("echo -e '%s/tp-base.vcf\n%s/tp-call.vcf' >a.%s.%s.tab"%(outFile,outFile,str(nux),str(self.RE_nu)))
        try:
            os.system("SURVIVOR merge %s 1000 1 1 1 0 30 %s.vcf"%(tabfile,tabfile))
            gs = os.popen("cat %s.vcf|grep -v \#|grep 'SUPP_VEC=11;'|cut -f 10-11|sed 's/\t/:/g'|cut -d\: -f 1,8,12,19"%(tabfile)).read().strip().split('\n')
        except BaseException as e:
            print('gsssss',len(gs),gs)
        
        for path in [self.baseFile,self.baseFile+'.tbi',self.commFile,self.commFile+'.tbi',tabfile,tabfile+'.vcf']:
            os.remove(path)
        bc1s,bc2s,lsvs,bc1Dict,bc2Dict,lsvDict = self.lsv_bc(gs)#计算偏差
        precisionDict, recallDict, F1Dict = self.lenDist(lenDictCommTP, lenDictCommAll, lenDictBaseTP, lenDictBaseAll)#计算不同长度SV的F1 precision recall
        return bc1s,bc2s,lsvs,lenDictCommTP,lenDictCommAll,lenDictBaseTP,lenDictBaseAll,precisionDict,recallDict,F1Dict,bc1Dict,bc2Dict,lsvDict
    
    def readSummary(self):
        file = open('%s/summary.txt'%self.outFile, 'r')
        js = file.read()
        dic = json.loads(js)
        file.close()
        return dic

    def base_comm(self):
        self.baseFile = "%s/tmp.base.%s.%s.vcf.gz"%(self.pwa_path,str(self.nux),str(self.RE_nu))
        self.commFile = "%s/tmp.comm.%s.%s.vcf.gz"%(self.pwa_path,str(self.nux),str(self.RE_nu))
        self.sort_tabix(self.VcfCommFile, self.commFile)
        self.sort_tabix(self.VcfBaseFile, self.baseFile)
        #self.bed = "tmp.bed.%s.%s.bed"%(str(self.nux),str(self.RE_nu))
        #try:
        #    self.getBedFile(self.VcfCommFile, self.bed)
        #except BaseException as e:
        #    print('dedeed',e)
        try:
            bc1s,bc2s,lsvs,lenDictCommTP,lenDictCommAll,lenDictBaseTP,lenDictBaseAll,precisionDict,recallDict,F1Dict,bc1Dict,bc2Dict,lsvDict = self.truvari(self.baseFile, self.commFile, self.outFile)
        except BaseException as e:
            print('aaaaa',e)
        
        self.SummaryDict = self.readSummary()
        precision_GT, recall_GT, F1_GT = self.SummaryDict["gt_precision"], self.SummaryDict["gt_recall"], self.SummaryDict["gt_f1"]
        precision_no_GT,recall_no_GT,F1_no_GT = self.SummaryDict["precision"], self.SummaryDict["recall"], self.SummaryDict["f1"]
        TP_GT_comm, TP_GT_base = self.SummaryDict['TP-call_TP-gt'], self.SummaryDict['TP-base_TP-gt']
        TP_no_GT_comm, TP_no_GT_base = self.SummaryDict['TP-call'], self.SummaryDict['TP-base']
        commSize, baseSize = self.SummaryDict['call cnt'], self.SummaryDict['base cnt']
        try:
            self.writeReslut(precision_GT,recall_GT,F1_GT,precision_no_GT,recall_no_GT,F1_no_GT,
                         TP_GT_comm,TP_GT_base,TP_no_GT_comm,TP_no_GT_base, commSize,
                         baseSize,bc1s,bc2s,lsvs,lenDictCommTP,lenDictCommAll,lenDictBaseTP,
                         lenDictBaseAll,precisionDict,recallDict,F1Dict,bc1Dict,bc2Dict,lsvDict)
        except BaseException as e:
            print('bbbbb',e)

    def run(self):
        self.base_comm()

if __name__ == '__main__':
    evalsv = EvaluationSV(VcfBaseFile=sys.argv[1], VcfCommFile=sys.argv[2], svType=sys.argv[3], outFile="tet", nux=1)
    evalsv.base_comm()
    
