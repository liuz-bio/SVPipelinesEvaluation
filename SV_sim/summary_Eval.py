import os
import sys

depth = ['5','10','15','25','30']
#svtype = ['BND','DEL','DUP','INS','INV']
#svtype = ['BND','DUP','INV']
svtype = ['DEL','INS','DUP','INV']
calltools = ['cutesv','debreak','delly','nanovar','pbsv','picky','sniffles','svim']
maptools = ['lordfast','minialign','minimap2','ngmlr','pbmm2','winnowmap']
plats = ['Nanopore','Pacbio']
typech = {'BND':'tra','DEL':'del','INS':'ins','DUP':'dup','INV':'inv'}
base = 'EvalOutFile'
baseSizeS = {'DUP':3328,'INV':44,'DEL':4871,'INS':7099}

def readsummary(summaryFile):
    summaryDict = {}
    for i in open(summaryFile, 'r'):
        i = i.strip().split(',')
        try:
            summaryDict[i[0]] = i[1]
        except:
            return summaryDict
    return summaryDict
    


for pi in plats:
    for mi in maptools:
        for ci in calltools:
            for di in depth:
                for re in [str(i) for i in  range(2,21)]:
                    commSizes = []
                    baseSizes = []
                    TP_no_GT_comms = []
                    TP_no_GT_bases = []
                    for si in svtype:
                        sampleDir = pi.lower() + '_' + typech[si] + '_' + di + 'x'#nanopore_tra_5x
                        #EvalOutFile/Nanopore/minimap2/cutesv/DEL/5/2/nanopore_del_5x/summary.txt
                        path = os.path.join(base,pi,mi,ci,si,di,re,sampleDir,'summary.txt')
                        baseSizes.append(baseSizeS[si])
                        if os.path.exists(path):
                            summaryDict = readsummary(path)
                            try:
                                commSizes.append(int(summaryDict['commSize']))
                            except:
                                commSizes.append(0)
                            try:
                                TP_no_GT_comms.append(int(summaryDict['TP_no_GT_comm']))
                            except:
                                TP_no_GT_comms.append(0)
                            try:
                                TP_no_GT_bases.append(int(summaryDict['TP_no_GT_base']))
                            except:
                                TP_no_GT_bases.append(0)
                            """
                            try:
                                commSizes.append(int(summaryDict['commSize']))
                                baseSizes.append(int(summaryDict['baseSize']))
                                TP_no_GT_comms.append(int(summaryDict['TP_no_GT_comm']))
                                TP_no_GT_bases.append(int(summaryDict['TP_no_GT_base']))
                            except:
                                continue
                            """
                    try:
                        precision = sum(TP_no_GT_comms)/float(sum(commSizes))
                    except ZeroDivisionError:
                        precision = 0

                    try: 
                        recall = sum(TP_no_GT_bases)/float(sum(baseSizes))
                    except ZeroDivisionError:
                        recall = 0
   
                    try:
                        F1 = 2*recall*precision/(precision+recall)
                    except ZeroDivisionError:
                        F1 = 0
                
                #print(pi,mi,ci,di,"baseSizes",baseSizes) 
                #print(pi,mi,ci,di,"commSizes",commSizes)
                #print(pi,mi,ci,di,"TP_no_GT_bases",TP_no_GT_bases)
                #print(pi,mi,ci,di,"TP_no_GT_comms",TP_no_GT_comms)
                    print(pi,mi,ci,di+'x',re,'precision',precision)
                    print(pi,mi,ci,di+'x',re,'recall', recall)
                    print(pi,mi,ci,di+'x',re,'F1', F1)
