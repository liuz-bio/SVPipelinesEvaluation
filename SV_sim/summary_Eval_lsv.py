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

def readsummary(summaryFile):
    axa = {'ID':0,'baseT':1,'commT':2,'baseAll':3,'commAll':4,'precision':5,'recall':6,'F1':7}
    #axa = {0:'ID', 1:'baseT', 2:'commT',3:'baseAll',4:'commAll', 5:'precision', 6:'recall', 7:'F1'}
    summaryDict = {}
    for i in open(summaryFile, 'r'):
        i = i.strip().split(',')
        if 'ID' in i:
            continue
        for ii in ['baseT', 'commT', 'baseAll', 'commAll', 'precision','recall','F1']:
            try:
                summaryDict[ii].append(i[axa[ii]])
            except:
                summaryDict[ii] = [i[axa[ii]]]
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
                        path = os.path.join(base,pi,mi,ci,si,di,re,sampleDir,'lsv_F1.txt')
                        #sampleDir = sa + '.' + di + 'x'#nanopore_tra_5x
                        #path = os.path.join(base,pi,mi,ci,si,di,re,sampleDir,'lsv_F1.txt')
                        if os.path.exists(path):
                            summaryDict = readsummary(path)
                            if len(summaryDict) == 0:
                                continue
                            for k1,v1 in summaryDict.items():
                                print('Sim',pi,mi,ci,si,di+'x',re,k1,','.join(v1))
