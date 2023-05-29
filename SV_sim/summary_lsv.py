import os
import sys

depth = ['5','10','15','25','30']
#svtype = ['BND','DEL','DUP','INS','INV']
#svtype = ['BND','DUP','INV']
svtype = ['DEL','INS','DUP','INS','INV']
calltools = ['cutesv','debreak','delly','nanovar','pbsv','picky','sniffles','svim']
maptools = ['lordfast','minialign','minimap2','ngmlr','pbmm2','winnowmap']
plats = ['Nanopore','Pacbio']
typech = {'BND':'tra','DEL':'del','INS':'ins','DUP':'dup','INV':'inv'}
base = 'EvalOutFile'

def readsummary(summaryFile):
    summaryDict = []
    for i in open(summaryFile, 'r').read().strip().split('\n')[1:]:
        i = i.strip().split(',')
        summaryDict.append(i[-1])
    return summaryDict
    


for pi in plats:
    for mi in maptools:
        for ci in calltools:
            for di in depth:
                for si in svtype:
                    sampleDir = pi.lower() + '_' + typech[si] + '_' + di + 'x'#nanopore_tra_5x
                    path = os.path.join(base,pi,mi,ci,si,di,sampleDir,'lsv.txt')
                    if os.path.exists(path):
                        summaryDict = readsummary(path)
                        #print(summaryDict)
                        tmp = ','.join(summaryDict)
                        print(pi,mi,ci,di,si,tmp)       
