import os
import sys

fileName = sys.argv[1]
group = sys.argv[2]

os.makedirs(group, exist_ok=True)

def getId(path):
    #/home/lz/Data_sequence/2020_5_11/HG0034/Nanopore/HG003.10x/minialign/vcf/cuteSV/HG003.10x.vcf
    ids = path.split('/')
    return '_'.join(ids[6:9]+[ids[10]]),ids[10]



def dealSVTYPE(calltool,svtype):
    if calltool == "svim" and svtype == "DUP:TANDEM":
        return DUP
    else:
        return svtype


for lx in open(fileName,'r'):
    lx =lx.strip()
    idx,calltool = getId(lx)
    data = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
    dataDict = {'head':[],'DEL':[],'INS':[],'INV':[],'DUP':[]}
      
    for line in open(lx,'r'):
        lines = line.strip().split('\t')
        #print(lines)
        if data in line:
            dataDict['head'].append(data+idx+'\n')
        elif '#' in line:
            dataDict['head'].append(line)
        elif '[' in lines[4] or ']' in lines[4]:
            continue
        elif lines[0] in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']:
            try:
                infs = [i for i in lines[7].split(';') if '=' in i]
                infD = dict(zip([i.split('=')[0] for i in infs],[i.split('=')[1] for i in infs]))
                svtype = dealSVTYPE(calltool,infD['SVTYPE'])
                dataDict[svtype].append(line)
            except:
                #print(infD)
                continue
        else:
            continue

    for k in ['DEL','INS','INV','DUP']:
        with open('%s/%s.%s.vcf'%(group,k,idx),'w') as out:
            for h in dataDict['head']:
                out.write(h)
            for v in dataDict[k]:
                out.write(v)

                
              
