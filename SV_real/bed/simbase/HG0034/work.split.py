import os
import sys
from collections import Counter

fileName = sys.argv[1]
group = sys.argv[2]
suppx = int(sys.argv[3])

Datas = {'DEL':[],'INS':[],'DUP':[],'INV':[]}

heads = []

nu_del,nu_ins,nu_dup,nu_inv = 0,0,0,0
tt = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
for lx in open(fileName,'r'):
    if '#' in lx:
        if tt in lx:
            heads.append(tt+'merge\n')
        else:
            heads.append(lx)
        continue
    lxs = lx.strip().split('\t')

    if '[' in lxs[4] or ']' in lxs[4]:
        continue

    infs = [i for i in lxs[7].split(';') if '=' in i]
    infDs = dict(zip([i.split('=')[0] for i in infs],[i.split('=')[1] for i in infs]))
    supp = int(infDs["SUPP"])
    svtype = str(infDs["SVTYPE"])
 
    geotyes = [i.split(':')[0] for i in lxs[9:] if i.split(':')[0]!='./.']
    b = dict(Counter(geotyes))
    c = dict(zip([ i/float(sum(b.values())) for i in b.values() ],[i for i in b.keys()]))
    print(c)
    if len(c)==0:
        continue

    if supp>=suppx and max(c.keys())>=0.7:
        try:
            if svtype == 'DEL':
                nu_del +=1
                ld = '\t'.join([lxs[0],lxs[1],'DEL.'+str(nu_del),'N','<DEL>',lxs[5],lxs[6],lxs[7],"GT:PSV:LN:DR",c[max(c.keys())]+":NaN:0:0,0"])
            elif svtype == 'INS':
                nu_ins +=1
                ld = '\t'.join([lxs[0],lxs[1],'INS.'+str(nu_ins),'N','<INS>',lxs[5],lxs[6],lxs[7],"GT:PSV:LN:DR",c[max(c.keys())]+":NaN:0:0,0"])
            elif svtype == 'DUP':
                nu_dup +=1
                ld = '\t'.join([lxs[0],lxs[1],'DUP.'+str(nu_dup),'N','<DUP>',lxs[5],lxs[6],lxs[7],"GT:PSV:LN:DR",c[max(c.keys())]+":NaN:0:0,0"])
            elif svtype == 'INV':
                nu_inv +=1
                ld = '\t'.join([lxs[0],lxs[1],'INV.'+str(nu_inv),'N','<INV>',lxs[5],lxs[6],lxs[7],"GT:PSV:LN:DR",c[max(c.keys())]+":NaN:0:0,0"])
            Datas[svtype].append(ld+'\n')
        except:
            continue

with open("%s.vcf"%group,'w') as out:
    for i in heads+Datas[group]:
        out.write(i)

