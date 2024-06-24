import os
import sys

samples = ['CHM13','HG00096','HG002','HG003','HG004','HG005','HG00512','HG006','HG007','NA12878']
#resultDict = {plat:{pipline:{svtype:{depth:{RE:{}}}}}}
resultDict = {}
for line in os.popen("cat three.all.vcf.info three.all.vcf.info.DUP_INS|cut -f 1,2,3,4,5|grep -v TRA|sort -u").read().strip().split("\n"):
    line = line.strip().split('\t')
    #print(line) Nanopore        HG003.10        DEL     minimap2        cutesv
    plat = line[0]
    call = line[4]
    maps = line[3]
    svtype = line[2]
    depth = line[1]
    pipline = maps+'-'+call
    for sa in samples:
        #sams = sa+'.'+depth+'x' EvalOutFile/Nanopore/minimap2/cutesv/DEL/HG004.25/2/HG004.25x/
        #sams = depth+'x'
        sams = sa+'.'+depth
        F1Dict = dict(zip(range(2,21),["0" for i in range(2,21)]))
        RecallDict = F1Dict.copy()
        precisionDict = F1Dict.copy()
        path = "/".join(["EvalOutFile",plat,maps,call,svtype,depth,'2',sams,"lsv_F1.txt"])
        #ID,baseT,commT,baseAll,commAll,precision,recall,F1
        #50-100,2460,2342,2949,3260,0.7184049079754601,0.8341810783316378,0.7719762848550846
        #print(path)
        if (plat == 'ONT') and  (sa not in ['CHM13','HG00096','HG002','HG003','HG004','HG00512','NA12878']):
            continue
        if (plat  == 'CCS') and  (sa not in ['CHM13','HG00096','HG002','HG003','HG004','HG005','HG00512','HG006','HG007','NA12878']):
            continue
        if (plat == 'CLR') and  (sa not in ['CHM13','HG002','HG003','HG004','HG005','HG00512','HG006','HG007']):
            continue
        if (maps == 'lra') and (call.lower() not in ['cutesv','cutesv2','debreak','delly','sniffles2','svim','svision']):
            continue
        if (maps == 'pbmm2') and (call.lower() not in ['cutesv','cutesv2','debreak','delly','nanosv','pbsv','picky','sniffles','sniffles2','svim','svision']):
            continue

        if os.path.exists(path) :
            lsvName = []
            dataDict = {}
            for lx in open(path,'r'):
                if 'baseT' in lx:
                    keys = lx.strip().split(',')
                    keys_len = len(keys)
                    continue
                lx = lx.strip().split(',')
                if lx[0] not in dataDict:
                    lsvName.append(lx[0])
                    dataDict[lx[0]] = {}
                    for i in range(1,keys_len):
                        dataDict[lx[0]][keys[i]] = lx[i]

            for k1 in keys[1:]:
                line = []
                for k2 in lsvName:
                    line.append(dataDict[k2][k1])
                print('\t'.join([plat,pipline,svtype,depth,sa,k1,','.join(line)]))
                
"""

    if plat in resultDict:
        if pipline in resultDict[plat]:
            if svtype in resultDict[plat][pipline]:
                if depth in resultDict[plat][pipline][svtype]:
                    resultDict[plat][pipline][svtype][depth] = {"F1":F1Dict, "recall":RecallDict, "precision":precisionDict}
                    print('存在重复')
                else:
                    resultDict[plat][pipline][svtype][depth] = {"F1":F1Dict, "recall":RecallDict, "precision":precisionDict}
            else:
                resultDict[plat][pipline][svtype] = {depth:{"F1":F1Dict, "recall":RecallDict, "precision":precisionDict}}
        else:
            resultDict[plat][pipline] = {svtype:{depth:{"F1":F1Dict, "recall":RecallDict, "precision":precisionDict}}}
    else:
        resultDict[plat] = {pipline:{svtype:{depth:{"F1":F1Dict, "recall":RecallDict, "precision":precisionDict}}}}


for plat,platV in resultDict.items():
    for pipline,piplineV in platV.items():
        for svtype,svtypeV in piplineV.items():
            for depth,depthV in svtypeV.items():
                for ID in ["F1","recall","precision"]:
                    print('\t'.join([plat,pipline,svtype,depth,ID,','.join([depthV[ID][i] for i in range(2,21)])]))
"""         
#Nanopore        10      TRA     minialign       cutesv
"""
precision_GT,0.3333333333333333
recall_GT,0.3333333333333333
F1_GT,0.3333333333333333
precision_no_GT,0.6666666666666666
recall_no_GT,0.5
F1_no_GT,0.5714285714285715
TP_GT_comm,1
TP_GT_base,1
TP_no_GT_comm,2
TP_no_GT_base,2
commSize,3
baseSize,4
"""
