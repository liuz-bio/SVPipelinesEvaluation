import os
import sys

def lsv_F1(fileName):
    lsvDict = {}
    for line in open(fileName,'r'):
        line = line.strip().split(',')
        if "ID" == line[0]:
            heads = line[1:]
            continue
        lsvDict[line[0]]=dict(zip(heads,line[1:]))
    print(fileName)
    return lsvDict
    

samples = ["HG002"]
#resultDict = {plat:{pipline:{svtype:{depth:{RE:{}}}}}}
resultDict = {}
for line in os.popen("cat vcf.info|cut -f 1,2,3,4,5|grep -v TRA").read().strip().split("\n"):
    line = line.strip().split('\t')
    #print(line)
    plat = line[0]
    call = line[4]
    maps = line[3]
    svtype = line[2]
    depth = line[1]
    pipline = maps+'-'+call
    for sa in [plat]:
        sams = sa.lower()+'_'+svtype.lower()+'_'+depth+'x' #EvalOutFile/Nanopore/minimap2/pbsv/DEL/15/11/nanopore_del_15x/summary.txt
        #dataDict = dict(zip(range(2,21),["0" for i in range(2,21)]))
        dataDict = dict(zip(range(2,21),[{}]*19))
        for idx in range(2,21):
            #EvalOutFile/Nanopore/lordfast/cutesv/DEL/15/2/nanopore_del_15x/summary.txt
            path = "/".join(["EvalOutFile",plat,maps,call,svtype,depth,str(idx),sams,"lsv_F1.txt"])
            if os.path.exists(path) :
                lsvDict = lsv_F1(path)
                for ka,va in lsvDict.items():
                    for kb,vb in va.items():
                        if ka in dataDict[idx]:
                            if kb in dataDict[idx][ka]:
                                dataDict[idx][ka][kb].append(vb)
                            else:
                                dataDict[idx][ka][kb] = [vb]
                        else:
                            dataDict[idx][ka]={kb:[vb]}
                          
        #for ID in ["F1","recall","precision"]:
        for re,v1 in dataDict.items():
            for k2,v2 in v1.items():
                for k3,v3 in v2.items():
                    print('\t'.join([plat,pipline,svtype,depth,sa,k3,k2,','.join(v3)]))
#        print('\t'.join([plat,pipline,svtype,depth,sa,"F1",','.join([F1Dict[i] for i in range(2,21)])]))
#        print('\t'.join([plat,pipline,svtype,depth,sa,"recall",','.join([RecallDict[i] for i in range(2,21)])]))
#        print('\t'.join([plat,pipline,svtype,depth,sa,"precision",','.join([precisionDict[i] for i in range(2,21)])]))
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
