import os
import sys

samples = ["HG002"]
#resultDict = {plat:{pipline:{svtype:{depth:{RE:{}}}}}}
resultDict = {}
for line in os.popen("cat three.all.vcf.info|cut -f 1,2,3,4,5|grep -v TRA").read().strip().split("\n"):
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
        sams = depth+'x'
        F1Dict = dict(zip(range(2,21),["0" for i in range(2,21)]))
        RecallDict = F1Dict.copy()
        precisionDict = F1Dict.copy()
        for idx in range(2,21):
            path = "/".join(["EvalOutFile",plat,maps,call,svtype,depth,str(idx),sams,"summary.txt"])
            #print(path)
            if os.path.exists(path) :
                for lx in open(path,'r'):
                    lx = lx.strip().split(',')
                    if lx[0] == "F1_no_GT":
                        F1Dict[idx] = lx[1]
                    elif lx[0] == "recall_no_GT":
                        RecallDict[idx] = lx[1]
                    elif lx[0] == "precision_no_GT":
                        precisionDict[idx] = lx[1]
        #for ID in ["F1","recall","precision"]:
        print('\t'.join([plat,pipline,svtype,depth,sa,"F1",','.join([F1Dict[i] for i in range(2,21)])]))
        print('\t'.join([plat,pipline,svtype,depth,sa,"recall",','.join([RecallDict[i] for i in range(2,21)])]))
        print('\t'.join([plat,pipline,svtype,depth,sa,"precision",','.join([precisionDict[i] for i in range(2,21)])]))
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
