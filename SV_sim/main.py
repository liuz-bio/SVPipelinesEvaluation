import os.path

from pbsvRead import pbsvRead
from cutesvRead import cutesvRead
from snifflesRead import snifflesRead
from debreakRead import debreakRead
from dellyRead import dellyRead
from nanosvRead import nanosvRead
from nanovarRead import nanovarRead
from svimRead import svimRead
from pickyRead import pickyRead
from svisionRead import svisionRead
from createSimulateSV import createSimulateSV
from StaticTruvari import StaticTruvari
from overlapVcf import overlapVcf
#from EvaluationSV import EvaluationSV
#from EvaluationSV_other import EvaluationSV
from EvaluationSV_truvari import EvaluationSV
import multiprocessing

def error(msg, *args):
    return multiprocessing.get_logger().error(msg, *args)


class LogExceptions(object):
    def __init__(self, callable):
        self.__callable = callable
        return

    def __call__(self, *args, **kwargs):
        try:
            result = self.__callable(*args, **kwargs)
        except Exception as e:
            # Here we add some debugging help. If multiprocessing's
            # debugging is on, it will arrange to log the traceback
            error(traceback.format_exc())
            # Re-raise the original exception so the Pool worker can
            # clean up
            raise
        # It was fine, give a normal answer
        return result
    pass

def funcSim(platform, depth, svType, vcfFile, outFile, ReferencePath, mapTool, svTool, RE_nu, nux):
    if svTool =='pbsv':
        SvToolVcf = pbsvRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                             outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    elif svTool == "cutesv":
        SvToolVcf = cutesvRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                               outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    elif svTool == "sniffles":
        SvToolVcf = snifflesRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    elif svTool == "debreak":
        SvToolVcf = debreakRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                   outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    elif svTool == "delly":
        SvToolVcf = dellyRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                   outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    elif svTool == "nanosv":
        SvToolVcf = nanosvRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                   outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    elif svTool == "nanovar":
        SvToolVcf = nanovarRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                   outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    elif svTool == "picky":
        SvToolVcf = pickyRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                   outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    elif svTool == "svim":
        SvToolVcf = svimRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                   outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    elif svTool == "svision":
        SvToolVcf = svisionRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                   outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    try:
        SvToolVcf.run()
    except Exception as e:
        print("Exception: ",e)

def readInfo(InfoName):
    multiprocessing.log_to_stderr()
    pool = multiprocessing.Pool(processes = 60)
    nux=0
    for ix in open(InfoName, 'r'):
        svTool, SvToolVcf = None, None
        ReferencePath = "/home/lz/Data_sequence/2021_4_2/visortest/simulation_pipline/SV/data/hs37d5.fa"
        ix = ix.strip().split('\t')
        if '#' in ix[0]:
            continue
        elif ix[4].lower() in ["pbsv","nanosv", "cutesv", "delly", "debreak", "svim", "sniffles", "picky", "nanovar","svision"]:
            #vcfFileName = ix[5].split('/')[-1]
            vcfFileName = ix[6]
            depth = ix[1]
            platform = ix[0]
            svType = ix[2].upper()
            if svType == "TRA":
                svType = "BND"
            mapTool = ix[3]
            vcfFile = ix[5]
            svTool = ix[4]
        else:
            continue

        #msg = "hello %d" %(i)
        for RE_nu in range(2,21):
            outFile = os.path.join("CallOutPath", platform, mapTool, svTool, svType, depth,str(RE_nu), vcfFileName)
            pool.apply_async(func=LogExceptions(funcSim), args=(platform, depth, svType, vcfFile, outFile, ReferencePath, mapTool, svTool, RE_nu, nux,))   
        nux+=1
    print("Mark~ Mark~ Mark~~~~~~~~~~~~~~~~~~~~~~")
    pool.close()
    pool.join()  
    print("Sub-process(es) done.")

def readSim(SimInfoName):
    Fil = []
    for line in open(SimInfoName,'r'):
        lines = line.strip().split('\t')
        svType = lines[0].upper()
        if svType == "TRA":
            svType = "BND"
        depth = lines[1]
        outFile  = os.path.join("SimOutPath", svType, "real_%s.vcf"%svType.lower())
        if outFile in Fil:
            continue
        Fil.append(outFile)
        SimulateBedFile = lines[3]
        ref = "/home/lz/Data_sequence/2021_4_2/visortest/simulation_pipline/SV/data/hs37d5.fa"
        print(SimulateBedFile)
        try:
            sim = createSimulateSV(SimulateBedFile=SimulateBedFile,
                               svType=svType,
                               outFile=outFile,
                               ref=ref, depth=depth)
        #try:
            sim.run()
        except BaseException as e:
            print(e)


def funcEval(simOutFile,callOutFile,svType,EvalOutFile,RE_nu,nux):
    try:
        evalSV = EvaluationSV(VcfBaseFile=simOutFile, VcfCommFile=callOutFile, svType=svType, outFile=EvalOutFile,RE_nu=RE_nu,nux=nux)
        evalSV.base_comm()
    except BaseException as e:
        print(e)

def evaluationSV(InfoName):
    multiprocessing.log_to_stderr()
    pool = multiprocessing.Pool(processes = 60)
    nux = 0
    for line in open(InfoName,'r'):
        lines = line.strip().split('\t')
        vcfFileName = lines[-1]
        sample = vcfFileName.split('.')[0]
        depth = lines[1]
        platform = lines[0]
        svType = lines[2].upper()
        if svType == "TRA":
            svType = "BND"
        mapTool = lines[3]
        svTool = lines[4]
        for RE_nu in range(2,21):
            callOutFile = os.path.join("CallOutPath", platform, mapTool, svTool, svType, depth,str(RE_nu), vcfFileName)
            simOutFile =os.path.join("SimOutPath", svType, "real_%s.vcf"%svType.lower())
            EvalOutFile = os.path.join("EvalOutFile",platform, mapTool, svTool, svType, depth,str(RE_nu),'.'.join(vcfFileName.split('.')[:-1]))
            os.makedirs(os.path.join("EvalOutFile",platform, mapTool, svTool, svType, depth,str(RE_nu)), exist_ok=True)
            pool.apply_async(func=LogExceptions(funcEval), args=(simOutFile,callOutFile,svType,EvalOutFile,str(RE_nu),str(nux),))
            nux+=1

    print("Mark~ Mark~ Mark~~~~~~~~~~~~~~~~~~~~~~")
    pool.close()
    pool.join()
    print("Sub-process(es) done.")





if __name__ == '__main__':
    readInfo("vcf.all.info")
    readSim("bed/Sim.info")
    evaluationSV("vcf.all.info")
