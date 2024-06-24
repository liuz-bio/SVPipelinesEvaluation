import os.path

from pbsvRead import pbsvRead
from cutesvRead import cutesvRead
from snifflesRead import snifflesRead
from sniffles2Read import sniffles2Read
from debreakRead import debreakRead
from dellyRead import dellyRead
from nanosvRead import nanosvRead
from nanovarRead import nanovarRead
from svimRead import svimRead
from pickyRead import pickyRead
from svisionRead import svisionRead
from createSimulateSV_Sample import createSimulateSV
from StaticTruvari import StaticTruvari
from overlapVcf import overlapVcf
#from EvaluationSV import EvaluationSV
#from EvaluationSV_other import EvaluationSV
from EvaluationSV_truvari import EvaluationSV
import multiprocessing
import traceback

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


def task_call_back(res):
    print(f'complete! {res}')

def task_err_call_back(err):
    print(f'出错啦~ error：{str(err)}')
    #print(err.__traceback__.tb_frame.f_globals['__file__'])
    #print(err.__traceback__.tb_lineno)


def funcSim(platform, depth, svType, vcfFile, outFile, ReferencePath, mapTool, svTool, RE_nu, nux):
    print(platform, depth, svType, vcfFile, outFile, ReferencePath, mapTool, svTool, RE_nu, nux)
    if svTool.lower() =='pbsv':
        SvToolVcf = pbsvRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                             outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    elif svTool.lower() == "cutesv":
        SvToolVcf = cutesvRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                               outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    elif svTool.lower() == "cutesv2":
        SvToolVcf = cutesvRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                               outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    elif svTool.lower() == "sniffles":
        SvToolVcf = snifflesRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    elif svTool.lower() == "sniffles2":
        SvToolVcf = sniffles2Read(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    elif svTool.lower() == "debreak":
        SvToolVcf = debreakRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                   outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    elif svTool.lower() == "delly":
        SvToolVcf = dellyRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                   outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    elif svTool.lower() == "nanosv":
        try:
            SvToolVcf = nanosvRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                   outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
        except Exception as e:
            print('ExceptionOOO: ',e)
    elif svTool.lower() == "nanovar":
        SvToolVcf = nanovarRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                   outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    elif svTool.lower() == "picky":
        SvToolVcf = pickyRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                   outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    elif svTool.lower() == "svim":
        SvToolVcf = svimRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                    outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    elif svTool == "svision":
        SvToolVcf = svisionRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                   outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,RE_nu=RE_nu,nux=nux)
    SvToolVcf.run()
    #try:
    #    SvToolVcf.run()
    #except Exception as e:
    #        print('ExceptionIII: ',vcfFile,e)


def readInfo(InfoName):
    multiprocessing.log_to_stderr()
    pool = multiprocessing.Pool(processes = 80)
    nux=0
    for ix in open(InfoName, 'r'):
        #print(ix)
        svTool, SvToolVcf = None, None
        ReferencePath = "/home/lz/Data_sequence/2023_11_14/Genome/hg38/genome.fa"
        ix = ix.strip().split('\t')
        if '#' in ix[0]:
            continue
        elif ix[4].lower() in ["pbsv", "nanosv","cutesv","cutesv2", "delly", "debreak", "svim", "sniffles", "sniffles2", "picky", "nanovar","svision"]:
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
            os.makedirs(os.path.join("CallOutPath", platform, mapTool, svTool, svType, depth,str(RE_nu)), exist_ok=True) 
            #print(platform, depth, svType, vcfFile, outFile, ReferencePath, mapTool, svTool, RE_nu, nux,)
            pool.apply_async(func=funcSim, args=(platform, depth, svType, vcfFile, outFile, ReferencePath, mapTool, svTool, RE_nu, nux,),callback=task_call_back, error_callback=task_err_call_back)
            #pool.apply_async(func=LogExceptions(funcSim), args=(platform, depth, svType, vcfFile, outFile, ReferencePath, mapTool, svTool, RE_nu, nux,))   #维持执行的进程总数为processes，当一个进程执行完毕后会添加新的进程进去
            nux+=1
    print("Mark~ Mark~ Mark~~~~~~~~~~~~~~~~~~~~~~")
    pool.close()
    pool.join()   #调用join之前，先调用close函数，否则会出错。执行完close后不会有新的进程加入到pool,join函数等待所有子进程结束
    print("Sub-process(es) done.")
"""
        if svTool =='pbsv':
            SvToolVcf = pbsvRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                               outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,nux=nux)
        elif svTool == "cutesv":
            SvToolVcf = cutesvRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,nux=nux)
        elif svTool == "sniffles":
            SvToolVcf = snifflesRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,nux=nux)
        elif svTool == "debreak":
            SvToolVcf = debreakRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                   outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,nux=nux)
        elif svTool == "delly":
            SvToolVcf = dellyRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                   outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,nux=nux)
        elif svTool == "nanosv":
            SvToolVcf = nanosvRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                   outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,nux=nux)
        elif svTool == "nanovar":
            SvToolVcf = nanovarRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                   outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,nux=nux)
        elif svTool == "picky":
            SvToolVcf = pickyRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                   outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,nux=nux)
        elif svTool == "svim":
            SvToolVcf = svimRead(platform=platform, depth=depth, svType=svType, vcfFile=vcfFile,
                                   outFile=outFile, ref=ReferencePath, mapTool=mapTool,svTool=svTool,nux=nux)
        SvToolVcf.run()
"""

def readSim(SimInfoName):
    Fil = []
    for line in open(SimInfoName,'r'):
        lines = line.strip().split('\t')
        svType = lines[0].upper()
        sample = lines[1]
        if svType == "TRA":
            svType = "BND"
        outFile  = os.path.join("SimOutPath", svType, sample, "real_%s.vcf"%svType.lower())
        if outFile in Fil:
            continue
        Fil.append(outFile)
        callSVFile = lines[2]
        ref = "/public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/genome.fa"
        sim = createSimulateSV(callSVFile=callSVFile,
                               svType=svType,
                               outFile=outFile,
                               ref=ref)
        sim.run()

def funcEval(simOutFile,callOutFile,svType,EvalOutFile,sample,pwa_path,RE_nu,nux):
    evalSV = EvaluationSV(VcfBaseFile=simOutFile, VcfCommFile=callOutFile, svType=svType, outFile=EvalOutFile,sample=sample,pwa_path=pwa_path,RE_nu=RE_nu,nux=nux)
    print(simOutFile,callOutFile,svType,EvalOutFile,sample,pwa_path,RE_nu,nux)
    try:
        evalSV.run()
    except Exception as e:
        print(e)

def evaluationSV(InfoName,pwa_path):
    multiprocessing.log_to_stderr()
    pool = multiprocessing.Pool(processes = 80)
    nux = 0
    for line in open(InfoName,'r'):
        lines = line.strip().split('\t')
        vcfFileName = lines[-1]
        sample = vcfFileName.split('.')[0]
        depth = lines[1]
        platform = lines[0]
        svType = lines[2].upper()
        if svType in ["TRA","BND"]:
            continue
            svType = "BND"
        mapTool = lines[3]
        svTool = lines[4]
        #callOutFile = os.path.join("CallOutPath",platform, mapTool, svTool, svType,depth,  vcfFileName)
        for RE_nu in range(2,21):
            callOutFile = os.path.join(pwa_path,"CallOutPath", platform, mapTool, svTool, svType, depth,str(RE_nu), vcfFileName)
            simOutFile =os.path.join(pwa_path,"SimOutPath", svType,sample, "real_%s.vcf"%svType.lower())
            EvalOutFile = os.path.join(pwa_path,"EvalOutFile",platform, mapTool, svTool, svType, depth,str(RE_nu), '.'.join(vcfFileName.split('.')[:-1]))
            os.makedirs(EvalOutFile, exist_ok=True)
            os.makedirs(os.path.join(pwa_path,"CallOutPath", platform, mapTool, svTool, svType, depth,str(RE_nu)), exist_ok=True)
            pool.apply_async(func=LogExceptions(funcEval), args=(simOutFile,callOutFile,svType,EvalOutFile,sample,pwa_path,str(RE_nu),'.'.join([platform,mapTool,svTool,svType,depth,str(nux)]),))
            nux+=1
        
    print("Mark~ Mark~ Mark~~~~~~~~~~~~~~~~~~~~~~")
    pool.close()
    pool.join()
    print("Sub-process(es) done.")

        #evalSV = EvaluationSV(VcfBaseFile=simOutFile, VcfCommFile=callOutFile, svType=svType, outFile=EvalOutFile)
        #evalSV.run()

if __name__ == '__main__':
    #readInfo("tt.info")#bed/vcf.info")
    readInfo("three.all.vcf.3.TRA.info")    #readInfo("three.all.vcf.infoi")
    #readInfo("three.nanosv.tab")
    #readInfo("three.no.pbsv.vcf.info")
    #readInfo("three.all.DUP.vcf.info")
    #readInfo("three.pbsv.vcf.info")
    #readSim("/public3/SVDataset_lz/backup/2023_11_14/SV_static/SV_new/bed/real_merge_vcf.bed")
    #evaluationSV("three.all.vcf.3.info","/public3/SVDataset_lz/backup/2023_11_14/SV_static/SV_new")
    #evaluationSV("/home/lz/Data_sequence/2020_5_11/SV_work/three_static/vcf.info")#bed/vcf.info")
    #readSim("bed/bed.info")
    #evaluationSV("/home/lz/Data_sequence/2021_4_2/SV_static/SV_new2/three_static/three.all.DUP.vcf.info","/home/lz/Data_sequence/2021_4_2/SV_static/SV_new2/three_static")
    #evaluationSV("/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/three.ngmlr.svim.vcf.info","/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static")
    #evaluationSV("/home/lz/Data_sequence/2020_5_11/SV_work/three_static/Pacbio.HG003.25x.pbmm2.sniffles.svim.info","/home/lz/Data_sequence/2020_5_11/SV_work/three_static")
    #static = StaticTruvari()
    #static.run()
    #over = overlapVcf()
    #over.run()

