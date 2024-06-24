import os
import sys
from concurrent.futures import ProcessPoolExecutor

def list_files(start_path):
    file_paths = []
    out = open("all.vcf.tab",'w')
    for root, dirs, files in os.walk(start_path):
        for file in files:
            file_path = os.path.join(root, file)
            if "Sim.TRA" in file_path:
                file_paths.append(file_path)
                out.write(file_path+'\n')
    out.close()
    return file_paths


def evalBND(path):
    #/public3/SVDataset_lz/backup/2023_11_14/SV_static/SV_sim/CallOutPath/CCS/lra/cutesv/BND/25x/2/Sim.TRA.25x.25x.vcf
    lx = path.strip().split("/")
    #CCS     25x     TRA     lra     cutesv2 /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/Pipelines/CCS/Sim.TRA.25x/lra/vcf/cuteSV2/Sim.TRA.25x.vcf       Sim.TRA.25x.25x.vcf
    plat = lx[8]
    depth = lx[12]
    svType = lx[11]
    maps = lx[9]
    call = lx[10]
    supp = lx[13]
    vcffile = path

    LASeR = "/public3/SVDataset_lz/backup/2023_11_14/Sim/bed/LASeR100.bed"
    bed = "/public3/SVDataset_lz/backup/2023_11_14/Sim/bed/sim_tra2_GRCh38.bed"
    summary_path = "/public3/SVDataset_lz/backup/2023_11_14/SV_static/SV_sim/EvalOutFile/"+plat+"/"+maps+"/"+call+"/"+svType+"/"+depth+"/"+supp+"/"+'.'.join(lx[-1].split('.')[:-1])
    os.makedirs(summary_path, exist_ok=True)
    summary = summary_path+"/summary.txt"
    
   #python eval_sim1.py BND /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/LASeR100.bed /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/sim_tra2_GRCh38.bed /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/Pipelines/CCS/Sim.TRA.25x/minimap2/vcf/cuteSV2/Sim.TRA.25x.vcf ./a.summary

    #/public3/SVDataset_lz/backup/2023_11_14/SV_static/SV_sim/EvalOutFile/CCS/lra/cutesv/DEL/25x/2/Sim.TRA.25x.25x/summary.txt
    cmd = ["python", "eval_sim1.py", "BND", LASeR, bed, vcffile, summary]
    print(' '.join(cmd))

def process_files_in_pool(lines, num_processes=4):
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        # 提交任务给进程池
        executor.map(evalBND, lines)

if __name__ == "__main__":
    start_path = "/public3/SVDataset_lz/backup/2023_11_14/SV_static/SV_sim/CallOutPath/"
    file_paths = list_files(start_path)
    #lines = os.popen("cat /public3/SVDataset_lz/backup/2023_11_14/SV_static/SV_sim/three.all.vcf.info|grep -P '\tTRA\t'|grep -v 'Sim.DEL'").read().strip().split("\n")

    process_files_in_pool(file_paths, num_processes=20)
