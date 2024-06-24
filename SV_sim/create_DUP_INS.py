import os
import sys

#CallOutPath/CCS/lra/cutesv/DUP/25x/2/HG003.25x.vcf

import os
from concurrent.futures import ProcessPoolExecutor

def list_files(start_path):
    file_paths = []
    out = open("all.vcf.tab",'w')
    for root, dirs, files in os.walk(start_path):
        for file in files:
            file_path = os.path.join(root, file)
            if "/INS/25x/" in file_path:
                file_paths.append(file_path)
                out.write(file_path+'\n')
    out.close()
    return file_paths


def modify_file(vcffile):
    head = "/public3/SVDataset_lz/backup/2023_11_14/SV_static/SV_new/Head.info"

    heads = open(head,'r').read().strip().split("\n")
    # 打开文件并读取内容
    with open(vcffile, 'r') as file:
        for lx in file.read().strip().split("\n"):
            if "#" not in lx:
                heads.append(lx)

    dup_vcffile = vcffile.replace("/INS/25x/","/DUP/25x/")

    with open(dup_vcffile, 'r') as file:
        for lx in file.read().strip().split("\n"):
            if "#" not in lx:
                heads.append(lx.replace('DUP','INS'))

    out_vcffile = vcffile.replace("/INS/25x/","/DUP_INS/25x/")
    os.makedirs('/'.join(out_vcffile.split('/')[:-1]),exist_ok=True)
    print(out_vcffile)
    #print('/'.join(out_vcffile.split('/')[:-1]))
    #print(dup_vcffile)
    # 写回文件
    with open(out_vcffile, 'w') as file:
        file.write('\n'.join(heads)+"\n")

def process_files_in_pool(file_paths, num_processes=4):
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        # 提交任务给进程池
        executor.map(modify_file, file_paths)

if __name__ == "__main__":
    folder_path = '/public3/SVDataset_lz/backup/2023_11_14/SV_static/SV_sim/CallOutPath'
    file_paths = list_files(folder_path)

    process_files_in_pool(file_paths, num_processes=50)
