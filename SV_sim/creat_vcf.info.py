import sys
tabFiles = sys.argv[1]
for line in open(tabFiles,'r'):
    lines = line.split('/')
    plat = lines[9]
    depth = "25x"#lines[7].split('_')[-1][:-1]
    maptool = lines[11]
    calltool = lines[13].lower()
    IDx = lines[10]+'.25x.vcf'
    for svtype in ['TRA','DEL','INS','DUP','INV']:
        if maptool=="lra" and calltool=="delly":
            continue
        print('\t'.join([plat,depth,svtype,maptool,calltool,line.strip(),IDx]))
#/home/lz/Data_sequence/2021_4_2/HG002/Nanopore/HG002.15x/winnowmap/vcf/cuteSV/HG002.15x.vcf
#/public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/Pipelines/CCS/Sim.DEL.INS.DUP.INV.25x/lra/vcf/cuteSV2/Sim.DEL.INS.DUP.INV.25x.vcf
