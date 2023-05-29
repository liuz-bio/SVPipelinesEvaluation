import sys
tabFiles = sys.argv[1]
for line in open(tabFiles,'r'):
    lines = line.split('/')
    plat = lines[7]
    depth = lines[8].split('_')[-1][:-1]
    maptool = lines[9]
    calltool = lines[11].lower()
    IDx = lines[8]+'.vcf'
    svtype = lines[8].split('_')[1].upper()
    print('\t'.join([plat,depth,svtype,maptool,calltool,line.strip(),IDx]))
#/home/lz/Data_sequence/2021_4_2/visortest/simulation_pipline/Nanopore/nanopore_dup_25x/minialign/vcf/pbsv/nanopore_dup_25x.vcf
