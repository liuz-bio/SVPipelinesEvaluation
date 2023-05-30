import sys
tabFiles = sys.argv[1]
for line in open(tabFiles,'r'):
    lines = line.split('/')
    plat = lines[6]
    depth = lines[7].split('_')[-1][:-1]
    maptool = lines[8]
    calltool = lines[10].lower()
    IDx = lines[7]+'.vcf'
    for svtype in ['TRA','DEL','INS','DUP','INV']:
        print('\t'.join([plat,depth,svtype,maptool,calltool,line.strip(),IDx]))
#/home/lz/Data_sequence/2021_4_2/HG002/Nanopore/HG002.15x/winnowmap/vcf/cuteSV/HG002.15x.vcf
