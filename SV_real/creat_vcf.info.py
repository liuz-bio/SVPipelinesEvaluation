import sys
tabFiles = sys.argv[1]
for line in open(tabFiles,'r'):
    lines = line.split('/')
    plat = lines[7]
    depth = "25x"#lines[7].split('_')[-1][:-1]
    maptool = lines[9]
    calltool = lines[11].lower()
    IDx = lines[8]+'.25x.vcf'
    for svtype in ['TRA','DEL','INS','DUP','INV']:
        if maptool=="lra" and calltool=="delly":
            continue
        print('\t'.join([plat,depth,svtype,maptool,calltool,line.strip(),IDx]))
#/home/lz/Data_sequence/2021_4_2/HG002/Nanopore/HG002.15x/winnowmap/vcf/cuteSV/HG002.15x.vcf
