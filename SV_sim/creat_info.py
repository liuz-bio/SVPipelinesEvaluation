import sys
import os

#lines = [ '/'+'/'.join(i.split('/')[1:]) for i in os.popen('cat ~/Data_sequence/2020_5_11/pbsv.raw.all.new.2.tab|grep -v lordfast|grep -v profile|grep -v "minialign/vcf/pbsv"|grep -v "minialign/vcf/picky"').read().strip().split('\n')]

lines = [ '/'+'/'.join(i.split('/')[1:]) for i in os.popen('cat ~/Data_sequence/2020_5_11/tmp.ll.tab |grep -v lordfast|grep -v profile|grep -v "minialign/vcf/pbsv"|grep -v "minialign/vcf/picky"').read().strip().split('\n')]
#Nanopore        10      INV     minimap2        cutesv  /home/lz/Data_sequence/2021_4_2/HG002/Nanopore/HG002.10x/minimap2/vcf/cuteSV/HG002.10x.vcf  HG002.10x.vcf

#/home/lz/Data_sequence/2020_5_11/NA12778/Nanopore/NA12778.10x/minialign/vcf/debreak/NA12778.10x.vcf
for lx in lines:
    lxs = lx.strip().split('/')
    plat = lxs[6]
    depth = lxs[-1].split('.')[1].replace('x','')
    maptool = lxs[8]
    caltool = lxs[10]
    vcfName = lxs[-1]
    for svtype in ['DEL','INS','DUP','INV']:
        print('\t'.join([plat,depth,svtype,maptool,caltool,lx,vcfName]))
       
    
