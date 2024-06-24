import os
import sys
from intervaltree import Interval, IntervalTree


def read_vcf(vcf_file):
    trees = {}
    with open(vcf_file, 'r') as f:
        for ix in f:
            if '#' in ix:
                continue
            i = ix.strip().split('\t')
            chrom = i[0]
            pos = int(i[1])
            tmp = [ii for ii in i[7].split(';') if '=' in ii]
            infos = dict(zip([ii.split('=')[0] for ii in tmp], [ii.split('=')[1] for ii in tmp]))
            end = int(infos['END'])+1
            infos['SEQ1'] = i[3]
            infos['SEQ2'] = i[4]
            infos['Line'] = ix
            trees.setdefault(chrom, IntervalTree()).addi(pos, end, infos)
    return trees

#sim DUP to INS
def compare(idx, trees_a, trees_b):
    num = 2919
    nu = 0
    nux = []
    for chrom, tree_a in trees_a.items():
        if chrom in trees_b:
            tree_b = trees_b[chrom]           
            for node in tree_a:
                start = node.begin
                end   = node.end
                svlen = int(node.data['SVLEN'])
                for i in tree_b.overlap(start, end):
                     if abs(svlen-int(i.data['SVLEN']))/min([svlen, int(i.data['SVLEN'])]) <0.7:
                         nux.append(node.data['Line'])
                         print(node.data['Line'].strip())
                         print(i.data['Line'].strip())
                         print("#############################################################################")
                         nu+=1
    print(idx, nu/2919, len(set(nux))/2919)

dup_file = "../SimOutPath/DUP/Sim/real_DUP.vcf"
ins_files = os.popen("ls /public3/SVDataset_lz/backup/2023_11_14/SV_static/SV_sim/CallOutPath/*/*/*/INS/25x/2/Sim.DEL.INS.DUP.INV.25x.25x.vcf").read().strip().split("\n")
dup_trees = read_vcf(dup_file)
for ins_file in ins_files:
    idx = '\t'.join(ins_file.split("/")[8:11])
    ins_trees = read_vcf(ins_file)
    compare(idx, dup_trees,ins_trees)
