import os

lines = [ '/'+'/'.join(i.split('/')[1:]) for i in os.popen("cat ../../pbsv.raw.all.tab|grep -v profile|grep -v lordfast|grep -v 'minialign/vcf/pbsv'|grep -v 'minialign/vcf/picky'|grep ' 0 '").read().strip().split('\n')]

for lx in lines[1:]:
    path = '/'.join(lx.split('/')[:-1])
    with open("%s/work.qsub"%path,'w') as out:
        out.write("#PBS -q default\n")
        out.write("#PBS -V\n")
        out.write("#PBS -N delly\n")
        out.write("#PBS -l nodes=pan01:ppn=5\n")
        out.write("#PBS -l mem=30G\n")
        out.write("#PBS -o %s\n"%path)
        out.write("#PBS -j oe\n")
        out.write("source /home/lz/anaconda3/etc/profile.d/conda.sh\n")
        out.write('export PATH="/home/lz/anaconda3/bin:$PATH"\n')
        out.write("source activate base\n")
        out.write("bash %s/work.sh\n"%path)
    print("qsub %s/work.qsub"%path)
