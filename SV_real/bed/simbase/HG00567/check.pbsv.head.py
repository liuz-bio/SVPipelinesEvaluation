import os
import sys

for filename in open(sys.argv[1],'r'):
    filename = filename.strip()
    if '#' not in open(filename,'r').read():
        with open('tmp.vcf','w') as out:
            out.write('\n'.join([i.strip() for i in open('Head.pbsv.info','r')])+'\n')
            for lx in open(filename,'r'):
                lx = lx.strip()
                if '#' not in lx:
                    out.write(lx+'\n')
        print(filename)
        os.system("cp tmp.vcf %s"%(filename))
        os.remove("tmp.vcf")
