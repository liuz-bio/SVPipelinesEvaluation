cat ../../../vcf.info |grep HG003|cut -f 6|sort -u >HG003.vcf.tab

cat ../../../vcf.info |grep HG004|cut -f 6|sort -u >HG004.vcf.tab

cat ../../../vcf.info |grep HG003|cut -f 6|sort -u|grep "15x" >HG003.vcf.15x.tab

cat ../../../vcf.info |grep HG004|cut -f 6|sort -u|grep "15x" >HG004.vcf.15x.tab

python work.merge.py HG003.vcf.15x.tab HG003.vcf.15x ;python work.merge.py HG004.vcf.15x.tab HG004.vcf.15x
