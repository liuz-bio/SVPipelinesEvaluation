#cat ../../../three.all.vcf.info |grep HG005|cut -f 6|sort -u >HG005.vcf.tab

#cat ../../../three.all.vcf.info |grep HG006|cut -f 6|sort -u >HG006.vcf.tab

python check.pbsv.head.py ../../../three.ll.tab

cat ../../../three.all.vcf.info |grep HG005|cut -f 6|sort -u|grep "25x" >HG005.vcf.25x.tab

cat ../../../three.all.vcf.info |grep HG006|cut -f 6|sort -u|grep "25x" >HG006.vcf.25x.tab

cat ../../../three.all.vcf.info |grep HG007|cut -f 6|sort -u|grep "25x" >HG007.vcf.25x.tab

python work.merge.py HG005.vcf.25x.tab HG005.vcf.25x 

python work.merge.py HG006.vcf.25x.tab HG006.vcf.25x

python work.merge.py HG007.vcf.25x.tab HG007.vcf.25x
