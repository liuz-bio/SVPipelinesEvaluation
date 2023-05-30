#cat ../../../three.all.vcf.info |grep HG003|cut -f 6|sort -u >HG003.vcf.tab

#cat ../../../three.all.vcf.info |grep HG004|cut -f 6|sort -u >HG004.vcf.tab

python check.pbsv.head.py ../../../three.ll.tab

cat ../../../three.all.vcf.info |grep HG003|cut -f 6|sort -u|grep "25x" >HG003.vcf.25x.tab

cat ../../../three.all.vcf.info |grep HG004|cut -f 6|sort -u|grep "25x" >HG004.vcf.25x.tab

python work.merge.py HG003.vcf.25x.tab HG003.vcf.25x ;python work.merge.py HG004.vcf.25x.tab HG004.vcf.25x
