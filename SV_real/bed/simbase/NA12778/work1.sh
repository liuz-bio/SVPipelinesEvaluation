#cat ../../../three.all.vcf.info |grep HG003|cut -f 6|sort -u >HG003.vcf.tab

#cat ../../../three.all.vcf.info |grep HG004|cut -f 6|sort -u >HG004.vcf.tab

python check.pbsv.head.py ../../../three.ll.tab

cat ../../../three.all.vcf.info |grep NA12778|cut -f 6|sort -u|grep "25x" >NA12778.vcf.25x.tab

ls ~/Data_sequence/2021_4_2/SV_static/SV_new/three_static/CallOutPath/*/*/*/*/NA12778.25/2/NA12778.25x.vcf |grep DEL >NA12778.vcf.25x/DEL.tab

echo -e "DEL\nINS\nDUP\nINV"|while read id;do(ls ~/Data_sequence/2021_4_2/SV_static/SV_new/three_static/CallOutPath/*/*/*/*/NA12778.25/2/NA12778.25x.vcf |grep $id >NA12778.vcf.25x/$id.tab);done
