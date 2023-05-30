echo -e "DEL\nINS\nDUP\nINV"|while read id;do(ls ~/Data_sequence/2021_4_2/SV_static/SV_new/three_static/CallOutPath/*/*/*/*/HG003.25/2/HG003.25x.vcf |grep $id >./$id.tab);done

SURVIVOR merge INV.tab 500 9 1 1 0 30 INV.tab.vcf
python ../work.split.py INV.tab.vcf INV 9

SURVIVOR merge DUP.tab 500 9 1 1 0 30 DUP.tab.vcf
python ../work.split.py DUP.tab.vcf DUP 9

SURVIVOR merge INS.tab 500 30 1 1 0 30 INS.tab.vcf
python ../work.split.py INS.tab.vcf INS 30

SURVIVOR merge DEL.tab 500 30 1 1 0 30 DEL.tab.vcf
python ../work.split.py DEL.tab.vcf DEL 30
