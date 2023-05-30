echo -e 'HG005\nHG006\nHG007'|while read ida;do(echo -e 'DEL\nINS\nDUP\nINV'|while read idb;do(echo -e $idb"\t"$ida"\t"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/simbase/HG00567/$ida.vcf.25x/$idb.vcf);done);done

echo -e 'HG003\nHG004'|while read ida;do(echo -e 'DEL\nINS\nDUP\nINV'|while read idb;do(echo -e $idb"\t"$ida"\t"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/simbase/HG0034/$ida.vcf.25x/$idb.vcf);done);done

echo -e 'NA12778'|while read ida;do(echo -e 'DEL\nINS\nDUP\nINV'|while read idb;do(echo -e $idb"\t"$ida"\t"/home/lz/Data_sequence/2021_4_2/SV_static/SV_new/three_static/bed/simbase/NA12778/$ida.vcf.25x/$idb.vcf);done);done
