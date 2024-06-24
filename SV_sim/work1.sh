cp ~/Data_sequence/2020_5_11/SV_work/three_static/*.py .

cp ~/Data_sequence/2020_5_11/SV_work/three_static/Head.info .

mkdir SimOutPath/ CallOutPath/ bed/ EvalOutFile/

ll /home/lz/Data_sequence/2023_11_14/working/Pipelines/*/*/*/vcf/svision/*/*.svision.s2.graph.vcf|perl -pe 's/^.*?\//\//g'|while read id;do(i=${id%/*};cp $id $i.vcf);done
ll /home/lz/Data_sequence/2023_11_14/working/Pipelines/*/*/*/vcf/svim/*/*/variants.vcf|perl -pe 's/^.*?\//\//g'|while read id;do(i=${id%/*};ii=${i%/*};cp $id $ii.vcf);done
ll /home/lz/Data_sequence/2023_11_14/working/Pipelines/*/*/*/vcf/nanovar/*/*/*nanovar.sort.nanovar.pass.vcf|perl -pe 's/^.*?\//\//g'|while read id;do(i=${id%/*};ii=${i%/*};cp $id $ii.vcf);done

ll /home/lz/Data_sequence/2023_11_14/working/Pipelines/*/*/*/vcf/picky/*.INS.xls.vcf|perl -pe 's/^.*?\//\//g'|awk -v OFS='/' -v FS='/' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13" "$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12" "$9".vcf"}'|while read id;do(a=(${id// / });echo "/usr/bin/cat <(grep \# "${a[0]}") <(/usr/bin/cat "${a[1]}"/*xls.vcf|grep -v \#)|awk '\$2>0'|vcf-sort >"${a[1]}"/"${a[2]} );done >a.sh

ll /home/lz/Data_sequence/2023_11_14/working/Pipelines/*/*/*/vcf/debreak/*/debreak.vcf|perl -pe "s/^.*?\//\//g"|awk -v FS='/' -v OFS='/' '{print "cp "$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14"  "$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13".vcf"}' >a.sh

#ll /home/lz/Data_sequence/2023_11_14/working/Pipelines/*/*/*/vcf/*/*.vcf|grep -vE 'minialign|lordfast|profile'|perl -pe 's/^.*?\//\//g'|grep -v "lra/vcf/delly"|grep -v "lra/vcf/nanovar"|grep  -v "lra/vcf/NanoSV" >three.ll.tab

ls /home/lz/Data_sequence/2023_11_14/working/Pipelines/*/*/*/vcf/*/*.vcf|grep -v "pbmm2/vcf/nanovar"|grep -v "lra/vcf/nanovar"|grep -v "lra/vcf/NanoSV"|grep -v profile|grep -v "lra/vcf/pbsv"|grep -v "lra/vcf/delly" >three.ll.tab

python creat_vcf.info.py three.ll.tab >three.all.vcf.info

#comm -23 <(ll /home/lz/Data_sequence/2023_11_14/working/Pipelines/*/*/*/vcf/*/work.sh|grep -v "pbmm2/vcf/nanovar"|grep -v "lra/vcf/nanovar"|grep -v "lra/vcf/NanoSV"|grep -v profile|grep -v "lra/vcf/pbsv"|cut -d\/ -f8-12|sort) <(ll /home/lz/Data_sequence/2023_11_14/working/Pipelines/*/*/*/vcf/*/*.vcf|grep -v "pbmm2/vcf/nanovar"|grep -v "lra/vcf/nanovar"|grep -v "lra/vcf/NanoSV"|grep -v profile|grep -v "lra/vcf/pbsv"|cut -d\/ -f8-12|sort)

#comm -23 <(ll /home/lz/Data_sequence/2023_11_14/working/Pipelines/*/*/*/vcf/*/work.sh|grep -v "pbmm2/vcf/nanovar"|grep -v "lra/vcf/nanovar"|grep -v "lra/vcf/NanoSV"|grep -v profile|grep -v "lra/vcf/pbsv"|grep -v "lra/vcf/delly"|cut -d\/ -f8-12|sort) <(ll /home/lz/Data_sequence/2023_11_14/working/Pipelines/*/*/*/vcf/*/*.vcf|grep -v "pbmm2/vcf/nanovar"|grep -v "lra/vcf/nanovar"|grep -v "lra/vcf/NanoSV"|grep -v profile|grep -v "lra/vcf/pbsv"|grep -v "lra/vcf/delly"|cut -d\/ -f8-12|sort)
