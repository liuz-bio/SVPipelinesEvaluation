echo -e "DEL\nINS\nINV\nDUP"|while read id;do(mkdir -p SimOutPath/$id/Sim/;cat Head.info <(cat /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/sim_h*.bed|sed "s/\t/_/g"|sort|uniq -c|sort -k1rn|sed "s/    //g"|sed "s/  //g"|sed "s/tandem duplication/DUP/g"|sed "s/deletion/DEL/g"|sed   "s/insertion/INS/g"|sed "s/inversion/INV/g"|grep $id|sed "s/_/\t/g"|sed "s/ /\t/g"|cut -f1,2,3,4,5,6|awk -v OFS='\t' '{if ($5=="INS") {print $1,$2,$3,$4,$5,length($6)} else {print $1,$2,$3,$4,$5,$4-$3}}'|sed "s/DEL\t/DEL\t-/g"|sed "s/^1\t/0\/1\t/g"|sed "s/^2\t/1\/1\t/g"|awk -v OFS='\t' '{print $2,$3,"Sim-"$5"."NR,"N","<"$5">",".","PASS","PRECISE;END="$4";SVTYPE="$5";SVLEN="$6,"GT:DR:DV",$1":5:5"}')|vcf-sort >SimOutPath/$id/Sim/real_$id.vcf);done

ls /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/Pipelines/*/*/*/vcf/svim/*/*/*.vcf|while read id;do(paste <(echo "cp "$id) <(echo $id|cut -d\/ -f1-15|sed "s/$/.vcf/g"));done >a

ll /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/Pipelines/*/*/*/vcf/*/work.sh|grep -v "lra/vcf/nanovar" |grep -v "lra/vcf/NanoSV"|grep -v "lra/vcf/delly"|grep -v "pbmm2/vcf/nanovar"|wc -l

ls /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/Pipelines/*/*/*/vcf/svision/*/*.svision.s2.graph.vcf|while read id;do(paste <(echo "cp "$id) <(echo $id|cut -d\/ -f1-15|sed "s/$/.vcf/g"));done >a

ll /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/Pipelines/*/*/*/vcf/*/*.vcf|grep -v profile|grep -v "lra/vcf/nanovar" |grep -v "lra/vcf/NanoSV"|grep -v "lra/vcf/delly"|grep -v "pbmm2/vcf/nanovar"|wc -l 

cat <(cat three.all.vcf.info|grep -v Sim.DEL|grep -P "\tTRA\t"|sort -u) <(cat three.all.vcf.info|grep Sim.DEL|grep -v TRA|sort -u) >three.all.vcf.new.info
