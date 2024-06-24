cat ../real_merge_vcf.bed|grep -v hg38.vcf|grep -vP "\tDEL\tHG002\t"|grep -vP "\tINS\tHG002\t"|grep -vP "\tDEL\tCHM13\t"|grep -vP "\tINS\tCHM13\t"|sed 's/\t/:/g'|while read id;do(i=(${id//:/ });cat ${i[2]}|grep -v \#|cut -f 1,2,8|sed 's/SUPP.*;END=//g'|cut -d\; -f1|awk -v OFS='\t' '{print $1,$2-500,$3+500}'|grep -v \_ |vcf-sort|python work.py > ${i[0]}.${i[1]}.bed);done

cat ../real_merge_vcf.bed|grep hg38.vcf|grep -vP "\tDEL\tHG002\t"|grep -vP "\tINS\tHG002\t"|grep -vP "\tDEL\tCHM13\t"|grep -vP "\tINS\tCHM13\t"|sed 's/\t/:/g'|while read id;do(i=(${id//:/ });cat ${i[2]}|grep -v \#|cut -f1,2,8|sed "s/PRECISE;END=//g"|cut -d\; -f1|awk -v OFS='\t' '{print $1,$2-500,$3+500}'|grep -v \_ |vcf-sort|python work.py > ${i[0]}.${i[1]}.bed);done

cat /public3/SVDataset_lz/backup/2023_11_14/SV_static/SV_new/test/CHM13/chm13v1.0_with38Y_to_GRCh38.dip.bed|grep -v \_ >DEL.CHM13.bed
cat /public3/SVDataset_lz/backup/2023_11_14/SV_static/SV_new/test/CHM13/chm13v1.0_with38Y_to_GRCh38.dip.bed|grep -v \_ >INS.CHM13.bed

cat /public3/SVDataset_lz/backup/2023_11_14/SV_static/SV_new/test/HG002/HG002_SVs_Tier1_v0.6.bed|grep -v 'chr[0-9,Un,XY]*_' >HG002_SVs_HG19_new_v0.6.bed
/public3/SVDataset_lz/backup/2023_11_14/SV_static/Data/bed/liftOver HG002_SVs_HG19_new_v0.6.bed /public3/SVDataset_lz/backup/2023_11_14/SV_static/Data/bed/hg19ToHg38.over.chain tmp_DEL.HG002.bed HG002.hg19.bed.unmap
cat tmp_DEL.HG002.bed|grep -v 'chr[0-9,Un,XY]*_' >DEL.HG002.bed
cat tmp_DEL.HG002.bed|grep -v 'chr[0-9,Un,XY]*_' >INS.HG002.bed
rm HG002.hg19.bed.unmap HG002_SVs_HG19_new_v0.6.bed tmp_DEL.HG002.bed


ls INS.*.bed|cut -d\. -f2|while read id;do(cat INS.$id.bed DUP.$id.bed|vcf-sort >DUP_INS.$id.bed);done

ls /public3/SVDataset_lz/backup/2023_11_14/SV_static/SV_new2/bed/bed_test/*.bed|while read id;do(a=${id##*\/};b=${a%.*};echo -e $b"\t"$id);done >bed.tab
