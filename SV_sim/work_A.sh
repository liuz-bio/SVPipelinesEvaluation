ls /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/Pipelines/*/*/*/vcf/*/*.vcf|grep -v profile|grep -v "lra/vcf/nanovar" |grep -v "lra/vcf/NanoSV"|grep -v "lra/vcf/delly"|grep -v "pbmm2/vcf/nanovar" >three.ll.tab
python creat_vcf.info.py three.ll.tab >three.all.vcf.info
