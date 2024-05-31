ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/lra/bam/Sim.TRA.25x.sort.bam /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/lra/bam/Sim.TRA.25x.picky.sort.bam
ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/lra/bam/Sim.TRA.25x.sort.bam.bai /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/lra/bam/Sim.TRA.25x.picky.sort.bam.bai
samtools sort -n /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/lra/bam/Sim.TRA.25x.picky.sort.bam|samtools view -Sh |/home/lz/software/sv_callers/Picky_demo/src/picky.pl sam2align >/public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/lra/vcf/picky/Sim.TRA.25x.align
cat /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/lra/vcf/picky/Sim.TRA.25x.align|/home/lz/software/sv_callers/Picky_demo/src/picky.pl callSV --fastq /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/fastq/ONT_R10.4_VISOR_LASeR_Sim_TRA.fq.gz --genome /public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/lra_ONT/genome.fa  --oprefix /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/lra/vcf/picky/Sim.TRA.25x
ls /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/lra/vcf/picky/*.xls|grep -v profile.xls|while read id;do(/home/lz/software/sv_callers/Picky_demo/src/picky.pl xls2vcf --re 2 --xls $id >$id.vcf);done
cat <(grep \# /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/lra/vcf/picky/*.profile.INS.xls.vcf) <(cat /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/lra/vcf/picky/*xls.vcf|grep -v \#)|vcf-sort >/public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/lra/vcf/picky/Sim.TRA.25x.vcf
