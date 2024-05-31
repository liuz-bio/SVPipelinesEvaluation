ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.TRA.25x/winnowmap/bam/Sim.TRA.25x.sort.bam /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.TRA.25x/winnowmap/bam/Sim.TRA.25x.cuteSV2.sort.bam
ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.TRA.25x/winnowmap/bam/Sim.TRA.25x.sort.bam.bai /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.TRA.25x/winnowmap/bam/Sim.TRA.25x.cuteSV2.sort.bam.bai
/home/lz/miniconda3/envs/cutesv/bin/cuteSV --genotype --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 -t 20 -s 2 -l 30 /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.TRA.25x/winnowmap/bam/Sim.TRA.25x.cuteSV2.sort.bam /public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/genome.fa /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.TRA.25x/winnowmap/vcf/cuteSV2/Sim.TRA.25x.vcf /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.TRA.25x/winnowmap/vcf/cuteSV2
