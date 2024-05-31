ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CCS/Sim.TRA.25x/lra/bam/Sim.TRA.25x.sort.bam /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CCS/Sim.TRA.25x/lra/bam/Sim.TRA.25x.cuteSV.sort.bam
ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CCS/Sim.TRA.25x/lra/bam/Sim.TRA.25x.sort.bam.bai /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CCS/Sim.TRA.25x/lra/bam/Sim.TRA.25x.cuteSV.sort.bam.bai
/home/lz/miniconda3/envs/cutesv1/bin/cuteSV --genotype --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 -t 20 -s 2 -l 30 /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CCS/Sim.TRA.25x/lra/bam/Sim.TRA.25x.cuteSV.sort.bam /public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/lra_CCS/genome.fa /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CCS/Sim.TRA.25x/lra/vcf/cuteSV/Sim.TRA.25x.vcf /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CCS/Sim.TRA.25x/lra/vcf/cuteSV
