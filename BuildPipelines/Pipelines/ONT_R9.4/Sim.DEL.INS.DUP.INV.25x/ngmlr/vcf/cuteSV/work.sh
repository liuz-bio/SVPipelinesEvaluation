ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.DEL.INS.DUP.INV.25x/ngmlr/bam/Sim.DEL.INS.DUP.INV.25x.sort.bam /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.DEL.INS.DUP.INV.25x/ngmlr/bam/Sim.DEL.INS.DUP.INV.25x.cuteSV.sort.bam
ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.DEL.INS.DUP.INV.25x/ngmlr/bam/Sim.DEL.INS.DUP.INV.25x.sort.bam.bai /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.DEL.INS.DUP.INV.25x/ngmlr/bam/Sim.DEL.INS.DUP.INV.25x.cuteSV.sort.bam.bai
/home/lz/miniconda3/envs/cutesv1/bin/cuteSV --genotype --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 -t 20 -s 2 -l 30 /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.DEL.INS.DUP.INV.25x/ngmlr/bam/Sim.DEL.INS.DUP.INV.25x.cuteSV.sort.bam /public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/genome.fa /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.DEL.INS.DUP.INV.25x/ngmlr/vcf/cuteSV/Sim.DEL.INS.DUP.INV.25x.vcf /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.DEL.INS.DUP.INV.25x/ngmlr/vcf/cuteSV
