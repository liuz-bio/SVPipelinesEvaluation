ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CCS/Sim.DEL.INS.DUP.INV.25x/pbmm2/bam/Sim.DEL.INS.DUP.INV.25x.sort.bam /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CCS/Sim.DEL.INS.DUP.INV.25x/pbmm2/bam/Sim.DEL.INS.DUP.INV.25x.cuteSV2.sort.bam
ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CCS/Sim.DEL.INS.DUP.INV.25x/pbmm2/bam/Sim.DEL.INS.DUP.INV.25x.sort.bam.bai /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CCS/Sim.DEL.INS.DUP.INV.25x/pbmm2/bam/Sim.DEL.INS.DUP.INV.25x.cuteSV2.sort.bam.bai
/home/lz/miniconda3/envs/cutesv/bin/cuteSV --genotype --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5 -t 20 -s 2 -l 30 /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CCS/Sim.DEL.INS.DUP.INV.25x/pbmm2/bam/Sim.DEL.INS.DUP.INV.25x.cuteSV2.sort.bam /public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/genome.fa /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CCS/Sim.DEL.INS.DUP.INV.25x/pbmm2/vcf/cuteSV2/Sim.DEL.INS.DUP.INV.25x.vcf /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CCS/Sim.DEL.INS.DUP.INV.25x/pbmm2/vcf/cuteSV2
