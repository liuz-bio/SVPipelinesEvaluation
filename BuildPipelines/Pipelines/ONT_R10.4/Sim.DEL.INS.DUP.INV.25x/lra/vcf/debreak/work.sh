ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.DEL.INS.DUP.INV.25x/lra/bam/Sim.DEL.INS.DUP.INV.25x.sort.bam /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.DEL.INS.DUP.INV.25x/lra/bam/Sim.DEL.INS.DUP.INV.25x.debreak.sort.bam
ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.DEL.INS.DUP.INV.25x/lra/bam/Sim.DEL.INS.DUP.INV.25x.sort.bam.bai /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.DEL.INS.DUP.INV.25x/lra/bam/Sim.DEL.INS.DUP.INV.25x.debreak.sort.bam.bai
/home/lz/software/sv_callers/DeBreak/debreak --bam /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.DEL.INS.DUP.INV.25x/lra/bam/Sim.DEL.INS.DUP.INV.25x.sort.bam --depth 25 --min_support 2 --thread 20 --min_size 30 -r /public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/lra_ONT/genome.fa -o /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.DEL.INS.DUP.INV.25x/lra/vcf/debreak/Sim.DEL.INS.DUP.INV.25x 
rm /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.DEL.INS.DUP.INV.25x/lra/vcf/debreak/Sim.DEL.INS.DUP.INV.25x/*.temp
cp /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.DEL.INS.DUP.INV.25x/lra/vcf/debreak/Sim.DEL.INS.DUP.INV.25x/debreak.vcf /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.DEL.INS.DUP.INV.25x/lra/vcf/debreak/Sim.DEL.INS.DUP.INV.25x/Sim.DEL.INS.DUP.INV.25x.vcf
cp /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.DEL.INS.DUP.INV.25x/lra/vcf/debreak/Sim.DEL.INS.DUP.INV.25x/debreak.vcf /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.DEL.INS.DUP.INV.25x/lra/vcf/debreak/Sim.DEL.INS.DUP.INV.25x.vcf
