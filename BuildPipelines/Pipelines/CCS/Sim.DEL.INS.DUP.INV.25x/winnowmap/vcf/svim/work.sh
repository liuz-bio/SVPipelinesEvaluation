ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CCS/Sim.DEL.INS.DUP.INV.25x/winnowmap/bam/Sim.DEL.INS.DUP.INV.25x.sort.bam /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CCS/Sim.DEL.INS.DUP.INV.25x/winnowmap/bam/Sim.DEL.INS.DUP.INV.25x.svim.sort.bam
ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CCS/Sim.DEL.INS.DUP.INV.25x/winnowmap/bam/Sim.DEL.INS.DUP.INV.25x.sort.bam.bai /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CCS/Sim.DEL.INS.DUP.INV.25x/winnowmap/bam/Sim.DEL.INS.DUP.INV.25x.svim.sort.bam.bai
/home/lz/anaconda3/bin/svim alignment --minimum_depth 2 --min_sv_size 30 /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CCS/Sim.DEL.INS.DUP.INV.25x/winnowmap/vcf/svim/Sim.DEL.INS.DUP.INV.25x/Sim.DEL.INS.DUP.INV.25x /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CCS/Sim.DEL.INS.DUP.INV.25x/winnowmap/bam/Sim.DEL.INS.DUP.INV.25x.svim.sort.bam /public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/genome.fa
