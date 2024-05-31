ngmlr -t 20 --rg-id Sim.DEL.INS.DUP.INV.25x --rg-sm Sim.DEL.INS.DUP.INV.25x --rg-pl Nanopore -x ont -r /public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/genome.fa -q /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/fastq/ONT_R10.4_VISOR_LASeR_Sim_DEL.DUP.INS.INV.fq.gz -o /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.DEL.INS.DUP.INV.25x/ngmlr/bam/Sim.DEL.INS.DUP.INV.25x.sam
sambamba view -f bam -S -o /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.DEL.INS.DUP.INV.25x/ngmlr/bam/Sim.DEL.INS.DUP.INV.25x.bam /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.DEL.INS.DUP.INV.25x/ngmlr/bam/Sim.DEL.INS.DUP.INV.25x.sam
bamaddrg -r Sim.DEL.INS.DUP.INV.25x -s Sim.DEL.INS.DUP.INV.25x -b /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.DEL.INS.DUP.INV.25x/ngmlr/bam/Sim.DEL.INS.DUP.INV.25x.bam >/public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.DEL.INS.DUP.INV.25x/ngmlr/bam/Sim.DEL.INS.DUP.INV.25x.rg.bam
samtools sort -@ 5 /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.DEL.INS.DUP.INV.25x/ngmlr/bam/Sim.DEL.INS.DUP.INV.25x.rg.bam -o /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.DEL.INS.DUP.INV.25x/ngmlr/bam/Sim.DEL.INS.DUP.INV.25x.sort.bam
samtools index /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.DEL.INS.DUP.INV.25x/ngmlr/bam/Sim.DEL.INS.DUP.INV.25x.sort.bam
