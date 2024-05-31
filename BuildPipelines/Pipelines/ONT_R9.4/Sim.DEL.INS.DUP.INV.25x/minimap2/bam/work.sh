/home/lz/miniconda3/bin/minimap2 /public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/genome.fa /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/fastq/ONT_R9.4_VISOR_LASeR_Sim_DEL.DUP.INS.INV.fq.gz -a -z 600,200 -x map-ont --MD -Y -o /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.DEL.INS.DUP.INV.25x/minimap2/bam/Sim.DEL.INS.DUP.INV.25x.sam -R '@RG\tID:Sim.DEL.INS.DUP.INV.25x\tSM:Sim.DEL.INS.DUP.INV.25x' -t 20
samtools view -@ 20 -bS /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.DEL.INS.DUP.INV.25x/minimap2/bam/Sim.DEL.INS.DUP.INV.25x.sam -o /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.DEL.INS.DUP.INV.25x/minimap2/bam/Sim.DEL.INS.DUP.INV.25x.bam
samtools sort -@ 5 /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.DEL.INS.DUP.INV.25x/minimap2/bam/Sim.DEL.INS.DUP.INV.25x.bam -o /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.DEL.INS.DUP.INV.25x/minimap2/bam/Sim.DEL.INS.DUP.INV.25x.sort.bam
samtools index /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.DEL.INS.DUP.INV.25x/minimap2/bam/Sim.DEL.INS.DUP.INV.25x.sort.bam
rm /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.DEL.INS.DUP.INV.25x/minimap2/bam/Sim.DEL.INS.DUP.INV.25x.sam /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R9.4/Sim.DEL.INS.DUP.INV.25x/minimap2/bam/Sim.DEL.INS.DUP.INV.25x.bam
