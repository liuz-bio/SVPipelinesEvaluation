ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/minimap2/bam/Sim.TRA.25x.sort.bam /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/minimap2/bam/Sim.TRA.25x.NanoSV.sort.bam
ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/minimap2/bam/Sim.TRA.25x.sort.bam.bai /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/minimap2/bam/Sim.TRA.25x.NanoSV.sort.bam.bai
/home/lz/anaconda3/bin/NanoSV -t 20 -b /public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/hg38_genome_sample.bed -o /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/minimap2/vcf/NanoSV/Sim.TRA.25x.vcf /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/minimap2/bam/Sim.TRA.25x.NanoSV.sort.bam
