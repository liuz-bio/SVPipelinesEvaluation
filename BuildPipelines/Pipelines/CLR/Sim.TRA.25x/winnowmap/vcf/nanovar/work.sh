ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CLR/Sim.TRA.25x/winnowmap/bam/Sim.TRA.25x.sort.bam /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CLR/Sim.TRA.25x/winnowmap/bam/Sim.TRA.25x.nanovar.sort.bam
ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CLR/Sim.TRA.25x/winnowmap/bam/Sim.TRA.25x.sort.bam.bai /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CLR/Sim.TRA.25x/winnowmap/bam/Sim.TRA.25x.nanovar.sort.bam.bai
/home/lz/miniconda3/envs/nanovar/bin/nanovar -x pacbio-clr -c 2 -t 20 -l 30 /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CLR/Sim.TRA.25x/winnowmap/bam/Sim.TRA.25x.nanovar.sort.bam /public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/genome.fa /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CLR/Sim.TRA.25x/winnowmap/vcf/nanovar/Sim.TRA.25x/Sim.TRA.25x
