ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CLR/Sim.TRA.25x/ngmlr/bam/Sim.TRA.25x.sort.bam /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CLR/Sim.TRA.25x/ngmlr/bam/Sim.TRA.25x.pbsv.sort.bam
ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CLR/Sim.TRA.25x/ngmlr/bam/Sim.TRA.25x.sort.bam.bai /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CLR/Sim.TRA.25x/ngmlr/bam/Sim.TRA.25x.pbsv.sort.bam.bai
/home/lz/anaconda3/bin/pbsv discover /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CLR/Sim.TRA.25x/ngmlr/bam/Sim.TRA.25x.pbsv.sort.bam /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CLR/Sim.TRA.25x/ngmlr/vcf/pbsv/Sim.TRA.25x.svsig.gz --tandem-repeats ~/Data_sequence/work_cuteSV/align/ref/pbsv/annotations/human_GRCh38_no_alt_analysis_set.trf.bed
/home/lz/anaconda3/bin/pbsv call --gt-min-reads 2 -m 30 -j 20 /public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/genome.fa /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CLR/Sim.TRA.25x/ngmlr/vcf/pbsv/Sim.TRA.25x.svsig.gz /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CLR/Sim.TRA.25x/ngmlr/vcf/pbsv/Sim.TRA.25x.vcf
