ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/minimap2/bam/Sim.TRA.25x.sort.bam /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/minimap2/bam/Sim.TRA.25x.svision.sort.bam
ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/minimap2/bam/Sim.TRA.25x.sort.bam.bai /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/minimap2/bam/Sim.TRA.25x.svision.sort.bam.bai
/usr/bin/docker run -i -u 1023:1002 --name=ONT_R10.4_minimap2_svision_Sim.TRA.25x --rm -v /public3/SVDataset_lz/backup/2023_11_14/:/public3/SVDataset_lz/backup/2023_11_14/ jiadongxjtu/svision:latest SVision -o /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/minimap2/vcf/svision/Sim.TRA.25x/ -b /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/ONT_R10.4/Sim.TRA.25x/minimap2/bam/Sim.TRA.25x.svision.sort.bam -m /public3/SVDataset_lz/backup/2023_11_14/SVision_old/model/svision-cnn-model.ckpt -g /public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/genome.fa -n Sim.TRA.25x -s 2 --graph --qname -t 20
