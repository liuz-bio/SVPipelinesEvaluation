ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CLR/Sim.TRA.25x/lra/bam/Sim.TRA.25x.sort.bam /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CLR/Sim.TRA.25x/lra/bam/Sim.TRA.25x.svision.sort.bam
ln -s /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CLR/Sim.TRA.25x/lra/bam/Sim.TRA.25x.sort.bam.bai /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CLR/Sim.TRA.25x/lra/bam/Sim.TRA.25x.svision.sort.bam.bai
/usr/bin/docker run -i -u 1023:1002 --name=CLR_lra_svision_Sim.TRA.25x --rm -v /public3/SVDataset_lz/backup/2023_11_14/:/public3/SVDataset_lz/backup/2023_11_14/ jiadongxjtu/svision:latest SVision -o /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CLR/Sim.TRA.25x/lra/vcf/svision/Sim.TRA.25x/ -b /public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/sim_pipelins/Pipelines/CLR/Sim.TRA.25x/lra/bam/Sim.TRA.25x.svision.sort.bam -m /public3/SVDataset_lz/backup/2023_11_14/SVision_old/model/svision-cnn-model.ckpt -g /public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/lra_CLR/genome.fa -n Sim.TRA.25x -s 2 --graph --qname -t 20
