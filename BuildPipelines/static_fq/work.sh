cat ../Sample.tab |awk '{print "nohup seqkit stat -j 10  -a -o "$2"."$3".tab "$5" &" }' >a.sh
cat /home/lz/Data_sequence/2023_11_14/working/Sample.tab|awk '{print "nohup seqkit stat -j 10  -a -o "$2"."$3".tab "$5" &" }' >b.sh
cat /public3/SVDataset_lz/backup/2023_11_14/HG002_chr20/Sample.tab|grep 25x |awk '{print "nohup seqkit stat -j 10  -a -o "$2"."$3".tab "$5" &" }' >c.sh
