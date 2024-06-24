le F1_recall_precision_geneType.tab|cut -f 2|sort -u|while read id;do(grep $id F1_recall_precision_geneType.tab|grep Pacbio|grep -P "\.15\t"|grep F1|grep DEL|cut -d\, -f 1|cut -f7|perl -pe "s/\n/,/g"|sed "s/,$/\n/g"|sed "s/^/$id,/g");done

le F1_recall_precision_geneType.tab|cut -f 2|sort -u|while read id;do(grep $id F1_recall_precision_geneType.tab|grep Pacbio|grep -P "\.15\t"|grep F1|grep INS|cut -d\, -f 1|cut -f7|perl -pe "s/\n/,/g"|sed "s/,$/\n/g"|sed "s/^/$id,/g");done

le F1_recall_precision_geneType.tab|cut -f 2|sort -u|while read id;do(grep $id F1_recall_precision_geneType.tab|grep Pacbio|grep -P "\.15\t"|grep F1|grep INV|cut -d\, -f 1|cut -f7|perl -pe "s/\n/,/g"|sed "s/,$/\n/g"|sed "s/^/$id,/g");done

le F1_recall_precision_geneType.tab|cut -f 2|sort -u|while read id;do(grep $id F1_recall_precision_geneType.tab|grep Pacbio|grep -P "\.15\t"|grep F1|grep DUP|cut -d\, -f 1|cut -f7|perl -pe "s/\n/,/g"|sed "s/,$/\n/g"|sed "s/^/$id,/g");done


cat F1_recall_precision.tab|grep -v summary.txt|cut -f 2|sort -u|while read id;do(grep $id F1_recall_precision.tab |grep -v summary.txt|grep Pacbio|grep -P "\.15\t"|grep F1|grep DEL|cut -d\, -f 1|cut -f7|perl -pe "s/\n/,/g"|sed "s/,$/\n/g"|sed "s/^/$id,/g");done

cat F1_recall_precision.tab|grep -v summary.txt|cut -f 2|sort -u|while read id;do(grep $id F1_recall_precision.tab |grep -v summary.txt|grep Pacbio|grep -P "\.15\t"|grep F1|grep INS|cut -d\, -f 1|cut -f7|perl -pe "s/\n/,/g"|sed "s/,$/\n/g"|sed "s/^/$id,/g");done

cat F1_recall_precision.tab|grep -v summary.txt|cut -f 2|sort -u|while read id;do(grep $id F1_recall_precision.tab |grep -v summary.txt|grep Pacbio|grep -P "\.15\t"|grep F1|grep INV|cut -d\, -f 1|cut -f7|perl -pe "s/\n/,/g"|sed "s/,$/\n/g"|sed "s/^/$id,/g");done

cat F1_recall_precision.tab|grep -v summary.txt|cut -f 2|sort -u|while read id;do(grep $id F1_recall_precision.tab |grep -v summary.txt|grep Pacbio|grep -P "\.15\t"|grep F1|grep DUP|cut -d\, -f 1|cut -f7|perl -pe "s/\n/,/g"|sed "s/,$/\n/g"|sed "s/^/$id,/g");done

paste <(cat F1_recall_precision.tab|grep -v summary.txt|cut -f 2|sort -u|while read id;do(grep $id F1_recall_precision.tab |grep -v summary.txt|grep Pacbio|grep -P "\.15\t"|grep F1|grep DEL|cut -d\, -f 1|cut -f7|awk '{sum+=$1} END {print sum/NR}'|sed "s/^/$id,/g");done) <(cat F1_recall_precision.tab|grep -v summary.txt|cut -f 2|sort -u|while read id;do(grep $id F1_recall_precision.tab |grep -v summary.txt|grep Pacbio|grep -P "\.15\t"|grep F1|grep INS|cut -d\, -f 1|cut -f7|awk '{sum+=$1} END {print sum/NR}'|sed "s/^/$id,/g");done)  <(cat F1_recall_precision.tab|grep -v summary.txt|cut -f 2|sort -u|while read id;do(grep $id F1_recall_precision.tab |grep -v summary.txt|grep Pacbio|grep -P "\.15\t"|grep F1|grep INV|cut -d\, -f 1|cut -f7|awk '{sum+=$1} END {print sum/NR}'|sed "s/^/$id,/g");done) <(cat F1_recall_precision.tab|grep -v summary.txt|cut -f 2|sort -u|while read id;do(grep $id F1_recall_precision.tab |grep -v summary.txt|grep Pacbio|grep -P "\.15\t"|grep F1|grep DUP|cut -d\, -f 1|cut -f7|awk '{sum+=$1} END {print sum/NR}'|sed "s/^/$id,/g");done) |sed 's/\t/,/g'|cut -d\, -f 1,2,4,6,8|sed "1i,DEL,INS,INV,DUP"
