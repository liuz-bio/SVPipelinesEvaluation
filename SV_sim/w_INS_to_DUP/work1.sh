cat INS_to_DUP.tab|grep -E "CCS|CLR|ONT"|grep -v \-|grep -vE "pbmm2\tnanovar"|grep -vP 'lra\tsniffles'|grep -vP "lra\tpicky"|grep -vP "lra\tdelly"|grep -vP "lra\tnanovar"|grep -vP "lra\tnanosv"|grep -vP 'lra\tpbsv'|awk -v OFS='\t' '{print $1,$2"-"$3, $5}'

cat DUP_to_INS.tab|grep -E "CCS|CLR|ONT"|grep -v \-|grep -vE "pbmm2\tnanovar"|grep -vP 'lra\tsniffles'|grep -vP "lra\tpicky"|grep -vP "lra\tdelly"|grep -vP "lra\tnanovar"|grep -vP "lra\tnanosv"|grep -vP 'lra\tpbsv'|awk -v OFS='\t' '{print $1,$2"-"$3, $5}'
