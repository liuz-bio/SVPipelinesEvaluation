ls Pacbio_minimap2_*.rule|cut -d\_ -f3|cut -d\. -f1|while read id;do(cp Pacbio_minimap2_$id.rule Pacbio_ngmlr_$id.rule;sed -i "s/minimap2/ngmlr/g" Pacbio_ngmlr_$id.rule);done

echo -e "Nanopore\nPacbio"|while read ida;do(echo -e "minialign\nwinnowmap\nBWA-MEM\nlordfast\nLAST"|while read idb;do(touch $ida'_'$idb'.rule');done);done
echo -e "Nanopore\nPacbio"|while read ida;do(echo -e "minialign\nwinnowmap\nBWA-MEM\nLAST\ngraphmap\nlordfast"|while read idb;do(echo -e "pbsv\ncuteSV\ndelly\nsvim\nsniffles\npicky\nnanovar\nNanoSV\ndebreak"|while read idc;do( touch $ida'_'$idb'_'$idc'.rule');done);done);done

echo -e "Nanopore\nPacbio"|while read ida;do(echo -e "minialign\nwinnowmap\nBWA-MEM\nMECAT\nLAST"|while read idb;do(echo $ida"="$idb"=debreak++ln::-s::bamPath/sampleId.sort.bam::bamPath/sampleId.svTools.sort.bam++ln::-s::bamPath/sampleId.sort.bam.bai::bamPath/sampleId.svTools.sort.bam.bai++/home/lz/software/sv_callers/DeBreak/debreak::--bam::bamPath/sampleId.sort.bam::--depth::SqDepth::--min_support::SReads::--thread::thread::--min_size::SvLength::-r::ref.fa::-o::outPath/sampleId::-p::sampleId" >$ida\_$idb\_debreak.rule);done);done

echo -e "Nanopore\nPacbio"|while read ida;do(echo -e "pbsv\ncuteSV\ndelly\nsvim\nsniffles\npicky\nnanovar\nNanoSV\ndebreak\nnpInv\nSVLR"|while read idc;do(echo -e "minialign\nwinnowmap\nBWA-MEM\nMECAT\nLAST"|while read idb;do(cat $ida\_minimap2_$idc.rule|sed "s/minimap2/$idb/g" > $ida\_$idb\_$idc.rule);done);done);done


echo -e "Nanopore\nPacbio"|while read ida;do(echo -e "minialign\nwinnowmap\nBWA-MEM\nLAST\ngraphmap\nlordfast"|while read idb;do(echo -e "pbsv\ncuteSV\ndelly\nsvim\nsniffles\npicky\nnanovar\nNanoSV\ndebreak"|while read idc;do( touch $ida'_'$idb'_'$idc'.rule');done);done);done
 
echo -e "Nanopore\nPacbio"|while read ida;do(echo -e "pbsv\ncuteSV\ndelly\nsvim\nsniffles\npicky\nnanovar\nNanoSV\ndebreak"|while read idc;do(echo -e "ngmlr\nminialign\nwinnowmap\nBWA-MEM\nLAST\ngraphmap\nlordfast"|while read idb;do(cat $ida\_minimap2_$idc.rule|sed "s/minimap2/$idb/g" > $ida\_$idb\_$idc.rule);done);done);done

echo -e "Nanopore\nPacbio"|while read ida;do(echo -e "pbsv\ncuteSV\ndelly\nsvim\nsniffles\npicky\nnanovar\nNanoSV\ndebreak"|while read idc;do( cat $ida\_minimap2_$idc.rule|sed "s/minimap2/pbmm2/g" >  $ida\_pbmm2_$idc.rule);done);done
