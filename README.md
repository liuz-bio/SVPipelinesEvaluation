## Comprehensive and deep evaluation of structural variation detection pipelines with third-generation sequencing data 

### created simulation reads

```shell
VISOR HACk -g hg38.fa -b sim_DEL.DUP.INS.INV_h1.bed sim_DEL.DUP.INS.INV_h2.bed -o VISOR_HACk_Sim_DEL.DUP.INS.INV 

VISOR HACk -g hg38.fa -b sim_TRA.bed -o VISOR_HACk_Sim_TRA 

echo -e "pacbio2016\npacbio2021\nnanopore2020\nnanopore2023"|while read id;do(VISOR LASeR -g -s VISOR_HACk_Sim_DEL.DUP.INS.INV -b LASeR100.bed -o $id\_VISOR_LASeR_Sim_DEL.DUP.INS.INV --coverage 25   --threads 80 --tag --read_type nanopore --error_model $id --qscore_model $id );done

echo -e "pacbio2016\npacbio2021\nnanopore2020\nnanopore2023"|while read id;do(VISOR LASeR -g -s VISOR_HACk_Sim_TRA -b LASeR100.bed -o $id\_VISOR_LASeR_Sim_ TRA --coverage 25   --threads 80 --tag --read_type nanopore --error_model $id --qscore_model $id );done
```



### created SVs detection pipelines 

```shell
cd BuildPipelines
# merge aligner and caller rule
cat Rule/*.rule >Rules.txt
# build SV detection pipelines
python PL.py Rules.txt

```

#### Pipelines directory structure

```shell
Pipelines/
├── CCS
│   ├── Sim.DEL.INS.DUP.INV.25x
│   │   ├── minimap2
│   │   │   ├── bam
│   │   │   │   └── work.sh
│   │   │   └── vcf
│   │   │       ├── cuteSV
│   │   │       │   └── work.sh
│   │   │       ├── cuteSV2
│   │   │       │   └── work.sh
│   │   │       ├── debreak
│   │   │       │   └── work.sh
│   │   │       ├── delly
│   │   │       │   └── work.sh
│   │   │       ├── NanoSV
│   │   │       │   └── work.sh
│   │   │       ├── nanovar
│   │   │       │   └── work.sh
│   │   │       ├── pbsv
│   │   │       │   └── work.sh
│   │   │       ├── picky
│   │   │       │   └── work.sh
│   │   │       ├── sniffles
│   │   │       │   └── work.sh
│   │   │       ├── sniffles2
│   │   │       │   └── work.sh
│   │   │       ├── svim
│   │   │       │   └── work.sh
│   │   │       └── svision
│   │   │           └── work.sh
│   │   ├── winnowmap
│   │   │   ├── bam
│   │   │   │   └── work.sh
│   │   │   └── vcf
│   │   │       ├── cuteSV
│   │   │       │   └── work.sh
│   │   │       ├── cuteSV2
│   │   │       │   └── work.sh
│   │   │       ├── debreak
│   │   │       │   └── work.sh
│   │   │       ├── delly
│   │   │       │   └── work.sh
│   │   │       ├── NanoSV
│   │   │       │   └── work.sh
│   │   │       ├── nanovar
│   │   │       │   └── work.sh
│   │   │       ├── pbsv
│   │   │       │   └── work.sh

......
```

### Analysis of pipelines SVs detection results

```shell
cd SV_sim/sim_real
python main.py
# readInfo： function standardizes the SV test results of all pipelines. 
# readSim： builds the Vcf file as the benchmark based on the bed file of the simulation data. 
# evaluationSV： evaluates precision, recall, F1 and other attributes of pipeline
```

