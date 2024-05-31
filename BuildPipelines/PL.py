import os
import sys

Platforms = ['Nanopore','Pacbio']

#"pbmm2" "blasr"
Map_tools = {'Nanopore':["pbmm2","ngmlr","lra","minimap2","winnowmap"],
             'Pacbio':["pbmm2","ngmlr","lra","minimap2","winnowmap"]}

Map_index ={'Other':'/public3/SVDataset_lz/backup/2023_11_14/Genome/hg38','lra':{'CCS':'/public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/lra_CCS','CLR':'/public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/lra_CLR','ONT':'/public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/lra_ONT'}}

type_steps = ['bam','vcf']

Call_tools = {'Nanopore':{'ONT':["cuteSV","cuteSV2","delly","NanoSV","nanovar","pbsv","sniffles","sniffles2","svim","picky","debreak","svision"]},
              'Pacbio':{'CLR':["cuteSV","cuteSV2","delly","NanoSV","nanovar","pbsv","sniffles","sniffles2","svim","picky","debreak","svision"],
                        'CCS':["cuteSV","cuteSV2","delly","NanoSV","nanovar","pbsv","sniffles","sniffles2","svim","picky","debreak","svision"]}}

Samples = {'Nanopore':{'ONT':{},'ONT_R10.4':{},'ONT_R9.4':{}},'Pacbio':{'CCS':{},'CLR':{}}}
for i in open('Sample.tab','r'):
    i=i.strip().split('\t')
    Samples[i[0]][i[1]][i[2]]=i[3:]

print(Samples)
'''
Samples = {'Nanopore':{'nanopore_del_5x':'/home/lz/Data_sequence/2021_4_2/visortest/simulation/nanopore_del_5x/sim.srt.bam',
                       'nanopore_ins_5x':'/home/lz/Data_sequence/2021_4_2/work/Pacbio/nanopore_ins_5x.fq'},
           'Pacbio':{'pacbio_del_5x':'/home/lz/Data_sequence/2021_4_2/visortest/simulation/pacbio_del_5x/sim.srt.bam',
                      'pacbio_ins_5x':'/home/lz/Data_sequence/2021_4_2/work/Pacbio/pacbio_ins_5x.fq'}}
'''

rule_info = {'sampleId':None,'ref.fa':'/public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/genome.fa','sample.fq':None,'ThreadNumber':20,'SReads':2,'SvLength':30,"svisionID":'svision'}


def replace_rule(rules):
    Orules = []
    for rule in rules:
        for k,v in rule_info.items():
            rule = rule.replace(str(k),str(v))
        rule = rule.replace('::',' ')
        Orules.append(rule)
    return Orules

Rules={}
for line in open(sys.argv[1],'r'):
    if ('#####' in line) or ('-----' in line):
        continue
    line = line.strip().split('++')
    Id = line[0]
    Value = line[1:]
    Rules[Id] = Value
print(Rules.keys())

#OutW = open('genomeRefIndex.sh','w')
#for MapI_k,MapI_v in Map_index.items():
#    rule_info['GenomeIndexPath'] = MapI_v
#    OutW.write('\n'.join(replace_rule(Rules['='.join([MapI_k,'Index'])]))+'\n')
#OutW.close()
#print(Rules)

outRules = {}
for platf in Platforms:
    if platf not in outRules:
        outRules[platf] = {}
        
    for paltData, samples in Samples[platf].items():
        if paltData.split("_")[0] not in outRules[platf]:
            outRules[platf][paltData.split("_")[0]] = {}
        for sampleId, sample in samples.items():
            rule_info['sampleId'] = sampleId
            rule_info['SqDepth']  = 25 #sampleId.split('.')[-1].replace('x','')
 
            if '.bam' in sample[-1]:
                rule_info['sample.fq'] = 'sample'
                rule_info['SampleID.fq'] = sample[-1]
                rule_info['DataRatio'] = sample[0]
                try:
                    bam2fq_rule = Rules['='.join([platf,paltData.plit("_")[0],'simulation'])]
                except:
                    continue
            else:
                rule_info['sample.fq'] = sample[-1]
                rule_info['DataRatio'] = sample[0]
 
            if  sampleId not in outRules[platf][paltData.split("_")[0]]:
                outRules[platf][paltData.split("_")[0]][sampleId] = {}
            
            for mapt in Map_tools[platf]:
  
                if mapt in Map_index.keys():
                    rule_info['GenomeIndexPath'] = Map_index[mapt][paltData.split("_")[0]]
                else:
                    rule_info['GenomeIndexPath'] = Map_index['Other']
 
                if  mapt not in outRules[platf][paltData.split("_")[0]][sampleId]:
                    outRules[platf][paltData.split("_")[0]][sampleId][mapt] = {}
 
                mapt_id = '='.join([platf,paltData.split("_")[0],mapt])
                try:
                    map_rule = Rules[mapt_id]
                except:
                    continue
            #print(map_rule)

                for step in type_steps:
                    if step=='bam':
                        bamPath = os.popen('pwd').read().strip()+'/'+'/'.join(['Pipelines',paltData,sampleId,mapt,step])
                        rule_info['bamPath'] = bamPath

                        if not os.path.exists(bamPath):
                            os.makedirs(bamPath)

                        outW = open('%s/work.sh'%bamPath,'w')
                        if rule_info['sample.fq'] == 'sample':
                            rule_info['sample.fq'] = bamPath+'/'+sampleId+'.fq'
                            outW.write('\n'.join(replace_rule(bam2fq_rule+map_rule+['rm::bamPath/sampleId.fq']))+'\n')
                            rule_info['sample.fq'] = 'sample'
                        else:
                            outW.write('\n'.join(replace_rule(map_rule))+'\n')
                        outW.close()

                    elif step=='vcf':
                        print("paltData:",paltData, platf,paltData.split("_")[0],mapt,sampleId)
                        for callt in Call_tools[platf][paltData.split("_")[0]]:
                            rule_info["svTools"] = callt
                            callt_id = '='.join([platf,paltData.split("_")[0],mapt,callt])
 
                            try:
                                callt_rule = Rules[callt_id]
                            except:
                                continue
                            if  True:#callt not in outRules[platf][paltData.split("_")[0]][sampleId][mapt]:
                                outPath = os.popen('pwd').read().strip()+'/'+'/'.join(["Pipelines",paltData,sampleId,mapt,step,callt])
                                print(outPath)
                                rule_info['outPath'] = outPath
                                rule_info['svisionID'] = '_'.join([paltData,mapt,callt,sampleId])

                                if not os.path.exists(outPath):
                                    os.makedirs(outPath)

                                outW = open('%s/work.sh'%outPath,'w')
                                outW.write('\n'.join(replace_rule(callt_rule))+'\n')
                                outW.close()
                                outRules[platf][paltData.split("_")[0]][sampleId][mapt][callt] = map_rule+callt_rule


