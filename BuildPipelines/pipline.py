import os
import sys

Platforms = ['Nanopore',"Pacbio"]

#"pbmm2" "blasr"
Map_tools = {'Nanopore':["pbmm2","ngmlr","minimap2","winnowmap",'lordfast'],
             'Pacbio':["pbmm2","lordfast","ngmlr","minimap2","winnowmap"]}

Map_index ={'Other':'/home/lz/Data_sequence/Ref/Other','lordfast':'/home/lz/Data_sequence/Ref/lordfast'}

type_steps = ['bam','vcf']

Call_tools = {'Nanopore':["cuteSV","delly","NanoSV","nanovar","pbsv","sniffles","svim","picky","debreak"],
              'Pacbio':["cuteSV","delly","NanoSV","nanovar","pbsv","sniffles","svim","picky","debreak"]}

Samples = {'Nanopore':{},'Pacbio':{}}
for i in open('Sample.tab','r'):
    i=i.strip().split('\t')
    Samples[i[0]][i[1]]=i[2]



rule_info = {'sampleId':None,'ref.fa':'/home/lz/Data_sequence/Ref/hs37d5.fa','sample.fq':None,'ThreadNumber':20,'SReads':2,'SvLength':30}


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

OutW = open('genomeRefIndex.sh','w')
for MapI_k,MapI_v in Map_index.items():
    rule_info['GenomeIndexPath'] = MapI_v
    OutW.write('\n'.join(replace_rule(Rules['='.join([MapI_k,'Index'])]))+'\n')
OutW.close()
#print(Rules)

outRules = {}
for platf in Platforms:
    if platf not in outRules:
        outRules[platf] = {}
        
    for sampleId, sample in Samples[platf].items():
        rule_info['sampleId'] = sampleId
        rule_info['SqDepth']  = sampleId.split('_')[-1][:-1] 

        if 'sim.srt.bam' in sample:
            rule_info['sample.fq'] = 'sample'
            rule_info['simulation.bam'] = sample
            try:
                bam2fq_rule = Rules['='.join([platf,'simulation'])]
            except:
                continue
        else:
            rule_info['sample.fq'] = sample
 
        if  sampleId not in outRules[platf]:
            outRules[platf][sampleId] = {}
            
        for mapt in Map_tools[platf]:
  
            if mapt in Map_index.keys():
                rule_info['GenomeIndexPath'] = Map_index[mapt]
            else:
                rule_info['GenomeIndexPath'] = Map_index['Other']

            if  mapt not in outRules[platf][sampleId]:
                outRules[platf][sampleId][mapt] = {}

            mapt_id = '='.join([platf,mapt])
            try:
                map_rule = Rules[mapt_id]
            except:
                continue
            #print(map_rule)

            for step in type_steps:
                if step=='bam':
                    bamPath = os.popen('pwd').read().strip()+'/'+'/'.join([platf,sampleId,mapt,step])
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
                    for callt in Call_tools[platf]:
                        rule_info["svTools"] = callt
                        callt_id = '='.join([platf,mapt,callt])

                        try:
                            callt_rule = Rules[callt_id]
                        except:
                            continue

                        if  callt not in outRules[platf][sampleId][mapt]:
                            outPath = os.popen('pwd').read().strip()+'/'+'/'.join([platf,sampleId,mapt,step,callt])
                            print(outPath)
                            rule_info['outPath'] = outPath

                            if not os.path.exists(outPath):
                                os.makedirs(outPath)

                            outW = open('%s/work.sh'%outPath,'w')
                            outW.write('\n'.join(replace_rule(callt_rule))+'\n')
                            outW.close()
                            outRules[platf][sampleId][mapt][callt] = map_rule+callt_rule

