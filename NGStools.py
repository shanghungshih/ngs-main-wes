# -*- coding: utf-8 -*-
"""
Created on Thu May 10 11:42:21 2018

@author: Shang-Hung, Shih
"""

import os
import re

###main docker run in deamon
def getID():
    yn = input('Show all patient ID of sub directory?(Y/n) (ex. y): ')
    if yn.strip().lower() == 'y':
        path = input('Please enter the directory name (ex. 66xWES): ')
        os.system('ls -lt %s > getID.txt' %(path))
        name = []
        with open('getID.txt', 'r') as f:
            for i in f.readlines():
                if i.split(' ')[-1:][0].replace('\n', '').endswith('T') is True:
                    name.append(i.split(' ')[-1:][0].replace('\n', '').replace('T', ''))
        out = ('%s' %(name))
        print('>>>>>> total patient of %s : %s' %(path, len(name)))
        print('>>>>>> recommanded minimum space for %s : %s T' %(path, len(name)*200/1024))
        os.system('df -h')
        print(out.replace('[', '').replace(']', '').replace("'", ""))
        os.system('rm getID.txt')

def getAllVCF():
    name = input('Please enter subproject: ')
    anno = os.getcwd()+'/annotation'
    cmd = ('cp data/%s/*/*vcf %s' %(name, anno))
    os.system(cmd)
        
def rmSAM():
    os.system('rm data/*/*/*sam data/*/*/*sai')
    print('>>>>>> all sam file have been removed.')

def rawdataRename():
    argv = input('Do you want to rename your rawdata directory?(Y/n) (ex. y): ')
    argv = argv.lower()
    if argv == 'y':
        subproject = input('Please enter folder to be modified (ex. 66xWES): ')
        total = input('Please enter name to be modified (ex. OC_712, OC_512, OC_345): ')
        prefix = input('Please enter head -h or tail -t (ex. -h OC_): ')

        total_1 = total.split(',')

        for patient in total_1:
            patient = patient.strip()

            if prefix.split(' ')[0].split('-')[1] == 'h':
                patient = patient.split(prefix.split(' ')[1])[1]
                normal = prefix.split(' ')[1]+patient+'N'
                tumor = prefix.split(' ')[1]+patient+'T'
                new_normal = patient+'N'
                new_tumor = patient+'T'

            elif prefix.split(' ')[0].split('-')[1] == 't':
                patient = patient.split(prefix.split(' ')[1])[0]
                normal = patient+'N'+prefix.split(' ')[1]
                tumor = patient+'T'+prefix.split(' ')[1]
                new_normal = patient+'N'
                new_tumor = patient+'T'

            os.system('mv %s/%s %s/%s' %(subproject, normal, subproject, new_normal))
            os.system('mv %s/%s %s/%s' %(subproject, tumor, subproject, new_tumor))
        print('>>>>>> %s have been rename.' %(total_1))

def enterData():
    subproject = input('Please enter your subproject name (ex. 66xWES): ').strip()
    patient = input('Please enter patient ID to add in analysis queue (ex. 599, 631, 632, 652): ').split(',')

    work = []
    for i in patient:
        work.append([subproject, i.strip()])
    work = tuple(work)
    return work

def ReferenceIndex(refPath, fasta, dbsnp, cosmic):
    if os.path.exists('ref_data/'+fasta+'.fasta.fa') is False:
        bwa(refPath, 'bwa index '+fasta+'.fasta')
    if os.path.exists('ref_data/'+fasta+'.fasta.fai') is False:
        samtools(refPath, 'samtools faidx '+fasta+'.fasta')
    if os.path.exists('ref_data/'+fasta+'.dict') is False:
        picard(refPath, 'java -jar picard.jar CreateSequenceDictionary R=data/'+fasta+'.fasta O=data/'+fasta+'.dict')
    if os.path.exists('ref_data/'+dbsnp+'.idx') is False:
        gatk4(refPath, 'gatk IndexFeatureFile -F data/'+dbsnp)

def NGSMainWESgetID(refPath, refGenome, dbsnp):
    with open('NGSMainWESid.txt', 'r') as f:
        for i in f.readlines():
            if i.startswith('CONTAINER') is False:
                mainID = i.split('        ')[0]
                mainBox = "sudo docker exec "+mainID+' '
                return mainID, mainBox

def NGSMainWES(refPath, scriptsPath, refGenome, dbsnp, cosmic):
    if os.path.exists(refPath+'/'+refGenome+'.fasta') is False or os.path.exists(refPath+'/'+dbsnp) is False or os.path.exists(refPath+'/'+cosmic) is False:
        print('>>>>>> deploying adgh456/ngs-main:wes in deamon...')
        cmd = 'sudo docker run -i -d -v /var/run/docker.sock:/var/run/docker.sock -v '+refPath+':/ref_data -v '+scriptsPath+':/scripts -i adgh456/ngs-main:wes_v2 '
        os.system(cmd)
        os.system('sudo docker ps -l > NGSMainWESid.txt')
        mainID, mainBox = NGSMainWESgetID(refPath, refGenome, dbsnp)
        print('>>>>>> download reference...')
        os.system(mainBox+"/bin/bash -c 'mv /box_ref_data/* /ref_data'")
        print('>>>>>> download scripts...')
        os.system(mainBox+"/bin/bash -c 'mv /box_scripts/* /scripts'")
        os.system('sudo docker rm -f '+mainID)
        ReferenceIndex(refPath, refGenome, dbsnp, cosmic)

    else:
        print('>>>>>> reference already exist')

###docker command
def afterQC(localPath, argv):
    cmd = 'sudo docker run --rm -v '+localPath+':/AfterQC/data -i adgh456/afterqc:afterqc_pypy2 '+argv
    os.system(cmd)
    os.system('gzip '+localPath+'/*good.fq')
    os.system('mv '+localPath+'/*html '+localPath+'/*good* ref_data')

def bwa(localPath, argv):
    cmd = 'sudo docker run --rm -v '+localPath+':/data -i adgh456/bwa:0.7.15 '+argv
    os.system(cmd)
    
def samtools(localPath, argv):
    cmd = 'sudo docker run --rm -v '+localPath+':/data -i adgh456/samtools:1.3.1 '+argv
    os.system(cmd)
    
def picard(localPath, argv):
    cmd = 'sudo docker run --rm -v '+localPath+':/data -i adgh456/picard:2.9.0-1-gf5b9f50-SNAPSHOT '+argv
    os.system(cmd)

def gatk3(localPath, argv):
    cmd = 'sudo docker run --rm -v '+localPath+':/data -i adgh456/gatk:gatk-v3.8-registered '+argv
    os.system(cmd)

def gatk4(localPath, argv):
    cmd = 'sudo docker run --rm -v '+localPath+':/gatk/data -i broadinstitute/gatk:4.0.4.0 '+argv
    os.system(cmd)

def gatk4CNV(localPath, argv):
    cmd = 'docker run --rm -v '+localPath+'/gatk-protected-1.0.0.0-alpha1.2.3/data -i adgh456/gatk:1.0.0.0-alpha1.2.3 java -jar /gatk-protected-1.0.0.0-alpha1.2.3/build/libs/gatk-protected.jar '+argv
    os.system(cmd)

def mutect(localPath, argv):
    cmd = 'docker run --rm -v '+localPath+':/opt/data -i docker.io/opengenomics/mutect '+argv
    os.system(cmd)
    
def muse(localPath, argv):
    cmd = 'docker run --rm -v '+localPath+':/data -i docker.io/opengenomics/muse '+argv
    os.system(cmd)

def somaticsniper(localPath, argv):
    cmd = 'docker run --rm -v '+localPath+':/opt/data -i quay.io/opengenomics/somatic-sniper:latest '+argv
    os.system(cmd)
    
def varscan2(localPath, argv):
    cmd = 'docker run --rm -v '+localPath+':/data -i quay.io/biocontainers/varscan:2.4.3--0 '+argv
    os.system(cmd)
    
def optitype(localPath, argv):
    cmd = 'docker run --rm -v '+localPath+':/data -i fred2/optitype '+argv
    os.system(cmd)

def meerkat(localPath, argv):
    cmd = 'docker run --rm -v '+localPath+':/app/ref_data -i adgh456/meerkat:v0.189 '+argv
    os.system(cmd)

def msisensor(localPath, argv):
    cmd = 'sudo docker run --rm -v '+localPath+':/data -i adgh456/msisensor:0.2 '+argv
    os.system(cmd)

###NGS analysis
def CheckDir(*path):
    p1 = path[0]+'/data'
    p2 = path[0]+'/ref_data'
    p3 = path[0]+'/data/'+path[1]
    p4 = p3+'/'+path[2]
    for i in (p1, p2, p3, p4):
        if os.path.exists(i) is False:
            print('>>>>>> '+i+' is not found in localhost, mkdir '+i+' automatically.')
            os.system('mkdir '+i)
        else:
            print('>>>>>> '+i+' is already exist.')
    
def MergeFastq(pathN, pathT, project, mainPath, storePath):
    print('>>>>>> merging [%s] fastq...' %(project))
    fastq_N1 = pathN+'/'+project+'N_1.fq.gz'
    fastq_T1 = pathT+'/'+project+'T_1.fq.gz'
    if os.path.exists(fastq_N1) is False or os.path.exists(fastq_T1) is False:
        ###pathN
        os.chdir(pathN)
        os.system('rm -f *json *html')
        os.system('ls -lt *fq.gz > merge.txt')
        with open('merge.txt', 'r') as f:
            fq_1 = ''
            fq_2 = ''
            n = []
            for i in f.readlines():
                n.append(i.split(re.findall(' [0-9][0-9]:[0-9][0-9]', i)[0]+' ')[1])
        for i in range(len(n)):
            if '1.clean.fq.gz' in n[i] or '1.fq.gz' in n[i]:
                fq_1 += n[i].strip()+' '
            elif '2.clean.fq.gz' in n[i] or '2.fq.gz' in n[i]:
                fq_2 += n[i].strip()+' '
        cmd_1 = 'cat '+fq_1+'> '+project+'N_1.fq.gz'
        cmd_2 = 'cat '+fq_2+'> '+project+'N_2.fq.gz'
        os.system(cmd_1)
        os.system(cmd_2)
        os.system('rm merge.txt')
        ###pathT
        os.chdir(pathT)
        os.system('rm -f *json *html')
        os.system('ls -lt *fq.gz > merge.txt')
        with open('merge.txt', 'r') as f:
            fq_1 = ''
            fq_2 = ''
            n = []
            for i in f.readlines():
                n.append(i.split(re.findall(' [0-9][0-9]:[0-9][0-9]', i)[0]+' ')[1])
        for i in range(len(n)):
            if '1.clean.fq.gz' in n[i] or '1.fq.gz' in n[i]:
                fq_1 += n[i].strip()+' '
            elif '2.clean.fq.gz' in n[i] or '2.fq.gz' in n[i]:
                fq_2 += n[i].strip()+' '
        cmd_1 = 'cat '+fq_1+'> '+project+'T_1.fq.gz'
        cmd_2 = 'cat '+fq_2+'> '+project+'T_2.fq.gz'
        os.system(cmd_1)
        os.system(cmd_2)
        os.system('rm merge.txt')
        
        os.chdir(mainPath)
                
def AfterQC(threshold, pathN, pathT, project, refPath, storePath):
    afterQC(pathN, 'pypy after.py -f0 -t0 -q '+str(threshold)+' -g /AfterQC/data -r /AfterQC/data -1 /AfterQC/data/'+project+'N_1.fq.gz -2 /AfterQC/data/'+project+'N_2.fq.gz')
    afterQC(pathT, 'pypy after.py -f0 -t0 -q '+str(threshold)+' -g /AfterQC/data -r /AfterQC/data -1 /AfterQC/data/'+project+'T_1.fq.gz -2 /AfterQC/data/'+project+'T_2.fq.gz')
    
def FastqtoSam(refPath, fasta, fqN1, fqN2, fqT1, fqT2, projectN, projectT):
    bwa(refPath, 'bwa aln -t 2 '+fasta+'.fasta '+fqN1+'.good.fq.gz > ref_data/'+fqN1+'.sai')
    bwa(refPath, 'bwa aln -t 2 '+fasta+'.fasta '+fqN2+'.good.fq.gz > ref_data/'+fqN2+'.sai')
    bwa(refPath, 'bwa sampe '+fasta+'.fasta '+fqN1+'.sai '+fqN2+'.sai '+fqN1+'.good.fq.gz '+fqN2+'.good.fq.gz > ref_data/'+projectN+'.sam')
    bwa(refPath, 'bwa aln -t 2 '+fasta+'.fasta '+fqT1+'.good.fq.gz > ref_data/'+fqT1+'.sai')
    bwa(refPath, 'bwa aln -t 2 '+fasta+'.fasta '+fqT2+'.good.fq.gz > ref_data/'+fqT2+'.sai')
    bwa(refPath, 'bwa sampe '+fasta+'.fasta '+fqT1+'.sai '+fqT2+'.sai '+fqT1+'.good.fq.gz '+fqT2+'.good.fq.gz > ref_data/'+projectT+'.sam')


def SamtoSortbam(refPath, projectN, projectT):
    picard(refPath, 'java -jar picard.jar SortSam I=data/'+projectN+'.sam O=data/'+projectN+'.sorted.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT')
    samtools(refPath, 'samtools index '+projectN+'.sorted.bam')
    picard(refPath, 'java -jar picard.jar SortSam I=data/'+projectT+'.sam O=data/'+projectT+'.sorted.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT')
    samtools(refPath, 'samtools index '+projectT+'.sorted.bam')
    
def MarkDuplicates(refPath, projectN, projectT):
    picard(refPath, 'java -jar picard.jar MarkDuplicates I=data/'+projectN+'.sorted.bam  O=data/'+projectN+'.sorted_dedup.bam REMOVE_DUPLICATES=true METRICS_FILE=rmdup.txt AS=true VALIDATION_STRINGENCY=LENIENT')
    picard(refPath, 'java -jar picard.jar AddOrReplaceReadGroups I=data/'+projectN+'.sorted_dedup.bam O=data/'+projectN+'.regroup.bam SORT_ORDER=coordinate RGID=illumina RGLB=library RGPL=iontorrent RGSM=sample RGPU=ATCG VALIDATION_STRINGENCY=LENIENT')
    samtools(refPath, 'samtools index '+projectN+'.regroup.bam')
    picard(refPath, 'java -jar picard.jar MarkDuplicates I=data/'+projectT+'.sorted.bam  O=data/'+projectT+'.sorted_dedup.bam REMOVE_DUPLICATES=true METRICS_FILE=rmdup.txt AS=true VALIDATION_STRINGENCY=LENIENT')
    picard(refPath, 'java -jar picard.jar AddOrReplaceReadGroups I=data/'+projectT+'.sorted_dedup.bam O=data/'+projectT+'.regroup.bam SORT_ORDER=coordinate RGID=illumina RGLB=library RGPL=iontorrent RGSM=sample RGPU=ATCG VALIDATION_STRINGENCY=LENIENT')
    samtools(refPath, 'samtools index '+projectT+'.regroup.bam')
    
def BaseRecalibrator(refPath, fasta, projectN, projectT, dbsnp):
    os.system('sudo apt install -y samtools')

    gatk4(refPath, 'gatk BaseRecalibrator -R data/'+fasta+'.fasta -I data/'+projectN+'.regroup.bam --known-sites data/'+dbsnp+' -O data/'+projectN+'.recal_data.table')
    gatk4(refPath, 'gatk ApplyBQSR -R data/'+fasta+'.fasta -I data/'+projectN+'.regroup.bam --bqsr-recal-file data/'+projectN+'.recal_data.table -O data/'+projectN+'.recal.bam')
    os.system('samtools view -H ref_data/'+projectN+".recal.bam | sed 's/SM:sample/SM:"+projectN+"/' | samtools reheader - ref_data/"+projectN+'.recal.bam > ref_data/'+projectN+'.recal_reheader.bam')
    samtools(refPath, 'samtools index '+projectN+'.recal_reheader.bam')
    
    gatk4(refPath, 'gatk BaseRecalibrator -R data/'+fasta+'.fasta -I data/'+projectT+'.regroup.bam --known-sites data/'+dbsnp+' -O data/'+projectT+'.recal_data.table')
    gatk4(refPath, 'gatk ApplyBQSR -R data/'+fasta+'.fasta -I data/'+projectT+'.regroup.bam --bqsr-recal-file data/'+projectT+'.recal_data.table -O data/'+projectT+'.recal.bam')
    os.system('samtools view -H ref_data/'+projectT+".recal.bam | sed 's/SM:sample/SM:"+projectT+"/' | samtools reheader - ref_data/"+projectT+'.recal.bam > ref_data/'+projectT+'.recal_reheader.bam')
    samtools(refPath, 'samtools index '+projectT+'.recal_reheader.bam')
    
def CreatePONforMutect2(refPath, fasta, projectN, projectT, project):
    gatk4(refPath, 'gatk Mutect2 -R data/'+fasta+'.fasta -I data/'+projectN+'.recal_reheader.bam -tumor '+projectN+' -O data/'+projectN+'_for_pon.vcf.gz')
    os.system('echo /gatk/data/'+projectN+'_for_pon.vcf.gz >> ref_data/normals_for_pon_vcf.args')

def Mutect2(refPath, fasta, projectN, projectT, project):
    gatk4(refPath, 'gatk Mutect2 -R data/'+fasta+'.fasta -I data/'+projectT+'.recal_reheader.bam -tumor '+projectT+' -I data/'+projectN+'.recal_reheader.bam -normal '+projectN+' -O data/'+project+'.somatic.vcf.gz')

def Mutect2_v3(refPath, fasta, projectN, projectT, project, cosmic, dbsnp):
    gatk3(refPath, 'gatk -T MuTect2 -R data/'+fasta+'.fasta -I:tumor data/'+projectT+'.recal_reheader.bam -I:normal data/'+projectN+'.recal_reheader.bam --cosmic data/'+cosmic+' --dbsnp data/'+dbsnp+' --contamination_fraction_to_filter 0.02 -o data/'+project+'.mutect2.vcf')

def MSIsensor(refPath, fasta, projectN, projectT, project, bed):
    if os.path.exists(refPath+'/'+fasta+'.microsatellites.list') is False:
        print('>>>>>> building '+fasta+'.microsatellites.list')
        msisensor(refPath, 'msisensor scan -d data/'+fasta+'.fasta -o data/'+fasta+'.microsatellites.list')
    msisensor(refPath, 'msisensor msi -d data/'+fasta+'.microsatellites.list -n data/'+projectN+'.recal_reheader.bam -t data/'+projectT+'.recal_reheader.bam -e data/'+bed+' -o data/msi.'+project)

def CheckVcf(refPath, subproject, project, storePath):
    if os.path.exists('ref_data/'+project+'.somatic.vcf.gz') is True:
        os.system('echo %s/%s >> good_report.txt' %(subproject, project))
        os.system('mv %s/%s* %s/msi.%s* %s/%s' %(refPath, project, refPath, project, storePath, project))
    elif os.path.exists('ref_data/'+project+'.somatic.vcf.gz') is False:
        os.system('echo %s/%s >> bad_report.txt' %(subproject, project))
        
def CreatePONforCNV(self):
    pass
    
def CNV(self):
    pass
    
def Phial(self):
    pass
