# -*- coding: utf-8 -*-
"""
Created on Thu May 10 11:07:44 2018

@author: Shang-Hung, Shih
"""

import os
import sys
import NGStools
import time
from multiprocessing import Pool

class somaticWES:
    
    """
    ===||| This is a Dockerize "somatic paired-WES(Normal-Tumor) pipeline" which is refer to GATK best practice. |||===
    """
    
    def __init__(self, mainPath, dataPath, project):
        
        ###dir to modify
        self.mainPath = mainPath
        self.dataPath = dataPath
        self.project = project
        
        ### reference   ###ftp://ftp.broadinstitute.org/bundle/
        self.refGenome = 'ucsc.hg19'
        self.dbsnp = 'dbsnp_138.hg19.vcf'
        self.cosmic = 'CosmicAllMutsHeaderSorted.vcf'
        
        ###dir no need to modify
        self.rawdataPath = mainPath+'/'+dataPath
        self.storePath = mainPath+'/data/'+dataPath
        self.refPath = mainPath+'/ref_data'
        self.scriptsPath = mainPath+'/scripts'
        self.project_N = project+'N'
        self.project_T = project+'T'
        self.fq1_name_N = project+'N_1'
        self.fq2_name_N = project+'N_2'
        self.fq1_name_T = project+'T_1'
        self.fq2_name_T = project+'T_2'
        
        ###container dir
        self.boxMainPath = '/'
        self.boxRefPath = '/ref_data'
        self.boxDataPath = '/data'
        self.boxScriptsPath = '/scripts'

###help(somaticWES)

######rawdataPreprocess
NGStools.rawdataRename()

######enterData
work = NGStools.enterData()

######multiprocessing
def work_log(work_data):
    startTime = time.time()

    p1 = somaticWES(os.getcwd(), work_data[0], work_data[1])
    print('>>>>>> start: subProject=%s , patient=%s.' %(work_data[0], work_data[1]))
    NGStools.CheckDir(p1.mainPath, p1.dataPath, p1.project)
    #NGStools.NGSMainWES(p1.refPath, p1.scriptsPath, p1.refGenome, p1.dbsnp)
    NGStools.ReferenceIndex(p1.refPath, p1.refGenome, p1.dbsnp, p1.cosmic)

    NGStools.MergeFastq(p1.rawdataPath+'/'+p1.project_N, p1.rawdataPath+'/'+p1.project_T, p1.project, p1.mainPath, p1.storePath)
    NGStools.AfterQC(30, p1.rawdataPath+'/'+p1.project_N, p1.rawdataPath+'/'+p1.project_T, p1.project, p1.refPath, p1.storePath)
    NGStools.FastqtoSam(p1.refPath, p1.refGenome, p1.fq1_name_N, p1.fq2_name_N, p1.fq1_name_T, p1.fq2_name_T, p1.project_N, p1.project_T)
    NGStools.SamtoSortbam(p1.refPath, p1.project_N, p1.project_T)
    NGStools.MarkDuplicates(p1.refPath, p1.project_N, p1.project_T)
    NGStools.BaseRecalibrator(p1.refPath, p1.refGenome, p1.project_N, p1.project_T, p1.dbsnp)
    NGStools.CreatePONforMutect2(p1.refPath, p1.refGenome, p1.project_N, p1.project_T, p1.project)
    NGStools.Mutect2(p1.refPath, p1.refGenome, p1.project_N, p1.project_T, p1.project)
    NGStools.Mutect2_v3(p1.refPath, p1.refGenome, p1.project_N, p1.project_T, p1.project, p1.cosmic, p1.dbsnp)
    NGStools.CheckVcf(p1.refPath, work_data[0], p1.project, p1.storePath)

    endTime = time.time()
    print('>>>>>> TotalTimeUse for subproject=%s , patient=%s : %s sec' %(work_data[0], work_data[1], endTime-startTime))


def pool_handler():
    p = Pool(15)
    p.map(work_log, work)

if __name__ == '__main__':
    pool_handler()



