# ngs-main-wes
This is a Dockerized "somatic paired-WES(Normal-Tumor) hg19 pipeline" which is refer to GATK best practice.

## Prepare env:
1. git clone to your main project directory (ex. OSCC)
2. put "ngs-main-wes.py" and "NGStools.py" in your main project directory (ex. OSCC)
3. prepare your own raw data in sub project directory (ex. 66xWES)

## Usage: 
#### python3 ngs-main-wes.py
#### Customizable：
1. genome coordinate：ucsc.hg19.fasta
2. dbSNP：dbsnp_138.hg19.vcf
3. COSMIC：CosmicAllMutsHeaderSorted.vcf
＊ multi-thread：p = Pool(15)


- rawdataRename: rename raw data directory (ex. OC_631N to 631N)
- enterData: input patient ID (for parallel processing)


- CheckDir: mkdir ref_data (for reference and intermediate file), data (for storage)
- NGSMainWES: copy reference (ucsc.hg19.fasta, dbsnp_138.hg19.vcf, 'CosmicAllMutsHeaderSorted.vcf) to ref_data
- ReferenceIndex: build index for reference with samtools


- MergeFastq: merge multiple fastq into one (N.fq1, N.fq2, T.fq1, T.fq2)
- AfterQC: trimming raw fastq data with Afterqc (auto trimming), and produce a QC report
- FastqtoSam: align fastq reads to reference genome with BWA, and produce sam file
- SamtoSortbam: transform sam file (sequence alignment map) to bam file (binary alignment map) with picard, and build index with samtools
- MarkDuplicates: MarkDuplicates and AddOrReplaceReadGroups with picard, and build index with samtools
- BaseRecalibrator: BaseRecalibrator and ApplyBQSR with GATK4, and build index with samtools
- CreatePONforMutect2: call variant using Mutect2 (tumor-only mode) without dbSNP and COSMIC with GATK4, and record patient to build Panel of Normal in normals_for_pon_vcf.args
- Mutect2: call variant using Mutect2 (Paired: Normal and Tumor) without dbSNP and COSMIC with GATK4
- Mutect2_v3: call variant using mutect2 (Paired: Normal and Tumor) with dbSNP and COSMIC with GATK4
- CheckVcf: record result (good_report.txt or bad_report.txt)
- CreatePONforCNV: create PON of CNV
- CNV: call copy number variant with GATK4 1.0.0.0.alpha
- Phial: clinical FDA drug relevence annotation
