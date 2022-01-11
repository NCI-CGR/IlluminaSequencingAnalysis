### Introduction (Background)
This is a project to prepare sliced BAM file for Cenk's team. 
  * Test a computational tool, Immunotyper, by using the large scale COVID dataset.
  * Original email
  
     *I am a senior investigator in the Cancer Data Science Lab at CCR. My lab has been developing Immunotyper, a computational tool to genotype and identify the alleleic composition of the  IGH region at the germ line. We are especially interested in the IGH-V genes to understand the germ line contribution to the antibody repertoire. We have been working on a subset of the NIAID COVID Consortium data with the goal of identifying distinct germline IGH alleles associated with type 1 interferon autoantibodies. We would like to expand the cohort to increase our statistical significance and I understand that the COVIDNET cohort is fairly large, but as per Dr. Holland’s email below these are typically not severe cases and no autoantibody analysis has been performed on them. Nevertheless they could still be quite useful as controls. In particular I would be very much interested in checking out if any of the autoantibody associated alleles are commonly found in this cohort especially among the older patients. If you would open to collaboration towards a paper we are putting together, primarily describing our computational method but also reporting on any association we may find between distinct alleles and autoantibodies it would be great if we could get access to the WGS data from this cohort.*

### Our initial schedule
```
(1) Amy, Vibha, Lisa, can you pull together a manifest of the first sets of samples that we sent to USU (the 96 low input + 384 standard)? 
   I know Clifton’s group is starting to deliver the set of 768 (first set, not second), but I believe that has just started so I don’t believe we need to consider those. 
(2) Once that manifest is prepared, Xin, can you slice out regions based on the attached .bed file from the final merged .bam files for these samples? 
   You’ll also need to add a file name to the manifest so that they can map files:subjects.
(3) Lisa, Meredith, you can see the requested phenotype info they would like to have below. Let me know what you can/cannot provide.
(4) Nathan, no rush, but need to set up a Globus transfer for this once Xin has files prepared.
```

### How to run it 
```
cd SourceCode
```

#### 1: Retrieve BAM from S3
```
load python/3.7
python3 GetFileFromS3.py
```
**Notice**: Since there is an issue in S3, we should always keep in mind three things below:
1. Never run to many jobs (e.g. > 100 jobs) simultaneously when you retrieving data from S3 to biowulf.
2. In the case of retrieving data from S3 to biowulf, always make sure each single job can be finished within 24 hours. 
3. In the case of retrieving data from S3 to biowulf, always make sure the number of concurrent transfers is around 20. 

Therefore, we need to deploy a crontab job to retreive data from S3 to biowulf. 
How to deploy crontab: 
```
cd /home/lix33/lxwg/Git/IlluminaSequencingAnalysis/Ad-hoc/Cenk/SourceCode && */30 * * * * . ~/.bashrc && cd /home/lix33/lxwg/Git/IlluminaSequencingAnalysis/Ad-hoc/Cenk/SourceCode && python3 /home/lix33/lxwg/Git/IlluminaSequencingAnalysis/Ad-hoc/Cenk/SourceCode/GetFileFromS3.py >> /data/COVID_ADHOC/Sequencing/COVID_WGS/Data/USU_First_Batch_low_std_96_361_18/debug/Log/sum.log 2>&1
```

#### 2: Slice Region from BAM
```
cd /home/lix33/lxwg/Git/IlluminaSequencingAnalysis/Ad-hoc/Cenk/SourceCode && bash SGFBam.wrapper.sh
```

#### 3: Update Excel
```
cd /home/lix33/lxwg/Git/IlluminaSequencingAnalysis/Ad-hoc/Cenk/SourceCode && python3 UpdateExcel.py
```

### Related Development History
1. Please check the github ticket below:
 * https://github.com/NCI-CGR/IlluminaSequencingAnalysis/issues/51
2. Please check the frogbug tickets below:
 * https://cgr-bugz.nci.nih.gov/f/cases/30077#BugEvent.372653

### Some Resources
1: Manifest file (from Lisa) and BED file (from Cenk)
  * Manifest file: COVIDwgs_Jan2022_NIAIDphenotypesLM.xlsx
    * Contains no duplicate sample
    * Contains 5 empty samples (cannot find their related BAM file from the dataset provided by USU)   
  * BED file: igh_orphons_GRCh38.bed
  * Check **Files** Dir

2: File system and working directory
  * We have been granted a new shared folder (file system) in biowulf to handle these ad-hoc request.
    * Name: COVID_ADHOC
    * Path: /data/COVID_ADHOC
    * Limitations: The maximum capacity of this folder is 30TB
    
  * Our working directory (All sequencing related project should go there):
    * /data/COVID_ADHOC/Sequencing

3: Retrieved Data from S3 to Biowulf
  * /data/COVID_ADHOC/Sequencing/COVID_WGS/Data/USU_First_Batch_low_std_96_361_18/debug

4: Results for Cenk's request
  * /data/COVID_ADHOC/Sequencing/COVID_WGS/ad-hoc/Cenk

