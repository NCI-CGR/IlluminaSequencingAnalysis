This is the project for generating customized QC and run in Biowulf

### Purpose 
1: Conduct alignment for the given non-standard (diffeent from our CGR flowcell) sample data (PE reads)

2: Generate QC report
- Similar style of report as our primary pipeline.

3: The BAM file and QC report will be used for the study of Covid 19

### Highlights

1: Everything should run in Biowulf (NOT CGR cluster)
- We need to use slurm rather than sge to submit jobs
- The running env has been changed
 - Module file system 

2: Should support high-through output dataset. 

3: Everything should be automatical. 

4: Should support Email Notification. 

### Implementation 
1: Fully reuse the existing code in primary pipeline 

2: Use the Flag system for automation 
- Can give us the maximum flexibility

3: Use python for the pipeline architecture. 

### Working directory (**Biowulf**)
1: Testing code
- /scratch/lix33/lxwg/SourceCode/CommonTools/CustomizedQC

2: Testing data
- /scratch/lix33/lxwg/Data/ad-hoc/WGS

3: Running results: 
- /scratch/lix33/lxwg/Data/ad-hoc/WGS/ProcessedData

### How to run the pipeline
1： Command line

```
Step 1 (Optional): Create seperate folder and put all of target delivered folders inside -- one time
   - This is optional but highly recommanded, since 
      - some new data may come in when you run the command line.
      - the file may not in sub-folder
      
Step 2: Reconstruct data  -- one time
python3 /home/lix33/lxwg/Git/IlluminaSequencingAnalysis/CustomizedQC/SourceCode/DataReconstruct.py /data/COVID_WGS/primary_analysis/Data/07_02_2021 /data/COVID_WGS/primary_analysis/COVID19/06_30_2021/ProcessedData

Step 3: Run Customized QC python sourcecode  -- crontab job
python3 /home/lix33/lxwg/Git/IlluminaSequencingAnalysis/CustomizedQC/SourceCode/CustomizedQC.py /data/COVID_WGS/primary_analysis/COVID19/06_30_2021/ProcessedData

Step 4 (Optional): Move data from biowulf to S3 (object storage system)
python3 /home/lix33/lxwg/Git/IlluminaSequencingAnalysis/CustomizedQC/SourceCode/ObjectStorage/Backup2S3.py
Also check: /home/lix33/lxwg/Test/slurm/object_storage/job.sh

Step 5: merge qc report together
python3 /home/lix33/lxwg/Git/IlluminaSequencingAnalysis/CustomizedQC/SourceCode/MergeQCReport.py
```

2: Output: 
- Step 1
   - softlink to the existing folder or fastq files
   - example: /data/COVID_WGS/primary_analysis/Data/07_26_2021

- Step 2
   - Flowcells which parsed from the existing fastq files  
   - example: /data/COVID_WGS/primary_analysis/COVID19/06_30_2021/ProcessedData
 
- Step 3
   - BAM files, QC report will be in each flowcell folder
   - example: /data/COVID_WGS/primary_analysis/COVID19/06_30_2021/ProcessedData
 
- Step 4
   - In vault: DCEG_COVID_WGS
   - you can use obj_ls and obj_df to check it

- Step 5
   - Once everything is all set, merge qc report together 
   - example: /data/COVID_WGS/primary_analysis/COVID19/06_30_2021/QCReport
 
### Naming rules
- Aligner name and reference verion will be appended as the suffix in final QC report
- Considering we may need to supprt multiple aligner and reference verison in the future, 2-level-sub-folder is located in all possible places, including 
 - ./Flag/BWA/v38
 - ./Log/BWA/v38
 - ./reports/BWA/v38
 - ./Sample/Sample_SC742232_TCACAGCA-AATGCCTC_L001/Flag/BWA/v38 (sample level flag)
 - ./Sample/Sample_SC742246_TGTTCGAG-GAGATACG_L001/Flag/BWA/v38 (sample level flag)
 - ./Sample/Sample_SC742247_AGAGGTTG-ACGATGAC_L001/Flag/BWA/v38 (sample level flag)
 - ./BAM/BWA/v38
 - ./ttemp/BWA/v38
