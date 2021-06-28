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
python3 ./CustomizedQC.py /scratch/lix33/lxwg/Data/ad-hoc/WGS
```

2: Output: 
- CustomizedQC will auto create a folder named "ProcessedData" in the folder of input file
 - All outputs will locate in the folder "ProcessedData"
- For example:
 - Final results: /scratch/lix33/lxwg/Data/ad-hoc/WGS/ProcessedData
 - BAM file: /scratch/lix33/lxwg/Data/ad-hoc/WGS/ProcessedData/BAM/BWA/v38
 - QC report: /scratch/lix33/lxwg/Data/ad-hoc/WGS/ProcessedData/QCReport-BWA-v38.csv
 
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
