## Introduction 
This project contains one complex pipeline, one framework, and multiple tools. 

1: Pipeline 
- COVID_2nd_Pipeline

2: Framework
- COVID_Auto_Framework

3: Tools
- Combine Report
- BAM contamination check
- S3 backup (toolkit)

For details please check the sections below

### Pipeline: COVID_2nd_Pipeline

1: Brief instruction
- The major bash code are inherited from "CCAD 2nd pipeline"
   - All irrelevant code has been moved into the folder "IrrelavantCode"  
   - The existing code has been modified to support the new changes, including new designed interface, automation, parallel computing and the new logic of COVID project. 
   - The new implemented python scripts are used to 
     - implement the new logic of how to handle USU raw data (fastq) 
     - replace the old code (too specific for CGR data)  
     - make the specific steps be easy to link with the whole 2nd pipeline automatically.   

2: Key steps
- MergeBAM 
- ReconstructBAMDir
- RecalibrateBAMFile
- ConstructBAMRecaliPerManifest
- CoverageReport
- PreCallingQCReport
- BAMContaminationCheck

3: Specific Path: 
- Build results   
   - /data/COVID_WGS/lix33/Test/2ndpipeline/Build
- BAM_new_incoming  
   - /data/COVID_WGS/lix33/Test/2ndpipeline/Data/BAM_new_incoming
- BAM_old_backups (we do not expect this folder)  
   - /data/COVID_WGS/lix33/Test/2ndpipeline/Data/BAM_old_backups 
- BAM_reformatted original  
   - /data/COVID_WGS/lix33/Test/2ndpipeline/Data/BAM_reformatted/BAM_original 
- BAM_reformatted recalibrated 
   - /data/COVID_WGS/lix33/Test/2ndpipeline/Data/BAM_reformatted/BAM_recalibrated - It also contains two subfolders, including
    -  "Flag": the flags related to "BAM_reformatted recalibrated" (from CCAD old 2nd pipeline) will be recorded here    - "WGS": the recalibrated BAM will be recorded here- Cluster job log files (recalibrated, coverage, precalling)  - /data/COVID_WGS/lix33/Test/2ndpipeline/Data/secondary_buf/cluster_job_logs- Coverage report and Pre-calling report  - /data/COVID_WGS/lix33/Test/2ndpipeline/Data/secondary_buf/coverage_report- tmp folder that will be used to save the tmp file generated during the calculation - /data/COVID_WGS/lix33/Test/2ndpipeline/Data/tmp
  ### Framework: COVID_Auto_Framework
1: Brief instruction- The old "CCAD 2nd pipeline" can only be run manually- The COVID 2nd pipeline contains many different steps, and it really takes  time if we run it step by step manually. - Implement it by using the "flag" mechanism.- Run with crontab job & fully automatic2: Folder structure: - Keytable:   - /data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Keytable- Crontab job log  - /data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Log- Retrieved raw data (from S3)  - /data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Data/BAM/Batch- Analysis Log and flag
  - /data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Data/BAM/Batch/std_input_01_200_keytable/Flag  - /data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Data/BAM/Batch/std_input_01_200_keytable/Log
3: Different steps and key flags (including both existing old flags set by CCAD 2nd pipeline and new flags defined by "COVID_Auto_Framework")- MergeBAM  - Flags    - /data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Data/BAM/Batch/std_input_01_200_keytable/Flag/BuildManifest    - /data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Data/BAM/Batch/std_input_01_200_keytable/Flag/MergeSubject  - Jobscript (auto-generated):    - /data/COVID_WGS/lix33/Test/2ndpipeline/Build/tmp/std_input_01_200_keytable/Script/MergeSubject  - Mimic Manifest file (auto-generated)    - /data/COVID_WGS/lix33/Test/2ndpipeline/Build/tmp/std_input_01_200_keytable/Manifest_Mimic.csv  - Logs     - /data/COVID_WGS/lix33/Test/2ndpipeline/Build/tmp/std_input_01_200_keytable/Log/MergeSample
- ReconstructBAMDir  - Flag     - /data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Data/BAM/Batch/std_input_01_200_keytable/Flag/ReconstructBAMDir  - Log     - /data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Data/BAM/Batch/std_input_01_200_keytable/Log/ReconstructBAMDir- RecalibrateBAMFile  - Flag     - /data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Data/BAM/Batch/std_input_01_200_keytable/Flag/RecalibrateBAM     - /data/COVID_WGS/lix33/Test/2ndpipeline/Data/BAM_reformatted/BAM_recalibrated (*.inqueue file)     - /data/COVID_WGS/lix33/Test/2ndpipeline/Data/BAM_reformatted/BAM_recalibrated/Flag  - Log     - /data/COVID_WGS/lix33/Test/2ndpipeline/Data/secondary_buf/cluster_job_logs/20211020- ConstructBAMRecaliPerManifest  - Flag  - Log- CoverageReport  - Flag  - Log  - Report File - PreCallingQCReport
  - Flag  - Log  - Report File- BAMContaminationCheck  - Flag  - Log  - Report File
4: How to run it - Crontab job- Example ```@hourly . /etc/bashrc && source /home/lix33/lxwg/Git/IlluminaSequencingAnalysis/COVID_Auto_Framework/SourceCode/Run.sh /data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Keytable/std_input/sub_keytable/std_input_01_200_keytable.csv >> /data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Log/$(date +\%Y\%m).txt 2>&1
```
### Tool: Combine Report1: Support combine 3 different types of QC report from different builds 2: 3 different types report includes: - Coverage report- Pre-calling QC report- BAM Contamination report
### Tool: BAM contamination check 1: the old code was implemented by perl, which is hard to be concatenated by other pipelines 2: Re-implemented it by using python and added a bunch of new features. 3: Fully automatic- Can be run separately - Can be run by using COVID_Auto_Framework. 
### Tool kit: Biowlf&S3 Backup&Retrieve1: Due to the storage limitation, we need to backup data from biowulf to S3 and retrieve data from S3 to biowulf. 2: Implemented several tools to do this jobs- Tool 1: Backup COVID 2nd results from biowulf to S3 - Tool 2: Retrieve BAM from S3 to biowulf - Tool 3: Backup COVID primary results from biowulf to S3
