## Discussion
Hi all,

1) Need to start addressing the WGS data that has been generated to date, 

2) review QC metrics for what we have been able to analyze on Biowulf/CCAD for some, and reconcile files received with samples sent.


3) We’ll follow up this first meeting with a larger one that will start to address data that is currently received and stored on GCP

4) Discuss potential use of Nvidia/parabricks for analysis, 

5) eventually loop in others to discuss the larger analytical plan post variant calling (including potential use for informing/improving imputation for some of our Hispanic GWAS efforts.

### For Question 1
Question 1: Need to start addressing the WGS data that has been generated to date,

#### 1: How many samples we alreay received

We currently received 4 batches of data

(1) Batch 1
  * a) Standard Input WGS	
  * b) Total number of Sample: 379
  * c) Received Date: July, 2021
  * d) Saved Location: Object Oriented Storage System，DCEG_COVID_WGS
    * Use obj_ls to check it:  
      ```
      obj_ls -v DCEG_COVID_WGS -h -m "*fastq.gz"| less -SN
      Notice: contains other fastq files, need check the code which part belongs to current batch
      ```
    * Screenshot example
    ![image](https://user-images.githubusercontent.com/11053933/162484986-436a5802-a960-4daa-a440-5b89d02cf9bf.png)

  * e) Analysis Status: Completed
  * f) Details 
    * i) 18 samples with 15 topoffs 
    * ii) 361 regular samples (no topoffs)
  * g) Keytable File
![image](https://user-images.githubusercontent.com/11053933/161784598-adbf4a53-09ae-4096-8d8f-84972d3b9bea.png)
	
(2) Batch 2
  * a) Low Input WGS	
  * b) Total number of Sample: 96
  * c) Received Date: July, 2021
  * d) Saved Location: Object Oriented Storage System，DCEG_COVID_WGS
    * Use obj_ls to check it:
    ```
    obj_ls -v DCEG_COVID_WGS -h -m "*PI8.Tagmentation.N96.00*fastq.gz"| less -SN
    Notice: contains other fastq files, need check the code which part belongs to current batch
    ```
    * Screenshot example
    
    ![image](https://user-images.githubusercontent.com/11053933/162485635-251b15d0-9464-45c3-8d53-fff7527e71d0.png)

    
  * e) Analysis Status: Completed
  * f) Keytable File
![image](https://user-images.githubusercontent.com/11053933/161785628-528e0e72-baf8-4a5b-831a-4abf7d3c8512.png)

(3) Batch 3
  * a) Std Input WGS	
  * b) Total number of Sample: 94
  * c) Received Date: November, 2021
  * d) Saved Location: GCP,  pI8.covnet.N94.fastq.00/
    * We downloaded these samples from GCP to biowulf and did analysis on Biowulf (NOT on GCP).
    * GCP location
    ```
    https://console.cloud.google.com/storage/browser/dceg-covnet-wgs-useast1;tab=objects?forceOnBucketsSortingFiltering=false&project=nih-nci-dceg-cgr&prefix=&forceOnObjectsSortingFiltering=false
    ```
  * e) Analysis Status: Completed
  * f) Keytable File

![image](https://user-images.githubusercontent.com/11053933/161797533-4c6df4ff-dfb8-4a9c-95ff-32fdf4bc2dd4.png)

(4) Batch 4
  * a) Std Input WGS	
  * b) Total number of Sample: 638
  * c) Received Date: January, 2022
  * d) Saved Location: GCP,   pI3.COVNET.batch3.N638.6topoffs.keytable.dec.16.202
  * e) Analysis Status: Nil
  * f) Keytable File
  
![image](https://user-images.githubusercontent.com/11053933/161799717-9800c274-a4b8-4ff1-98e4-aa33cc9f196b.png)


#### 2: The working directory in biowulf
(1) All works are on biowulf

(2) Details
  * Working Directory: /data/COVID_WGS/lix33
  * Ad-hoc requirement Directory: /data/COVID_ADHOC/Sequencing/COVID_WGS
    * Dr Cenk's research lab. 
  * Include Three parts:
    * Customized QC (primary pipeline, fastq -> BAM)
    * 2nd pipeline (merge sample, bam contamination check, etc)
    * COVID Automation framework (make all steps be processed automatically)
  * Git Repo
    * Custmized QC: https://github.com/NCI-CGR/IlluminaSequencingAnalysis/tree/main/CustomizedQC
    * 2nd pipeline: https://github.com/NCI-CGR/IlluminaSequencingAnalysis/tree/main/COVID_2nd_Pipeline
    * COVID Auto Framework: https://github.com/NCI-CGR/IlluminaSequencingAnalysis/tree/main/COVID_Auto_Framework
  * How to understand Customized QC and 2nd pipeline (documentation) 
    * https://github.com/NCI-CGR/IlluminaSequencingAnalysis/blob/main/Note/Project/COVID/COVID_2nd_pipeline.md
    * https://github.com/NCI-CGR/IlluminaSequencingAnalysis/blob/main/Note/Project/COVID/COVID_2nd_pipeline_addendum.md
    * https://github.com/NCI-CGR/IlluminaSequencingAnalysis/blob/main/Note/Project/COVID/COVID_primary_pipeline.md

#### 3: The object oriented system that was used for this project
(1) Vault: DCEG_COVID_WGS

![image](https://user-images.githubusercontent.com/11053933/161803128-a0e4bb42-c69d-40f7-8e9e-b28e124f3e81.png)



### For Question 2
#### 1: How many types of QC reports do we have currently.
We typically has two different phases of COVID upstream analysis, including Custmized QC and 2nd pipeline 

(1) Custmized QC
  * Classify the raw fastq files with the same flwocell name together, and generate differnet flowcells.
  * Generate ONE QC report FOR EACH FLOWCELL (the same as our primary pipeline) 
  * Please check the attachment for details (Customzied_QC_Report.zip)
 
 [QCReport-BWA-v38.csv](https://github.com/NCI-CGR/IlluminaSequencingAnalysis/files/8420219/QCReport-BWA-v38.csv)

  * For all QC report please check MS team
 
 ![image](https://user-images.githubusercontent.com/11053933/161804939-d18e86eb-9d26-4167-8830-70f131f67bca.png)

 

(2) 2nd pipeline
  * Generate THREE QC reports, including 
    *  Coverage QC Report (txt)
    *  Pre Calling QC Report (txt)
    *  BAM Contamination Report (csv)
  * Please check the attachment for details (2nd_QC_Report.zip)
 
[2nd_QC_Report.zip](https://github.com/NCI-CGR/IlluminaSequencingAnalysis/files/8420210/2nd_QC_Report.zip)

  * All 2nd QC reports (3 types) have been saved in object oriented storage system.
    * For BAM contamination QC report (CSV)
    ```
    obj_ls -v DCEG_COVID_WGS -h -m "*SumContaminationReport.csv"| less -SN
    ```
    ![image](https://user-images.githubusercontent.com/11053933/162489838-87de18d2-307c-4681-bf6f-493655a0df16.png)
    
    *  For Coverage QC Report (txt)
      * Comparing with primary pipeline: has additional statistic info of 1x, 5x, 50x
    ```
    obj_ls -v DCEG_COVID_WGS -h -m "*coverage_report_*.txt"| less -SN
    ```
    ![image](https://user-images.githubusercontent.com/11053933/162490379-2e2a36d6-0110-4cfc-aff6-937f53e4063f.png)
    
    *  For pre calling QC Report (txt)
    ```
    obj_ls -v DCEG_COVID_WGS -h -m "*pre_calling_qc_report_*.txt"| less -SN
    ```
    ![image](https://user-images.githubusercontent.com/11053933/162490652-fa279bf6-2ed5-4c21-b44d-1df503965960.png)


### For Question 3 
We are talking about the data from Batch 4
  * a) Std Input WGS
  * b) Total number of Sample: 638
  * c) Received Date: January, 2022
  * d) Saved Location: GCP, pI8.covnet.N94.fastq.00/
  * e) Analysis Status: Nil
  * f) Keytable File

![image](https://user-images.githubusercontent.com/11053933/161799717-9800c274-a4b8-4ff1-98e4-aa33cc9f196b.png)

### For Question 4
Some questions 
1. Any similar jobs has been done? (GenScan)?
2. Any small demo code availabe for testing? 
3. If bioinfomatics tools are availabe in google cloud? 
4. If the cloud has some issues as the zombie job in HPC? Any experiences to detect and re-submit jobs in this case?  
5. Cost 
    * 1) Due to the cost of cloud storage, do we have the plan to backup file and delete files periodly in cloud?
    * 2) Any other concern the things need to pay attention to avoid the unnecessary cost in cloud? 

### For Question 5
Any comments from Jia and Wen?
