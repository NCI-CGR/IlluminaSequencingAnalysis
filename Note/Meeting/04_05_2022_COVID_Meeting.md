## Discussion
### Original Email
Hi all,

1) Need to start addressing the WGS data that has been generated to date, 

2) review QC metrics for what we have been able to analyze on Biowulf/CCAD for some, and reconcile files received with samples sent.


3) We’ll follow up this first meeting with a larger one that will start to address data that is currently received and stored on GCP

4) Discuss potential use of Nvidia/parabricks for analysis, 

5) eventually loop in others to discuss the larger analytical plan post variant calling (including potential use for informing/improving imputation for some of our Hispanic GWAS efforts.

### For Question 1
Question 1: Need to start addressing the WGS data that has been generated to date,

#### 1: How many samples we alreay received

We currently received 4 + 1 batches of data

(1) Batch 1
  * a) Standard Input WGS	
  * b) Total number of Sample: 379
  * c) Received Date: July, 2021
  * d) Saved Location: Object Oriented Storage System，DCEG_COVID_WGS
    * Use obj_ls to check it:  
      ```
      obj_ls -v DCEG_COVID_WGS -h -m "*fastq.gz"| less -SN
      Notice: contains other fastq files, need to check the output to know which part belongs to current batch
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
  * d) Saved Location: GCP, pI3.COVNET.batch3.N638.6topoffs.keytable.dec.16.202
  * e) Analysis Status: Nil
  * f) Keytable File
  
![image](https://user-images.githubusercontent.com/11053933/161799717-9800c274-a4b8-4ff1-98e4-aa33cc9f196b.png)


(4) Batch 5 (newly added)
   * a) It looks like we just got batch 5 on April 7,8, 2022. 
   * b) The name of Batch: pI3.COVNET.NCI.batch4.N1512.01/
   ![image](https://user-images.githubusercontent.com/11053933/163228922-7b8885e9-7f17-4570-bad9-529c9b060791.png)

   * b) Details
      *  i) Number of sample: 1512
   * c) Question
      * i) one fastq outside sub folders: dceg-covnet-wgs-useast1/I3-18752_S14_L001_R1_001.fastq.gz
      * ii)missing keytable csv file

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
Question 2: review QC metrics for what we have been able to analyze on Biowulf/CCAD for some, and reconcile files received with samples sent.
#### 1: How many types of QC reports do we have currently.
We typically has two different phases of COVID upstream analysis, including Custmized QC and 2nd pipeline 

![image](https://user-images.githubusercontent.com/11053933/163254417-ed95ca78-fb79-4e09-9d02-4af3ba52a67c.png)


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
       * Comparing with primary pipeline: contains additional statistic info, including 1x, 5x, 50x
    ```
    obj_ls -v DCEG_COVID_WGS -h -m "*coverage_report_*.txt"| less -SN
    ```
    ![image](https://user-images.githubusercontent.com/11053933/162490379-2e2a36d6-0110-4cfc-aff6-937f53e4063f.png)
    
    *  For pre calling QC Report (txt)
       * overlap some info with primary pipeline QC, but contains many new info.
    ```
    obj_ls -v DCEG_COVID_WGS -h -m "*pre_calling_qc_report_*.txt"| less -SN
    ```
    ![image](https://user-images.githubusercontent.com/11053933/162490652-fa279bf6-2ed5-4c21-b44d-1df503965960.png)


### For Question 3 
Question 3: We’ll follow up this first meeting with a larger one that will start to address data that is currently received and stored on GCP

#### Current Status
1. Currently both Batch 3 and Batch 4 were stored in GCP
2. For Batch 3 (Completed)
   * a) Std Input WGS	
   * b) Total number of Sample: 94
   * c) Received Date: November, 2021
   * d) Saved Location: GCP,  pI8.covnet.N94.fastq.00/
      * We downloaded these samples from GCP to biowulf and did analysis on Biowulf (NOT on GCP).
      * GCP location
    ```
    https://console.cloud.google.com/storage/browser/dceg-covnet-wgs-useast1;tab=objects?forceOnBucketsSortingFiltering=false&project=nih-nci-dceg-cgr&prefix=&forceOnObjectsSortingFiltering=false&pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))
    ```
   * e) Analysis Status: Completed
   * f) Keytable File

![image](https://user-images.githubusercontent.com/11053933/161797533-4c6df4ff-dfb8-4a9c-95ff-32fdf4bc2dd4.png)

3. We are talking about the data from Batch 4 (This is our target)
   * a) Std Input WGS
   * b) Total number of Sample: 638
   * c) Received Date: January, 2022
   * d) Saved Location: GCP, pI3.COVNET.batch3.N638.6topoffs.keytable.dec.16.202
   * e) Analysis Status: Nil
   * f) Keytable File

![image](https://user-images.githubusercontent.com/11053933/161799717-9800c274-a4b8-4ff1-98e4-aa33cc9f196b.png)

4. Similar thing for Batch 5 (pI3.COVNET.NCI.batch4.N1512.01)

### For Question 4
Question 4: Discuss potential use of Nvidia/parabricks for analysis,
#### Some questions/discussions
1. Regarding parabricks
   * Optimize the existing popular bioinfomatics tools by fully taking advantage of GPU cards.
   * Any benchmark comparsion? 
   * Based on our CGR in-house experience, how about its performance for real case.
   
2. Regarding using Nvidia/parabricks for analysis
   * Any similar jobs has been done? (GenScan)?
   * Any small demo code availabe for testing? 
  
3. Regarding using GCP for the implementation of the whole calculation logic 
   * If bioinfomatics tools are all availabe in google cloud? 
   * Do we need to deploy our our env in cluster. Or we can reused some existed & well-buit env in cloud (e.g. conda).
   * If the cloud has some issues? 
      * e.g. the zombie job in HPC? Any experiences to detect and re-submit jobs for this case?   
   
4. Cost 
    * Due to the cost of cloud storage, do we have the plan to backup file and delete files periodly in cloud?
    * Any other concerns/strategies to avoid the unnecessary cost in cloud? 
  
5. Reference
   * Parabricks: https://docs.nvidia.com/clara/parabricks/v3.5/text/software_overview.html

### For Question 5
Question 5: eventually loop in others to discuss the larger analytical plan post variant calling (including potential use for informing/improving imputation for some of our Hispanic GWAS efforts
#### Discussions
1. Any comments from Jia and Wen?
2. Do we need to make our Pacbio SV pipeline suppot cloud eventually? 
