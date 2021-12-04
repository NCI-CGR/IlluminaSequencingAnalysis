10-25 to 11-25 Go

### 1: For data analysis 
-------------------------
1: The analysis of all subjects for the first and the second batches deliveried from USU are all set. 
   (1) All recalibreated BAMs (target BAM) have been generated 
   (2) All types types of QC report for each subjects have been generated 

2: All raw fastq files  (31 TB) have been moved from biowulf to S3. 


3: Some finished BAM are still located in biowulf   
   (1) back up to S3 
   (2) since wen may use it for downstream analysis, will discuss with her when is the best time to remove it from biowulf 

-------------------------
*******************************

### 2: Trouble shooting 

-------------------------
1: Due to the code error (one additional space), we failed to back up the recalibrated BAM for 96 low input sample before. 
   * That's why we have to do it again. 
   
2: bizzard performance issue of GATK 
  (1) Some jobs is 3 or 10 times slower than other. 
  (2) Same command line but different system output.
    
 
3: mess to re-run the failed and missing and new subject together. 
  (1) Out corret pipeline do not support run differnet batch togethr (result will be mixed together). 
  However, 
  (1) For finished batch,they may miss some result  (need to regenerated, only run part of pipeline)
  (2) For current running batch: most of them are done, and only one to subject are still running 
  (3) We may still have one more batch (keytable) need to analysis 
  
  Since we have enough storage right now, we plan to run these together to save time. 
  
  Manually make some script to make these parallel run (multiple differnet batches) be possible. 
  
  Will update current pipeline to support this feature in the future. 

--------------------------
*****************************


### 3: Pacbio WGS SV pipeline 

-------------------------    
1: Install some latest version of some python packages to run it 
2: Check the code 
3: Discuss with the Pac bio tech support to solve some issues of running the pipeline. 
   (1) we should follow the rule to name the folder and BAM file, otherwise pipelien cannot regonize it. 
   (2) The rule is wrote by a complex regular expression  (hard to be read)

### 4: Debug the additional issues. 

------------------------------
*****************************


### 5: Others 

--------------------------
1: RNA QC Pipelin: make it support multiple project ID in one flowcell 
2: Some permission issue fro kristie. 
3: Data transfer in GCP


---------------------------
### 6: Remainign works 
1: Submit code, finished docs. 
2: Move data to S3
3: Waiting Clifton finishe the data transfer 
4: Write code to retrieve data from GCP to biowulf (back up fro biowulf to GCP)
5: Analysis the new delivered data from clifton. 
6: Continue Pacbio WGS SV pipeline.


---------------------------
### 7: Question
1: Which research groups are the major groups that will use COIVD USU data for their research project?
2: Which research project they plan to do respectively? 
3: If it is possible for us to partially attent it (looks like interesting)? 
