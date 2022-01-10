### Introduction
This is a project to prepare sliced BAM file for Cenk's team. 

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
*/30 * * * * . /etc/bashrc && cd /home/lix33/lxwg/Git/IlluminaSequencingAnalysis/Ad-hoc/Cenk/SourceCode && python3 /home/lix33/lxwg/Git/IlluminaSequencingAnalysis/Ad-hoc/Cenk/SourceCode/GetFileFromS3.py >> /data/COVID_ADHOC/Sequencing/COVID_WGS/Data/USU_First_Batch_low_std_96_361_18/debug/Log/sum.log 2>&1
```


#### 2: Slice Region from BAM
```
bash SGFBam.wrapper.sh
```

#### 3: Update Excel
```
python3 UpdateExcel.py
```

### Related Development History
Please check the ticket below:
https://github.com/NCI-CGR/IlluminaSequencingAnalysis/issues/51
