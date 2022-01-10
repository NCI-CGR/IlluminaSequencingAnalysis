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

#### 2: Slice Region from BAM
```
bash SGFBam.wrapper.sh
```

#### 3: Update Excel
```
python3 UpdateExcel.py
```
