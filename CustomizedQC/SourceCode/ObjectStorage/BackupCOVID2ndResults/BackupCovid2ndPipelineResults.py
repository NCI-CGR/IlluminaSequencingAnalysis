'''
Input: 
1: Keytable Name
2: Array: bash dir list

Notice: The folder structure is somethings like below:

>>>>>>>>
Root Path in S3: 
 lix33/COVID19/UpstreamAnalysis/PostPrimaryRun/Data/10_11_2021/low_input_01_96_keytable
 
 -------We need to copy the file below from biowulf to s3
 /data/COVID_WGS/lix33/Test/2ndpipeline/Build
 
 /data/COVID_WGS/lix33/Test/2ndpipeline/Data/BAM_reformatted/BAM_original/WGS
 /data/COVID_WGS/lix33/Test/2ndpipeline/Data/BAM_reformatted/BAM_recalibrated/WGS
 
 /data/COVID_WGS/lix33/Test/2ndpipeline/Data/secondary_buf/cluster_job_logs
 /data/COVID_WGS/lix33/Test/2ndpipeline/Data/secondary_buf/coverage_report
 /data/COVID_WGS/lix33/Test/2ndpipeline/Data/secondary_buf/PRE_QC
 
 
 /data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Data/BAM/Batch/low_input*/Flag
 /data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Data/BAM/Batch/low_input*/Log
 /data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Data/BAM/Batch/low_input*/Report
 
 
 /data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Keytable/low_input
<<<<<<<<<
'''

import os
import sys
import re
import time
import subprocess


DATE = "10_27_2021"
KEYTABLEName = "std_input_01_361_keytable"
KEYWORD = "std_input"

class ClsDirSet:
    def __init__(self):
        self.strS3Root = ""
        
        self.strBiowulfBuidDir = ""
        self.strBiowulfBAMOrgDir = ""
        self.strBiowulfBAMRecalibrateDir = ""
        
        self.strBiowulfSecondBufLogsDir = ""
        self.strBiowulfSecondBufReportDir = ""
        self.strBiowulfSecondBufPreQCDir = ""
                
        self.arryBiowulfUpstreamFlag = []
        self.arryBiowulfUpstreamLog = []
        self.arryBiowulfUpstreamReport = []        
        
        self.strBiowulfUpstreamKeyTableDir = ""
        
        #-> Use it to remove data from biowulf
        self.arryBatchDir = []
    
    def Init(self):
        #->S3 root
        self.strS3Root = "lix33/COVID19/UpstreamAnalysis/PostPrimaryRun/Data/" + DATE + "/" + KEYTABLEName
        
        #-> Build Dir
        # --> this may contain softlink -> need to zip first
        self.strBiowulfBuidDir = "/data/COVID_WGS/lix33/Test/2ndpipeline/Build"
        self.strBiowulfBAMOrgDir = "/data/COVID_WGS/lix33/Test/2ndpipeline/Data/BAM_reformatted/BAM_original/WGS"
        self.strBiowulfBAMRecalibrateDir = " /data/COVID_WGS/lix33/Test/2ndpipeline/Data/BAM_reformatted/BAM_recalibrated/WGS"
        
        #-> SecondBuf Dir
        self.strBiowulfSecondBufLogsDir = "/data/COVID_WGS/lix33/Test/2ndpipeline/Data/secondary_buf/cluster_job_logs"
        self.strBiowulfSecondBufReportDir = "/data/COVID_WGS/lix33/Test/2ndpipeline/Data/secondary_buf/coverage_report"
        self.strBiowulfSecondBufPreQCDir = "/data/COVID_WGS/lix33/Test/2ndpipeline/Data/secondary_buf/PRE_QC"
        
        
        strBiowulfUpstreamBatchDir = "/data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Data/BAM/Batch"
        CMD = "find " + strBiowulfUpstreamBatchDir + " -maxdepth 1 -type d -iname '" + KEYWORD + "*'"
        #print(CMD)
        vDir = subprocess.getoutput(CMD).split('\n')
        self.arryBatchDir = vDir 
        #print(vDir)
        for strDir in vDir:
            CMD = "find " + strDir + " -mindepth 1 -maxdepth 1 -type d"
            vSubDir = subprocess.getoutput(CMD).split('\n')
            #print(vSubDir)
            for strSubDir in vSubDir:
                if "/Flag" in strSubDir:
                    self.arryBiowulfUpstreamFlag.append(strSubDir)
                    continue
                if "/Report" in strSubDir:                
                    self.arryBiowulfUpstreamReport.append(strSubDir)
                    continue
                if "/Log" in strSubDir:                
                    self.arryBiowulfUpstreamLog.append(strSubDir)
                    continue        
                
        print(self.arryBiowulfUpstreamFlag)
        print(self.arryBiowulfUpstreamReport)
        print(self.arryBiowulfUpstreamLog)        

        self.strBiowulfUpstreamKeyTableDir = "/data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Keytable/" + KEYWORD     
        
        print("DirSet Init has been finished!")
        print()        
        
    def Upload(self):
        BackupDirFromBiowulf2S3(self.strS3Root, self.strBiowulfBuidDir, "/data/COVID_WGS/lix33/Test/")
        BackupDirFromBiowulf2S3(self.strS3Root, self.strBiowulfBAMOrgDir, "/data/COVID_WGS/lix33/Test/")
        BackupDirFromBiowulf2S3(self.strS3Root, self.strBiowulfBAMRecalibrateDir, "/data/COVID_WGS/lix33/Test/")
        
        BackupDirFromBiowulf2S3(self.strS3Root, self.strBiowulfSecondBufLogsDir, "/data/COVID_WGS/lix33/Test/")
        BackupDirFromBiowulf2S3(self.strS3Root, self.strBiowulfSecondBufReportDir, "/data/COVID_WGS/lix33/Test/")
        BackupDirFromBiowulf2S3(self.strS3Root, self.strBiowulfSecondBufPreQCDir, "/data/COVID_WGS/lix33/Test/")
        
        for strFlagDir in self.arryBiowulfUpstreamFlag:
            BackupDirFromBiowulf2S3(self.strS3Root, strFlagDir, "/data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Data/BAM/")
        
        for strLogDir in self.arryBiowulfUpstreamLog:
            BackupDirFromBiowulf2S3(self.strS3Root, strLogDir, "/data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Data/BAM/")
        
        for strReportDir in self.arryBiowulfUpstreamReport:
            BackupDirFromBiowulf2S3(self.strS3Root, strReportDir, "/data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Data/BAM/")
            
        BackupDirFromBiowulf2S3(self.strS3Root, self.strBiowulfUpstreamKeyTableDir, "/data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/")

    def Remove(self):
        print("Remove data from Biowulf -->", "\n")
        
        print("******* Remove BiowulfBuidDir, BiowulfBAMOrgDir and BiowulfBAMRecalibrateDir **********")
        RemoveDirFromBiowulf(self.strBiowulfBuidDir)
        RemoveDirFromBiowulf(self.strBiowulfBAMOrgDir)
        RemoveDirFromBiowulf(self.strBiowulfBAMRecalibrateDir)
        
        print("******* Remove BiowulfSecondBufLogsDir, BiowulfSecondBufReportDir and BiowulfSecondBufPreQCDir **********")
        RemoveDirFromBiowulf(self.strBiowulfSecondBufLogsDir)
        RemoveDirFromBiowulf(self.strBiowulfSecondBufReportDir)
        RemoveDirFromBiowulf(self.strBiowulfSecondBufPreQCDir)
        
        print("******* Remove UpstreamAnalysisBatchDirs **********")
        for strDir in self.arryBatchDir:
            RemoveDirFromBiowulf(strDir)
        

def RemoveDirFromBiowulf(strDir):
    if not os.path.exists(strDir):
        print("Error: Path do not exist!")
        return 1
    
    print("Removing Dir:", strDir)
    CMD = "rm -r " + strDir
    os.system(CMD)
    print("Removing Done!", "\n")
            
def BackupDirFromBiowulf2S3(strS3Root, strOrgDir, strCutoffPrefix, bContainSoftLink=False):
    for subdir, dirs, files in os.walk(strOrgDir):
        if len(files) == 0:
            continue
        elif "/bam_location/WGS" in subdir:
            strZipFile = subdir.split("bam_location")[0] + "bam_location_softlink.zip"
            if os.path.exists(strZipFile):
                continue
            else:
                CMD = "zip --symlinks -r " + strZipFile + " " + subdir.split("/WGS/")[0]
                #print(CMD)
                os.system(CMD)
                #Copy it to S3
                strPrefix = strS3Root + "/" + os.path.dirname(strZipFile).split(strCutoffPrefix)[1] + "/"
                #CMD = "obj_put -v DCEG_COVID_WGS -p " + "\"" + strPrefix + "\" " + strZipFile + " --dry-run -V" # Dry run
                CMD = "obj_put -v DCEG_COVID_WGS -p " + "\"" + strPrefix + "\" " + strZipFile + " -V" # Wet run
                os.system(CMD)                 
                continue
        else:
            for file in files:
                strPrefix = strS3Root + "/" + subdir.split(strCutoffPrefix)[1] + "/"
                #CMD = "obj_put -v DCEG_COVID_WGS -p " + "\"" + strPrefix + "\" " + os.path.join(subdir, file) + " --dry-run -V" # Dry run
                CMD = "obj_put -v DCEG_COVID_WGS -p " + "\"" + strPrefix + "\" " + os.path.join(subdir, file) + " -V" # Wet run
                os.system(CMD)
                #print(os.path.join(subdir, file))
                #Copy it to S3
                    
    # strPrefix = strObjPrefix + os.path.dirname(strFile).replace(strRootDir, '') + "/"        
    # #CMD = "obj_put -v DCEG_COVID_WGS -p " + "\"" + strPrefix + "\" " + strFile + " --dry-run -V" # Dry run
    # CMD = "obj_put -v DCEG_COVID_WGS -p " + "\"" + strPrefix + "\" " + strFile + " -V" # wet run
    

def BackUpToS3(objDirSet):
    objDirSet.Upload()
    print("BackUpToS3 is all set!")


def RemoveFromBiowulf(objDirSet):
    objDirSet.Remove()
    print("RemoveFromBiowulf is all set!")

def main():
    objDirSet = ClsDirSet()
    objDirSet.Init()
    
    BackUpToS3(objDirSet)
    #RemoveFromBiowulf(objDirSet)


if __name__ == "__main__":    
    print("====== Calculation: Backup COVID 2nd Pipeline Results ======", flush=True)
    print()
    start = time.time()
    rValue = main()
    end = time.time()    
    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    print()
    print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds), flush=True)
    print("====== End: Backup COVID 2nd Pipeline Results ======", flush=True)
    
    
    