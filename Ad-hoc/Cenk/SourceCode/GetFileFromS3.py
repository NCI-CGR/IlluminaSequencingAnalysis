'''
1: Check if the samples defined in excel can be found in S3
2: Download the existing samples from S3
3: Output a list of missing samples
'''
import subprocess
import sys
import os

DESTDir = "/data/COVID_ADHOC/Sequencing/COVID_WGS/Data/USU_First_Batch_low_std_96_361_18"
LOGDir = DESTDir + "/Log" 

class ClsSample:
    def __init__(self):
        self.strCGRID = ""
        self.strBAM = ""
        self.bDup = False
    
    def Init(self, strCGRID):
        self.strCGRID = strCGRID
    
    def RetrieveDataFromS3(self):
        # 1: Check if BAM is existed
        if self.strBAM == "":
            print("Error: BAM file is not existed in S3!", " ---> Sample ID:", self.strBAM)
            return
        
        if not os.path.exists(DESTDir):
            print("Error: Dest Dir is not existed!")
            return
        
        if self.bDup:
            print("Warning: duplicate sample, skip it!", " ---> Sample ID:", self.strBAM)
            return
        
        strFileName = os.path.basename(self.strBAM)
        strBashScript = ("bash " + "RetrieveFileFromS3.sh " + " '" + self.strBAM + "'" + 
                                                              " '" + strFileName + "'" +
                                                              " '" + DESTDir + "'")
        
        # 2: Submit jobs
        strNumCore = "1"
        strNodeNum = "1"
        strJobName = "RTS3-" + self.strCGRID
                     
        strRunningTime = "3-00:00:00"
        strStdOut = LOGDir + "/" + strFileName.split('.')[0] + ".std.out"
        strStdErr = LOGDir + "/" + strFileName.split('.')[0] + ".std.err"        
                        
        CMD = ("sbatch " + "--ntasks=" + strNumCore + " " +  
                           "--nodes=" + strNodeNum + " " +
                           "--job-name=" + strJobName + " " + 
                           "--time=" + strRunningTime + " " + 
                           "--output=" + strStdOut + " " +
                           "--error=" + strStdErr + " " + 
                           "--wrap=\"" + strBashScript + "\"")
        
        print(CMD)
        os.system(CMD)                                
        print("Retrieving BAM ->", self.strCGRID, "\n")   

def InitSamples(vSample):
    CMD = "awk -F ',' '{print $1}' ../Files/COVIDwgs_Jan2022_NIAIDphenotypesLM_Comma_Delimited.csv | tail -n +2"
    strIDList = subprocess.getoutput(CMD)
    vIDList = strIDList.split('\n')
    vSample.clear()
    for strID in vIDList:
         objSample = ClsSample()
         objSample.Init(strID)
         vSample.append(objSample)
    
    #Check if there is any duplicate samples
    dictIDList = {}
    for sample in vSample:
        if sample.strCGRID in dictIDList:
            dictIDList[sample.strCGRID] += 1
            sample.bDup = True/data/COVID_ADHOC/Sequencing/COVID_WGS/Data/USU_First_Batch_low_std_96_361_18/Log
        else:
            dictIDList[sample.strCGRID] = 1
    
    iDupNum = 0
    for key in dictIDList:
        if dictIDList[key] > 1:
            print(key, ",", dictIDList[key])
            iDupNum += 1
    print(" ********", "\n", "iDupNum:", iDupNum, "\n", "********", "\n")
        

def MatchSampleWithS3Archive(vSample):
    #Get S3 recalibrated BAM list from S3
    CMD = "obj_ls -v DCEG_COVID_WGS -m '*recalibrated*.bam' | tail -n +2 | awk '{print $NF}'"
    strBAMList = subprocess.getoutput(CMD)
    vBAMList = strBAMList.split('\n')
    dictBAMList = {}
    for strBAM in vBAMList:
        strCGRID = os.path.basename(strBAM).split('_')[1]
        if not strCGRID in dictBAMList:
              dictBAMList[strCGRID] = strBAM
    for sample in vSample:
        if sample.strCGRID in dictBAMList:
            sample.strBAM = dictBAMList[sample.strCGRID]
    
    # Check the number of Sample that do not contain BAM Path
    iEmptyBAM = 0
    for sample in vSample:
        if sample.strBAM == "":
            print(sample.strCGRID)
            iEmptyBAM += 1
    print(" ********", "\n", "iEmptyBAM:", iEmptyBAM, "\n", "********", "\n")
    
def RetrieveDataFromS3(vSample):
    for sample in vSample:
        sample.RetrieveDataFromS3()
        #return

def main():
    #strExcel= sys.argv[1]
    vSample = []
    InitSamples(vSample)
    MatchSampleWithS3Archive(vSample)
    RetrieveDataFromS3(vSample)

if __name__ == "__main__":
    main()