'''
1: Extract sample from excel 
2: Get related BAM from local dir
3: Run samtools view
'''
import subprocess
import sys
import os

BAMDir = "/data/COVID_ADHOC/Sequencing/COVID_WGS/Data/USU_First_Batch_low_std_96_361_18"
LOGDir = BAMDir + "/Log"

DESTDir = "/data/COVID_ADHOC/Sequencing/COVID_WGS/ad-hoc/Cenk/TargetBAM"
DESTLogDir = DESTDir + "/Log"


class ClsSample:
    def __init__(self):
        self.strCGRID = ""
        self.strBAM = ""
        self.bDup = False
    
    def Init(self, strCGRID):
        self.strCGRID = strCGRID
    
    def GetTargetReadsByBed(self):
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
        
        # 2: Submit jobs
        strNumCore = "12"
        strNodeNum = "1"
        strJobName = "Target-" + self.strCGRID
        
        strFileName = os.path.basename(self.strBAM).split('.')[0] + ".target.bam"
        
        # self.strBAM = "/data/COVID_WGS/lix33/Test/2ndpipeline/Data/BAM_reformatted/BAM_recalibrated/WGS/WGS_SC689824_I3-97832.bam"
        # strFileName = "WGS_SC689824_I3-97832.target.bam"
        
        strBashScript = ("bash SliceTargetReadsByBed.sh" + " '" + self.strBAM + "'" + 
                                                           " '" + strFileName + "'" +
                                                           " '" + DESTDir + "'" + 
                                                           " '" + strNumCore + "'")
                     
        strRunningTime = "3-00:00:00"
        strStdOut = DESTLogDir + "/" + strFileName.split('.')[0] + ".std.out"
        strStdErr = DESTLogDir + "/" + strFileName.split('.')[0] + ".std.err"        
                        
        CMD = ("sbatch " + "--ntasks=" + strNumCore + " " +  
                           "--nodes=" + strNodeNum + " " +
                           "--job-name=" + strJobName + " " + 
                           "--time=" + strRunningTime + " " + 
                           "--output=" + strStdOut + " " +
                           "--error=" + strStdErr + " " + 
                           "--wrap=\"" + strBashScript + "\"")
        
        print(CMD)
        os.system(CMD)                                
        print("Get Target BAM ->", self.strCGRID, "\n") 
        

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
            sample.bDup = True
        else:
            dictIDList[sample.strCGRID] = 1
    
    iDupNum = 0
    for key in dictIDList:
        if dictIDList[key] > 1:
            print(key, ",", dictIDList[key])
            iDupNum += 1
    print(" ********", "\n", "iDupNum:", iDupNum, "\n", "********", "\n")

def FindBAMInLocal(vSample):
    CMD = "find " + BAMDir + " -type f -iname '*.bam'"
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

def GetTargetReadsByBed(vSample):
    for sample in vSample:
        sample.GetTargetReadsByBed()
        #return

def main():
    #strExcel= sys.argv[1]
    vSample = []
    InitSamples(vSample)
    FindBAMInLocal(vSample)
    GetTargetReadsByBed(vSample)

if __name__ == "__main__":
    main()