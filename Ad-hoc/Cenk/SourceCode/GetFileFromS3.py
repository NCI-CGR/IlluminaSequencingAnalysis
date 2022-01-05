'''
1: Check if the samples defined in excel can be found in S3
2: Download the existing samples from S3
3: Output a list of missing samples
'''
import subprocess
import sys
import os

class ClsSample:
    def __init__(self):
        self.strCGRID = ""
        self.strBAM = ""        
    
    def Init(self, strCGRID):
        self.strCGRID = strCGRID 


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

def main():
    #strExcel= sys.argv[1]
    vSample = []
    InitSamples(vSample)
    MatchSampleWithS3Archive(vSample)

if __name__ == "__main__":
    main()