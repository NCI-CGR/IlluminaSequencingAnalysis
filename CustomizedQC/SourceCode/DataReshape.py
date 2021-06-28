'''
The purpose of this code is convert the data structure into 
the data structure which is compatible with our pipeline
'''
import sys
import os
import subprocess

class ClsSample:
    def __init__(self):
        self.strReads1 = ""
        self.strReads2 = ""
        self.strName = ""
        self.strLane = ""
        self.strBarcode = ""
    
    def AddReads(self, strReads, strIndex):
        if strIndex == "R1":
            self.strReads1 = strReads
        elif strIndex == "R2":
            self.strReads2 = strReads
    
    def CheckIdentical(self, strSampleName, strLane, strBarcode):
        if (self.strName == strSampleName and 
            self.strLane == strLane and 
            self.strBarcode == strBarcode):
            return True
        else:
            return False
    
    def Init(self, strSampleName, strBarcode, strLane, strReads, strReadsIndex):
        self.strName = strSampleName
        self.strBarcode = strBarcode
        self.strLane = strLane
        self.AddReads(strReads, strReadsIndex)
    
    def CheckValidation(self):
        if self.strReads1 == "" or self.strReads2 == "":
            print("Error: Invalid Sample -> ", self.strName, ",", self.strBarcode, ",", self.strLane)
            print("Reads 1:", self.strReads1)
            print("Reads 2:", self.strReads2)
            return False
        else:
            return True
    
    def Reshape(self, strOutputDir):        
        # 1: Create Sample Folder        
        strFolderName = "Sample_" + self.strName + "_" + self.strLane + "_" + self.strBarcode
        strFolderPath = strOutputDir + "/" + strFolderName 
        CMD = "mkdir -p " + strFolderPath
        os.system(CMD)
        
        #2: Create soft link for reads file
        # 1) reads 1
        strReadsName = self.strName + "_" + self.strBarcode + "_" + self.strLane + "_" + "R1_001.fastq.gz"
        strReadsPath = strFolderPath + "/" + strReadsName 
        CMD = "ln -s " + self.strReads1 + " " + strReadsPath
        os.system(CMD)   
             
        # 2) reads 2
        strReadsName = self.strName + "_" + self.strBarcode + "_" + self.strLane + "_" + "R2_001.fastq.gz"
        strReadsPath = strFolderPath + "/" + strReadsName 
        CMD = "ln -s " + self.strReads2 + " " + strReadsPath
        os.system(CMD)
        
def AddSample(vSample, strReads):
    # 1: Parse current reads
    #Check the first line of Reads
    #print(strReads)
    CMD = "less " + strReads + " | head -n 1"
    strReadsInfo = subprocess.getoutput(CMD)
    #print(strReadsInfo)
    vItem = strReadsInfo.split(':')
    strSampleName = vItem[2]
    strBarcode = vItem[-1].replace('+', '-')
    strLane = os.path.basename(strReads).split('.')[0].split('_')[-3]
    strReadsIndex = os.path.basename(strReads).split('.')[0].split('_')[-2]
 
    # 2: Check current saved sample
    bFind = False
    for sample in vSample:
        if sample.CheckIdentical(strSampleName, strLane, strBarcode):
            #update the existing sample
            sample.AddReads(strReads, strReadsIndex)
            bFind = True
            break
    
    if not bFind:
        objSample = ClsSample()
        objSample.Init(strSampleName, strBarcode, strLane, strReads, strReadsIndex)
        vSample.append(objSample)
        
def CheckValidation(vSample):
    bValid = True
    for sample in vSample:
        if not sample.CheckValidation():
            bValid = False
            break
    return bValid

def Reshape(strRootOutputDir, vSample):
    # Create output Dir
    if not os.path.exists(strRootOutputDir):
        CMD = "mkdir -p " + strRootOutputDir
        os.system(CMD)
    
    # Generate related Sample Folder and Read files (softlink) 
    for sample in vSample:
        sample.Reshape(strRootOutputDir)
    
    print()
    print("Total number of Samples:", len(vSample))
    
def main():
    # Get arguments and check if the folder path is available
    strRootDir = sys.argv[1]
    strOutputDir = sys.argv[2]
    if not os.path.exists(strRootDir):
        print("Error: Path does not exist!")
        return
    
    #Scan fastq files and get sample (reads pair)
    vSample = []
    for file in os.listdir(strRootDir):
        if file.endswith(".fastq.gz"):
            #Add Sample into Sample list
            AddSample(vSample, strRootDir + "/" + file)
    
    #Check if the samples we parsed are valid
    bValid = CheckValidation(vSample)
    if not bValid:
        print("Some samples in current dataset are invalid, please check the error info!")
        return
    
    # Reshape data if everything is OK -> Go
    strRootFolderName = os.path.basename(strRootDir)
    strRootOutputDir = strOutputDir + "/" + strRootFolderName  
    Reshape(strRootOutputDir, vSample)

if __name__ == "__main__":
    main()
    