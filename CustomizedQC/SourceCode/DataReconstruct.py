'''
The purpose of this code is convert the data structure into 
the data structure which is compatible with our pipeline
'''
import sys
import os
import subprocess
from datetime import date
import gzip

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

class ClsSample:
    def __init__(self):
        self.strReads1 = ""
        self.strReads2 = ""
        self.strSampleName = ""
        self.strFlowcellName = ""
        self.strLane = ""
        self.strBarcode = ""        
    
    def AddReads(self, strReads, strIndex):
        if strIndex == "R1":
            self.strReads1 = strReads
        elif strIndex == "R2":
            self.strReads2 = strReads
    
    def CheckIdentical(self, strFlowcellName, strSampleName, strLane, strBarcode):
        if (self.strFlowcellName == strFlowcellName and 
            self.strSampleName == strSampleName and 
            self.strLane == strLane and 
            self.strBarcode == strBarcode):
            return True
        else:
            return False
    
    def Init(self, strFlowcellName, strSampleName, strBarcode, strLane, strReads, strReadsIndex):
        self.strFlowcellName = strFlowcellName
        self.strSampleName = strSampleName
        self.strBarcode = strBarcode
        self.strLane = strLane
        #print("Add Reads:", strFlowcellName)
        self.AddReads(strReads, strReadsIndex)
        
    
    def CheckValidation(self):
        if self.strReads1 == "" or self.strReads2 == "":
            print("Error: Invalid Sample -> ", self.strSampleName, ",", self.strBarcode, ",", self.strLane)
            print("Reads 1:", self.strReads1)
            print("Reads 2:", self.strReads2)
            return False
        else:
            return True
    
    def Reconstruct(self, strCASAVA):
        strOutputDir = strCASAVA + "/" + self.strLane + "/Project_COVID19"
        CMD = "mkdir -p " + strOutputDir
        os.system(CMD)  
        
        # 1: Create Sample Folder        
        strFolderName = "Sample_" + self.strSampleName + "_" + self.strBarcode + "_" + self.strLane 
        strFolderPath = strOutputDir + "/" + strFolderName 
        CMD = "mkdir -p " + strFolderPath
        os.system(CMD)
        
        #2: Create soft link for reads file
        # 1) reads 1
        strReadsName = self.strSampleName + "_" + self.strBarcode + "_" + self.strLane + "_" + "R1_001_HQ_paired.fastq.gz"
        strReadsPath = strFolderPath + "/" + strReadsName 
        CMD = "ln -s " + self.strReads1 + " " + strReadsPath
        os.system(CMD)   
             
        # 2) reads 2
        strReadsName = self.strSampleName + "_" + self.strBarcode + "_" + self.strLane + "_" + "R2_001_HQ_paired.fastq.gz"
        strReadsPath = strFolderPath + "/" + strReadsName 
        CMD = "ln -s " + self.strReads2 + " " + strReadsPath
        os.system(CMD)

class ClsFlowcell:
    def __init__(self):
        self.strName = ""
        self.vSample = []
        self.strTimeStamp = ""
    
    def AddSample(self, objSample):
        self.vSample.append(objSample)
    
    def Init(self, objSample):
        self.strName = objSample.strFlowcellName
        self.vSample.clear()
        self.vSample.append(objSample)
        #--> Get Tiem Stamp
        todays_date = date.today()
        strYear = str(todays_date.year)
        strMonth = str(todays_date.month)
        if int(strMonth) < 10: 
            strMonth = "0" + strMonth
        self.strTimeStamp = strYear + strMonth   
        #<--
        
    def Reconstruct(self, strRootOutputDir):
        # 1: Check if this flowcell is existed (folder contains the flowcell's name)
        for subDir in os.listdir(strRootOutputDir):
            strSubDirPath = strRootOutputDir + "/" + subDir
            if os.path.isdir(strSubDirPath) and self.strName in subDir:
                return
        
        # 2: Create Flowcell Folder
        strFlowcellDir = strRootOutputDir + "/" + self.strTimeStamp + "_XXXXXX_XXXX_" + self.strName
        CMD = "mkdir -p " + strFlowcellDir
        os.system(CMD)        
        
        # 3: create folder to put samples
        strCASAVA = strFlowcellDir + "/CASAVA"
        CMD = "mkdir -p " + strFlowcellDir
        os.system(CMD)
        
        # 4: Put sample into CASAVA folder
        for sample in self.vSample:
            sample.Reconstruct(strCASAVA)
        
def AddSample(vSample, strReads):
    # 1: Parse current reads
    #Check the first line of Reads
    #print(strReads)
    #(1) Get flowcell name and barcode from reads content!
    strReadsInfo = ""
    if is_gz_file(strReads):
        with gzip.open(strReads, 'rb') as f1:
            for x in f1:                
                strLine = x.decode("utf-8")
                strReadsInfo = strLine.split('\n')[0] #  this is very  tricky!! Important!
                break
    else:        
        CMD = "less " + strReads + " | head -n 1"
        #print(CMD)
        strReadsInfo = subprocess.getoutput(CMD)
    
    #print(strReadsInfo)
    vItem = strReadsInfo.split(':')
    strFlowcellName = vItem[2]
    strBarcode = vItem[-1].replace('+', '-')
    
    #(2) Get SampleName, ReadsIndex and Lane Index from reads' name
    vItem = os.path.basename(strReads).split('.')[0].split('_')
    strSampleName = vItem[0]
    strReadsIndex = vItem[-2] 
    strLane = vItem[-3]
    #print("strLane:", strLane)
    
    # 2: Check current saved sample
    bFind = False
    for sample in vSample:
        if sample.CheckIdentical(strFlowcellName, strSampleName, strLane, strBarcode):
            #update the existing sample
            #print("Add to existing sample")
            sample.AddReads(strReads, strReadsIndex)
            bFind = True
            break
    
    if not bFind:
        #print("Create a new sample")
        objSample = ClsSample()
        objSample.Init(strFlowcellName, strSampleName, strBarcode, strLane, strReads, strReadsIndex)
        vSample.append(objSample)
    
    #print("---")
        
def CheckValidation(vSample):
    bValid = True
    for sample in vSample:
        if not sample.CheckValidation():
            bValid = False
            break
    return bValid

def GetFlowcells(vFlowcell, vSample):
    vFlowcell.clear()
    for sample in vSample:
        bFind = False
        for flowcell in vFlowcell:
            if sample.strFlowcellName == flowcell.strName:
                flowcell.AddSample(sample)
                bFind = True
                break
        if not bFind:
            objFlowcell = ClsFlowcell()
            objFlowcell.Init(sample)            
            vFlowcell.append(objFlowcell)    

def Reconstruct(strRootOutputDir, vFlowcell):
    # Create output Dir
    if not os.path.exists(strRootOutputDir):
        CMD = "mkdir -p " + strRootOutputDir
        os.system(CMD)
        
    for flowcell in vFlowcell:
        # Generate related Sample Folder and Read files (softlink)
        flowcell.Reconstruct(strRootOutputDir)    
    print()
    print("Total number of flowcells:", len(vFlowcell))
    print("----", '\n', "Flowcell Details (Flowcell Name + Sample Number)")
    for flowcell in vFlowcell:
        print(flowcell.strName, "->", len(flowcell.vSample))

def main():
    # Get arguments and check if the folder path is available
    strRootDir = sys.argv[1]
    strOutputDir = sys.argv[2]
    if not os.path.exists(strRootDir):
        print("Error: Path does not exist!")
        return
    
    #Scan fastq files and get sample (reads pair)
    vSample = []
    for subDir in os.listdir(strRootDir):
        strSubDirPath = strRootDir + "/" + subDir
        if os.path.isdir(strSubDirPath): 
            for file in os.listdir(strSubDirPath):
                if file.endswith(".fastq.gz"):
                    #print(file)
                    #Add Sample into Sample list
                    AddSample(vSample, strSubDirPath + "/" + file)
    
    #Check if the samples we parsed are valid
    print("CheckValidation -->")
    bValid = CheckValidation(vSample)
    if not bValid:
        print("Some samples in current dataset are invalid, please check the error info!")
        return
    
    #Cluster these sample by using their flowcell ID
    print('\n', "GetFlowcells -->")
    vFlowcell = []
    GetFlowcells(vFlowcell, vSample)
    
    # Reshape data if everything is OK -> Go
    print('\n', "Reconstruct -->")
    Reconstruct(strOutputDir, vFlowcell)

if __name__ == "__main__":
    main()
    