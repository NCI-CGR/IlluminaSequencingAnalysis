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

def func_match(s1, s2):
    return sum(s1[i] != s2[i] for i in range(len(s1))) <= 1

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
    
    def Reconstruct(self, strFolderPath):
        # Create soft link for reads file
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
        

class ClsSubject:
    def __init__(self):
        self.strSampleName = ""
        self.strFlowcellName = ""
        self.strLane = ""
        self.strBarcode = ""
        self.vSample = []
    
    def Init(self, objSample):
        self.vSample.clear()
        self.strSampleName = objSample.strSampleName
        self.strFlowcellName = objSample.strFlowcellName
        self.strLane = objSample.strLane
        self.strBarcode = objSample.strBarcode
        self.vSample.append(objSample)
            
    def UpdateLaneIndex(self, strCurSampleLaneIndex):
        if int(self.strLane[-1]) > int(strCurSampleLaneIndex[-1]):
            self.strLane = strCurSampleLaneIndex     
    
    def CheckIdentical(self, strSampleName, strBarcode):
        if (self.strSampleName == strSampleName and             
            func_match(self.strBarcode, strBarcode)):
            return True
        else:
            return False 
    
    def Reconstruct(self, strCASAVA):
        strOutputDir = strCASAVA + "/" + self.strLane + "/Project_COVID19"
        CMD = "mkdir -p " + strOutputDir
        os.system(CMD)  
        
        # 1: Create Sample Folder        
        strFolderName = "Sample_" + self.strSampleName + "_" + self.strBarcode + "_" + self.strLane 
        strFolderPath = strOutputDir + "/" + strFolderName 
        CMD = "mkdir -p " + strFolderPath
        os.system(CMD)
        
        for sample in self.vSample:
            sample.Reconstruct(strFolderPath)
        

class ClsFlowcell:
    def __init__(self):
        self.strName = ""
        self.vSample = []
        self.strTimeStamp = ""
        self.vSubject = []
        # parse ID and name
        self.dictIDName = {} 
    
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
    
    def PrepareSubject(self):
        self.vSubject.clear()
        for sample in self.vSample:
            bFind = False
            for subject in self.vSubject:
                if subject.CheckIdentical(sample.strSampleName, sample.strBarcode):
                    subject.vSample.append(sample)
                    subject.UpdateLaneIndex(sample.strLane)
                    bFind = True
            if not bFind:
                #-> It is a new subject
                objSubject = ClsSubject()
                objSubject.Init(sample)
                self.vSubject.append(objSubject)
            
        
    def Reconstruct(self, strRootOutputDir):
        # 1: Check if this flowcell is existed (folder contains the flowcell's name)
        for subDir in os.listdir(strRootOutputDir):
            strSubDirPath = strRootOutputDir + "/" + subDir
            if os.path.isdir(strSubDirPath) and self.strName in subDir:
                print("WARNING: Flowcell has already existed (skip):", strSubDirPath)
                return
        
        # 2: Create Flowcell Folder
        strFlowcellDir = strRootOutputDir + "/" + self.strTimeStamp + "_XXXXXX_XXXX_" + self.strName
        CMD = "mkdir -p " + strFlowcellDir
        os.system(CMD)        
        
        # 3: create folder to put samples
        strCASAVA = strFlowcellDir + "/CASAVA"
        CMD = "mkdir -p " + strFlowcellDir
        os.system(CMD)
        
        # 4: Prepare subject -> Go!        
        self.PrepareSubject()
                
        # 5: Put sample into CASAVA folder
        for subject in self.vSubject:
            subject.Reconstruct(strCASAVA)

        
def AddSample(vSample, strReads):
    # 1: Parse current reads
    #Check the first line of Reads
    #print(strReads)
    #(1) Get flowcell name and barcode from reads content!    
    strReadsInfo = ""
    # -> obtain the first "iMaxReadsNum" reads and choose the most common barcode. (this is because that the barcode may have 1 bp different)
    #    Notice: we use the first 5 reads before, however, it is not enough based on the real data, then we are using the first 101 reads to reduce the bias. 
    vLine = []
    iStartReadsNum = 100 # we start to count the reads if the reads index > "iStartReadsNum" 
                         # motivation: we plan to skip the first some parts of reads, since the barcode in 
                         #             these reads may contain many uncertain things "e.g. contain 'N'")
    iMaxReadsNum = 101 # total number of reads we need to collect
    #print(strReads)
    if is_gz_file(strReads):
        i = 1
        targetLine = 1
        iCount = 0
        with gzip.open(strReads, 'rb') as f1:            
            for x in f1:                
                strLine = x.decode("utf-8")
                if i == targetLine:
                    iCount += 1
                    if iCount > iStartReadsNum:
                        vLine.append(strLine.split('\n')[0]) #  this is very tricky!! Important!
                    targetLine = targetLine + 4
                i += 1 
                if len(vLine) == iMaxReadsNum:
                    break
    else:        
        CMD = "less " + strReads + " | awk 'NR % 4 == 1' | head -n " + str(iStartReadsNum + iMaxReadsNum) + " | tail -n " + str(iMaxReadsNum)
        #os.system(CMD)
        #print(subprocess.getoutput(CMD))
        vLine = subprocess.getoutput(CMD).split('\n')
    
    #print(vLine)
    # Get the most common barcode
    dicBarcode = {}
    for strLine in vLine:
        bFind = False
        strBarcode = strLine.split(':')[-1]
        for key in dicBarcode:
            if key == strBarcode:
                dicBarcode[key] += 1
                bFind = True
                break
        if not bFind:
            dicBarcode[strBarcode] = 1
            
    #print(dicBarcode)
    listSortedBarcode = sorted(dicBarcode.items(), key=lambda x: x[1], reverse=True)
    #print(dicSortedBarcode)
    #print(dicBarcode)
    strMostCommonBarcode = listSortedBarcode[0][0]
    #print(strMostCommonBarcode)
    for strLine in vLine:
        if strLine.split(':')[-1] == strMostCommonBarcode:
            strReadsInfo = strLine
            break    
    
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
    
    #Print the info of flowcell
    print("Number of Flowcell:", len(vFlowcell)) 
    for flowcell in vFlowcell:
        print(flowcell.strName, "-->", len(flowcell.vSample))        

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
        print(flowcell.strName, "->", len(flowcell.vSubject))

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
    
    #  
    
    # Reshape data if everything is OK -> Go
    print('\n', "Reconstruct -->")
    Reconstruct(strOutputDir, vFlowcell)

if __name__ == "__main__":
    main()
    