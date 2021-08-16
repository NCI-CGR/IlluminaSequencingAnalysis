'''
Merge different QC report into a big summary file
'''
import sys
import os
import subprocess

class ClsFlowcell:
    def __init__(self):
        self.strDir = ""
        self.strQCReport = ""
        self.strName = ""
        self.strDate = ""
    
    def Init(self, strDir):
        # 1: Set dir         
        self.strDir = strDir
        
        # 2: Find file 
        CMD = "find " + strDir + " -maxdepth 1 -type f -name \"*.csv\""        
        #print(CMD)
        strFastqList = subprocess.getoutput(CMD)
        if CMD == "":
            returm 
        vCSVList = strFastqList.split('\n') 
        self.strQCReport = vCSVList[0]
        
        # 3: Parse Name and Date
        strFolderName = os.path.basename(strDir)
        self.strName = strFolderName.split('_')[-1]
        self.strDate = strFolderName.split('_')[0] 
               
def CheckValidation(vFlowcell):
    if len(vFlowcell) == 0:
        print("Error: No flowcell be found!")
        return 1
    for flowcell in vFlowcell:
        if flowcell.strQCReport == "":
            print("Error flowcell:", flowcell.strDir)
            return 2
    return 0

def MergeQCReport(strDir, vFlowcell):
    # Get suffix of sum report
    strSuffix = os.path.basename(vFlowcell[0].strQCReport)
    if '_' in strSuffix:
        strSuffix = strSuffix.split('_')[1] 
         
    strMergedFile = strDir + "/" + "Sum_" + strSuffix
    
    #1: Remove old file first 
    if os.path.exists(strMergedFile):
        CMD = "rm -f " + strMergedFile
        os.system(CMD)
    
    #2: Generate merged report
    # (1) Get fields
    CMD = "head -n 1 " + vFlowcell[0].strQCReport
    strFields = subprocess.getoutput(CMD)
    strFields = "FlowcellName,Date," + strFields
     
    CMD = "echo " + "\"" + strFields + "\" >> " + strMergedFile
    os.system(CMD)  
    
    for flowcell in vFlowcell:
        # Using readlines()
        with open(flowcell.strQCReport) as fp:
            strLine = fp.readline()
            iCount = 0
            while strLine:    
                #print(strLine)            
                if iCount != 0:                
                    CMD = "echo " + "\"" + flowcell.strName + "," + flowcell.strDate + "," + strLine.strip() + "\" >> " + strMergedFile
                    os.system(CMD)
                iCount += 1
                strLine = fp.readline()            
    
    print("Merging result:", strMergedFile)
    print("All Set!")           

def main():
    strProcessedDir = sys.argv[1]
    
    if not os.path.exists(strProcessedDir):
        print("Error: Folder is not existed:", strProcessedDir)

        
    vFlowcell = []
    for subDir in os.listdir(strProcessedDir):
        strSubDirPath = strProcessedDir + "/" + subDir
        if os.path.isdir(strSubDirPath):
            objFlowcell = ClsFlowcell()
            objFlowcell.Init(strSubDirPath)  
            vFlowcell.append(objFlowcell)
    
    if CheckValidation(vFlowcell) != 0:
        print("Error: QC report is missing!")
    
    MergeQCReport(strProcessedDir, vFlowcell)    
    
if __name__ == "__main__":
    main()