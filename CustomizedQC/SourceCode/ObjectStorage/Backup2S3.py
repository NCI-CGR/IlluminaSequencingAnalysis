'''
1: Back up the finished flowcells in biowulf target root dir to Object Storage system
2: Scan all flowcells in target root dir
3: How to check if a flowcell is all set
  (1) Done flag is existing
  (2) QC report has been generated
'''

import subprocess
import sys
import os

def BackupFlowcell2S3(strFlowcellDir, strRootDir, strObjPrefix):
    print("Start to back up flowcell:", strFlowcellDir, " ---->")
    # Step1: Check if such flowcell is already existed in s3 -> Do it later
    # Step2: Backup
    CMD = "find " + strFlowcellDir + " -type f"    
    strFileList = subprocess.getoutput(CMD)
    vFile = strFileList.split('\n')
    for strFile in vFile:
        if strFile == "":
            continue
        
        strPrefix = strObjPrefix + os.path.dirname(strFile).replace(strRootDir, '') + "/"        
        #CMD = "obj_put -v DCEG_COVID_WGS -p " + "\"" + strPrefix + "\" " + strFile + " --dry-run -V" # Dry run
        CMD = "obj_put -v DCEG_COVID_WGS -p " + "\"" + strPrefix + "\" " + strFile + " -V" # wet run
        print(CMD)
        iReturn = os.system(CMD)
        if iReturn != 0:
            print("\n", ">>>>>> Error Occured!")
            print(CMD, "\n", "<<<<<<<", "\n") 
            exit()
    print("\n", "**********")
    print("Current Flowcell is all set:", strFlowcellDir)
    print("**********", "\n")
            
def RemoveFlowcellFromBiowulf(strFlowcellDir):
    print("Removing finished flowcell -> ")
    if os.path.exists(strFlowcellDir):
        CMD = "rm -r " + strFlowcellDir
        iReturn = os.system(CMD)
        if iReturn != 0:
            print("\n", ">>>>>> Error Occured!")            
            print(CMD, "\n", "<<<<<<<", "\n") 
            exit()
    print("\n", "**********")
    print("Flowcell has been removed from Biowulf:", strFlowcellDir)
    print("<<<<<<<<<<<<<<", "\n")

def main():
    strDataDir = "/data/COVID_WGS/primary_analysis/COVID19/06_30_2021/ProcessedData"
    strRootDir = "/data/COVID_WGS/primary_analysis/COVID19"
    strObjPrefix = "lix33/COVID19/USU/Data"
        
    #Get all of files 
    CMD = "find " + strDataDir + " -mindepth 1 -maxdepth 1 -type d"
    strFlowcellList = subprocess.getoutput(CMD)
    vFlowcell = strFlowcellList.split('\n')
    for strFlowcellDir in vFlowcell:
        if strFlowcellDir == "":
            continue
        if not os.path.exists(strFlowcellDir):
            continue
        
        # Check done flag and QC report to see if this flowcell is all set!
        # 1: Check flag
        bAllDone = False
        strDoneFlag = strFlowcellDir + "/Flag/BWA/v38/flag.all.done"
        if os.path.exists(strDoneFlag):
            bAllDone = True
        
        # 2: Check QC report
        bQCReportGenerated = False
        CMD = "find " + strFlowcellDir + " -maxdepth 1 -type f -iname \"*QCReport-BWA-v38.csv\""
        strQCReportList = subprocess.getoutput(CMD)
        if strQCReportList != "":
            vQCReport = strQCReportList.split('\n')
            if len(vQCReport) == 1 and os.path.exists(vQCReport[0]):
                bQCReportGenerated = True
        
        if bAllDone and bQCReportGenerated:
            BackupFlowcell2S3(strFlowcellDir, strRootDir, strObjPrefix)
            RemoveFlowcellFromBiowulf(strFlowcellDir)
                
    print("All Set!")

if __name__ == "__main__":    
    main()