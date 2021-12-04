import os
import sys
import subprocess

REPORTType = { "1": "Coverage_Report",
               "2": "Pre_Calling_QC_Report",
               "3": "BAM_Contamination_Check_Report"}

REPORTSuffix = {"1": ".txt",
                "2": ".txt",
                "3": ".csv"}

REPORTKeywords = { "1": "coverage_report",
                   "2": "pre_calling_qc_report",
                   "3": "ContaminationReport"}

def main():
    # Check
    strReportType = sys.argv[1]
    strReportDir = sys.argv[2]
    strName = sys.argv[3]
    
    if strReportType not in REPORTType.keys():
        print("Error: Wrong Key value!")
        print("Please check the list below: ")
        print("\"1\": \"Coverage_Report\"", "\n", 
              "\"2\": \"Pre_Calling_QC_Report\"", "\n",
              "\"3\": \"BAM_Contamination_Check_Report\"", "\n")
        return 1
    
    if not os.path.exists(strReportDir):
        print("Error: Report Dir does not exist!")
        return 1
    
    suffix = REPORTSuffix[strReportType]
    strFileName = strName.split('.')[0] + "_" + REPORTType[strReportType] + suffix
    strSumReportFile = strReportDir + "/" + strFileName
    if os.path.exists(strSumReportFile):
        print("Warning: everything was all set! No further action needed!")
        print("File location:", strSumReportFile)
        return 1
    
    # Find all files
    CMD = "find " + strReportDir + " -maxdepth 2 -type f -iname " + "'*" + REPORTKeywords[strReportType] + "*'" 
    strFileList = subprocess.getoutput(CMD)
    if strFileList == "":
        print("Warning: No file detected! No further action needed!")        
        return 1
    
    vFileList = strFileList.split('\n')
     
    bHeadFilled = False
    for strFile in vFileList:
        if not bHeadFilled:
            CMD = "cat " + strFile + " > " + strSumReportFile
            os.system(CMD)
            bHeadFilled = True 
        else:
            CMD = "awk 'NR > 1' " + strFile + " >> " + strSumReportFile
            os.system(CMD)
    
    # remove x permission 
    CMD = "chmod -x " + strSumReportFile
    os.system(CMD)
    
    print("Good: Everythign is all set!")
    print("Report type  :", REPORTType[strReportType])
    print("File Name    :", strFileName)
    print("File Location:", strSumReportFile)
    return 0

if __name__ == "__main__":
    main()
    