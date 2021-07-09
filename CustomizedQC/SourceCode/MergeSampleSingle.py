'''
1: Check if current sample need to be merged
2: Merge the related sample
3: manage the flags for each sample
'''

import os
import sys
import subprocess

FLAGMultiSample = "flag.multi.sample"

def main():
    strSampleDir = sys.argv[1]
    strWorkingFlag = sys.argv[2]
    strDoneFlag = sys.argv[3]
    
    # Check how many softlinks in sample dir
    CMD = "find " + strSampleDir + " -maxdepth 1 -type l -name \"*.fastq.gz\""
    strFastqList = subprocess.getoutput(CMD)
    vFastqList = strFastqList.split('\n')
    
    if len(vFastqList) == 2:
        # This is the normal case, and we do not need to do anything
        # 1: delete working flag
        if os.path.exists(strWorkingFlag):
            CMD = "rm " + strWorkingFlag
            os.system(CMD)
        # 2: Create done flag
        CMD = "touch " + strDoneFlag
        os.system(CMD)
        print("This is the normal case -> Update working and done flag directly!")
        print("--> All Set")
    elif len(vFastqList) > 2 and len(vFastqList) % 2 == 0:
        # Collect all soft link and group them by using reads index
        strR1KW = "R1_001_HQ_paired.fastq.gz"
        strR2KW = "R2_001_HQ_paired.fastq.gz"
        
        vReads1 = []
        vReads2 = []
        
        for file in vFastqList:
            if strR1KW in file:
                vReads1.append(file)
                continue
            if strR2KW in file:
                vReads2.append(file)
                continue
        
        # Combine all fastq files for reads1:
        vReads1.sort()
        strDir = os.path.dirname(vReads1[0])
        strTmpFile = os.path.basename(vReads1[0])
        strNewName = strTmpFile.split('_')[0] + "-sum_" + '_'.join(strTmpFile.split('_')[1:])          
        strSumFile = strDir + "/" + strNewName  
        CMD = "cat"
        for reads in vReads1:     
            CMD += " " + reads
        CMD += " > " + strSumFile
        os.system(CMD)
        print("SumFile(Reads 1):", strSumFile)
                
        # Combine all fastq files for reads2:          
        vReads2.sort()
        strDir = os.path.dirname(vReads2[0])
        strTmpFile = os.path.basename(vReads2[0])
        strNewName = strTmpFile.split('_')[0] + "-sum_" + '_'.join(strTmpFile.split('_')[1:])          
        strSumFile = strDir + "/" + strNewName  
        CMD = "cat"
        for reads in vReads2:     
            CMD += " " + reads
        CMD += " > " + strSumFile
        os.system(CMD)
        print("SumFile(Reads 2):", strSumFile)
        
        # update flag
        # 1: create the special flag to let people know we combined fastq files
        CMD = "touch " + os.path.dirname(strDoneFlag) + "/" + FLAGMultiSample
        os.system(CMD)
        
        # 2: delete working flag
        if os.path.exists(strWorkingFlag):
            CMD = "rm " + strWorkingFlag
            os.system(CMD)
        # 3: Create done flag
        CMD = "touch " + strDoneFlag
        os.system(CMD)                        
        print("--> All Set")
        
    else:
        print("Error: no sample files or odd number of samples been found!")


if __name__ == "__main__":
    main()