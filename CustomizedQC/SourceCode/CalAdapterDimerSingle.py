#Notice this is calculation is for a single sample
#This script is called by primary pipeline

import gzip
import os
import sys
import re
import time

def IsDimerReads(reads):
    #Skip the reads with length > 50 (the adapter dimer reads should be with length 35(33 + 1 + 1))     
    lenReads = len(reads)
    if lenReads > 50:
        return False        
    bDimerReads = bool(re.match('^[N]+$', reads.upper()))    
    return bDimerReads

#Main Function ------
def main():        
    #1: Get arguments     
    sampleFolderPath = sys.argv[1]
    resultFilePath = sys.argv[2]
    
    #2: Check if the sample is valid
    if not os.path.isdir(sampleFolderPath):
        print("Error (file path doesn't exist): ", sampleFolderPath, flush=True)
        return
    
    #3: Check if the calculation result is valid          
    if os.path.isfile(resultFilePath):
        print("Dimer Calculation has been finished previously. No more action required.", flush=True)
        return   
    
    #3: Collect all valid fastq files under sampleRootPath
    arryTrimFile = []    
    arryReadsNum = []
    arryDimerNum = []
    arryIgnoreFile = ["singleton", "UMI"] #["singleton", "HQ_paired", "UMI"]    
    
    for r, d, f in os.walk(sampleFolderPath):
        for file in f:            
            if '.fastq.gz' in file:
                if any(x in file for x in arryIgnoreFile):
                    continue
                else:
                    arryTrimFile.append(os.path.join(r, file))
                    
    #4:Check each fastq    
    for fqFile in arryTrimFile:
        i = 1
        targetLine = 2        
        iCurReadsNum = 0
        iCurDimerNum = 0
        #print(fqFile)
        with gzip.open(fqFile, 'rb') as f1:            
            for x in f1:                
                strReads = x.decode("utf-8") #remember --> decode first                
                if i == targetLine:        
                    iCurReadsNum += 1
                    targetLine = targetLine + 4                                        
                    #if bool(re.match('^[N]+$', strReads.upper())): #check if only contain N
                    if IsDimerReads(strReads):                        
                        iCurDimerNum += 1                                 
                i = i + 1
        arryReadsNum.append(iCurReadsNum)
        arryDimerNum.append(iCurDimerNum)
    
    #5: Check result
    readsSum = sum(arryReadsNum)
    dimerSum = sum(arryDimerNum)
    dimerRatio = 'Nil'
    if readsSum != 0: 
        dimerRatio = '{:.5%}'.format(dimerSum / readsSum) 
    
    #Print result and save it to log
    print("Result: \t Num_of_reads \t Num_of_Dimer \t Dimer_Ratio", file=open(resultFilePath, "a"), flush=True)
    print(readsSum, "\t", dimerSum, "\t", dimerRatio, file=open(resultFilePath, "a"), flush=True)
    print("Calculation is finished successfully! Please check the result:", resultFilePath, flush=True)
                    
# Require 2 arguments, including: sampleFolderPath, tmpResultFilePath
if __name__ == "__main__":
    print("====== Calculation: Adapter Dimer Ratio ======", flush=True)
    start = time.time()
    rValue = main()
    end = time.time()    
    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds), flush=True)
    print("====== End: Adapter Dimer Ratio ======", flush=True)
    
    