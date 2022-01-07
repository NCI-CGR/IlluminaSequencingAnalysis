'''
1: Update CSV based on the finished BAM 
2: Transfer CSV to excel
3: Notify the related people
'''

SOURCEDir = "/data/COVID_ADHOC/Sequencing/COVID_WGS/ad-hoc/Cenk/TargetBAM"
DESTDir = "/data/COVID_ADHOC/Sequencing/COVID_WGS/ad-hoc/Cenk"

import pandas as pd
import subprocess
import sys
import os

def Test():
    df = pd.read_csv('../Files/COVIDwgs_Jan2022_NIAIDphenotypesLM_Comma_Delimited.csv')
    df["BAM_File_Name"] = ""
    for i, row in df.iterrows():
        df.at[i, 'BAM_File_Name'] = str(i) + "_lxwg"
        print(df.at[i, 'LIMS_Sample_ID'])    
    #print(df)
    df.to_excel('../Files/tmp.xlsx')
    
def UpdateExcel():
    if not os.path.exists(SOURCEDir) or not os.path.exists(DESTDir):
        print("Error: No such Directory!", " --> ", SOURCEDir, "\n", " --> ", DESTDir)
        return 1
    
    # 1: Get BAM from current
    CMD = "find " + SOURCEDir + " -type f -iname " + "*.bam"    
    strBAMList = subprocess.getoutput(CMD)
    vBAMList = strBAMList.split('\n')
    
    # Build dictionary based on current BAM List
    print(vBAMList)
    dictBAM = {}
    for strBAM in vBAMList:
        strFileName = os.path.basename(strBAM)
        strCGRID = strFileName.split('_')[1]                 
        if strCGRID not in dictIDList:
            dictBAM[strCGRID] = strFileName
    
    # Update Excel data framework
    df = pd.read_csv('../Files/COVIDwgs_Jan2022_NIAIDphenotypesLM_Comma_Delimited.csv')
    df["BAM_File_Name"] = ""
    for i, row in df.iterrows():
        strCGRID = df.at[i, 'LIMS_Sample_ID']
        if strCGRID in dictBAM:
            df.at[i, 'BAM_File_Name'] = dictBAM[strCGRID]
            
    # Export Excel
    df.to_excel(DESTDir + "/COVIDwgs_Jan2022_NIAIDphenotypesLM_Append_BAM.csv")

def main():
    UpdateExcel()    

if __name__ == "__main__":
    main()