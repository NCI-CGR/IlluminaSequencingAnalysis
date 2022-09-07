'''
Generate the sample list for the delivered sample
'''
import subprocess
import sys
import os

BAMFolder = "/data/COVID_ADHOC/Sequencing/COVID_WGS/ad-hoc/Cenk/TargetBAM_TCR_IGK"
DELIVERSampleListFile = "/data/COVID_ADHOC/Sequencing/COVID_WGS/ad-hoc/Cenk/UpdatedExcel/TargetBAM_TCR_IGK/Delivered_sample_list_09_07_2022.txt"

def GenerateStatisticSampleList():
    strFilePath = BAMFolder
    
    CMD = 'find ' + strFilePath + " -iname '*.bam'"
    vLine = subprocess.getoutput(CMD).split('\n')
    strFileName = DELIVERSampleListFile
    osf = open(strFileName, "w")
    strRow = "SequenceType\tReadsType\tUSU_ID\tCGR_ID\n"
    osf.write(strRow) 
    for strLine in vLine:
        vKeywords = os.path.basename(strLine).split('.')[0].split('_')
        strRow = ""
        if len(vKeywords) == 3:
            strRow = vKeywords[0] + "\t" + "Std" + "\t" + vKeywords[2] + "\t" + vKeywords[1] + "\n"
        elif len(vKeywords) == 4:
            strRow = vKeywords[0] + "\t" + "Std+Topoff" + "\t" + vKeywords[2] + "+" + vKeywords[3]  + "\t" + vKeywords[1] + "\n"
        osf.write(strRow)
    osf.close()
def main():
    GenerateStatisticSampleList()    

if __name__ == "__main__":
    main()
