'''
Input: keytable
e.g.
/data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Keytable/low_input/sub_keytable/low_input_01_50_keytable.csv

Output:
the running results after 2nd pipeline be finished!
'''

import sys
import os

DIRCustomizedQC = "/home/lix33/lxwg/Git/IlluminaSequencingAnalysis/CustomizedQC/SourceCode"
DIR2ndPipeline = "/home/lix33/lxwg/Git/IlluminaSequencingAnalysis/COVID_2nd_Pipeline/SourceCode"

DIRBAMRoot = "/data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Data/BAM/Batch"


class ClsPhase:
    def __init__(self):
        self.strCmd = ""

def MergeBAM(strKeytable):
    # Pre check if everything is all set! -->
    strScript = DIRCustomizedQC + "/Tools/MergeSubject/MergeSubject.py"
    CMD = "python " + strScript + " " + strKeytable
    os.system(CMD)
    #<-- 

def main():
    strKeytable = sys.argv[1]
    if not os.path.exists(strKeytable):
        print("Error: Keytable does not exist! -->", strKeytable)
        return
    # Step 1: Merge BAM file
    MergeBAM(strKeytable)
    
    # Step 2:
    
    # Step 3:
    
    # Step 4:
    
    # Step 5:     
    
    # Step 6:
    
    # Step 7:
    
    

if __name__ == "__main__":
    main()