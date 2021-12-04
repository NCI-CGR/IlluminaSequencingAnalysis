'''
1: Input key_table: folder
/data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/keytable/low_input/sub_keytable/low_input_01_10_keytable.csv 

2: Output Data: folder
/data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Data/BAM/Batch

Purpose:
1: retrieve data from S3 to Biowulf
'''

import os
import sys


DIRBAMRoot = "/data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Data/BAM/Batch"


class ClsSample:
    def __init__(self):
        self.strUSUID = ""
        self.strCGRID = ""
        self.strFlowcellID = ""
        self.strBarcode = ""
        self.bTopoff = False
    
    def InitByKTLine(self, strLine):
        vItem = strLine.split(',')
        #print(vItem)
        self.strFlowcellID = vItem[0].split('_')[-1][1:]
        self.strUSUID = vItem[1]
        self.strCGRID = vItem[2].split('-')[-1]
        if len(vItem) == 4 and vItem[3] == "topoff":
            self.bTopoff = True
    
    def RetrieveBAM(self, strCurBAMDir):
        strLogDir = strCurBAMDir + "/Log/RetrieveBAM"
        strFlagDir = strCurBAMDir + "/Flag/RetrieveBAM"
        
        if not os.path.exists(strLogDir):
            CMD = "mkdir -p " + strLogDir
            os.system(CMD)
        
        if not os.path.exists(strFlagDir):
            CMD = "mkdir -p " + strFlagDir
            os.system(CMD)
        
        outputPrefix = "retrieveBAM_" + self.strUSUID + "_" + self.strCGRID + "_" + self.strFlowcellID
        strFlagWorking = strFlagDir + "/" + outputPrefix + ".working"
        strFlagDone = strFlagDir + "/" + outputPrefix + ".done"
        
        if os.path.exists(strFlagWorking):
            print("It is still running! -->", outputPrefix)
            return
        elif os.path.exists(strFlagDone):
            print("It has been finished! -->", outputPrefix)
            return
        else:
            CMD = "touch " + strFlagWorking
            os.system(CMD) 
        
        # Start to prepare job for current sample
        strCurDir = os.path.dirname(os.path.abspath(__file__))
        strScript = strCurDir + "/RetrieveUSUBAM.sh"
        #strScript = strCurDir + "/RetrieveUSU2ndBAMOrg.sh"

        # CMD = ("obj_ls -v DCEG_COVID_WGS -h -m " + 
        #        "*" + self.strFlowcellID + "*" + self.strUSUID +  "*dedup_nophix.??? " + 
        #        "|awk 'NR>1 {print $8}' | " + 
        #        "while read line; do obj_get -v DCEG_COVID_WGS ${line} -D /low_input_01_20_keytable -p -V --strip 6 --dry-run; done")
        # print(CMD)
        strBashScript = ("bash " + strScript + " " +
                                               self.strFlowcellID + " " +
                                               self.strUSUID + " " +
                                               self.strCGRID + " " +
                                               strCurBAMDir + " " + 
                                               strFlagWorking + " " + 
                                               strFlagDone)        
        
        #Submit jobs ->
        strNumCore = "1"
        strNodeNum = "1"
        strJobName = "RTV-" + self.strUSUID + "-" + self.strCGRID + "-" + self.strFlowcellID
        
             
        strRunningTime = "1-00:00:00"
        strStdOut = strLogDir + "/" + outputPrefix + ".std.out"
        strStdErr = strLogDir + "/" + outputPrefix + ".std.err"        
                        
        CMD = ("sbatch " + "--ntasks=" + strNumCore + " " +  
                           "--nodes=" + strNodeNum + " " +
                           "--job-name=" + strJobName + " " + 
                           "--time=" + strRunningTime + " " + 
                           "--output=" + strStdOut + " " +
                           "--error=" + strStdErr + " " + 
                           "--wrap=\"" + strBashScript + "\"")
        
        print(CMD)
        os.system(CMD)                                
        print("Retrieving BAM ->", self.strUSUID, "\n")        
    
def GetSampleFromKeytable(strKeyTable, vSample):
    vSample.clear()
    file = open(strKeyTable, 'r')    
    iIndex = 0
    while True:        
        strline = file.readline()
        if not strline:
            break;
        if iIndex == 0:
            iIndex += 1
            continue
        else:
            objSample = ClsSample()
            objSample.InitByKTLine(strline.strip())
            vSample.append(objSample)
            iIndex += 1
    file.close()


def main():
    strKeyTable = sys.argv[1]
    if not os.path.exists(strKeyTable):
        print("Error: Keytable does NOT exist!", strKeyTable)
        return
    
    strCurBAMDir = DIRBAMRoot + "/" + os.path.basename(strKeyTable).split('.')[0]
    if not os.path.exists(strCurBAMDir):
        CMD = "mkdir -p " + strCurBAMDir
        print(CMD)
        os.system(CMD)
    
    vSample = []
    GetSampleFromKeytable(strKeyTable, vSample)
    
    for sample in vSample:
        sample.RetrieveBAM(strCurBAMDir)


if __name__ == "__main__":
    main()





