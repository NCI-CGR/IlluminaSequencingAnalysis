'''
1: Cluster Sample types of QC 
   a) LowIput
   b) StdInput   
2: Add CGR ID based on the keytable info between USU ID and CGR ID
3: Add Top off info
'''
import os
import sys
import subprocess

class ClsSample:
    def __init__(self):
        self.strUSUID = ""
        self.strCGRID = ""
        self.strFlowcellID = ""
        self.strBarcode = ""
        self.bTopoff = False
        self.strQCInfo = ""
        self.strLaneIndex = ""
    
    def CheckSame(self, objSample):
        bSame = False
        if self.strUSUID == objSample.strUSUID and self.strFlowcellID == objSample.strFlowcellID:
            bSame = True
        return bSame

    def InitBySampleDir(self, samnpleDir):
        # Get Flowcell ID
        strFlowcellDirName = samnpleDir.split('/')[-5]
        strFlowcellID = strFlowcellDirName.split('_')[-1]
        self.strFlowcellID = strFlowcellID
        
        # Get remaining info    
        strDirName = os.path.basename(samnpleDir)
        vItem = strDirName.split('_')        
        self.strUSUID = vItem[1]
        self.strBarcode = vItem[2]        
        self.strLaneIndex = vItem[3]
        #self.Print()
    
    def InitByKTLine(self, strline):
        vItem = strline.split(',')            
        self.strFlowcellID = vItem[0].split('_')[-1][1:]
        self.strUSUID = vItem[1]
        self.strCGRID = vItem[2].split('-')[-1]
        if len(vItem) == 4 and vItem[3] == "topoff":
            self.bTopoff = True
        #self.Print()
    
    def InitByCSVLine(self, strLine, strFlowcellID):
        self.strFlowcellID = strFlowcellID
        vItem = strLine.split(',')
        self.strUSUID = vItem[0]
        self.strBarcode = vItem[1]
        self.strLaneIndex = "L00" + vItem[2]
        self.strQCInfo = strLine            
        #self.Print()
            
    def Print(self):
        print("********")
        print("strUSUID     :", self.strUSUID)
        print("strCGRID     :", self.strCGRID)
        print("strFlowcellID:", self.strFlowcellID)
        print("strBarcode   :", self.strBarcode)
        print("bTopoff      :", self.bTopoff)
        print("strQCInfo    :", self.strQCInfo)
        print("strLaneIndex :", self.strLaneIndex)
        print("********", '\n')

class ClsQCReport:
    def __init__(self):
        self.vSample = []
        self.strFields = ""
        self.strRoot = ""
    
    def GenerateQC4KT(self, strType):
        strQCReport = self.strRoot + "/" + strType + ".Sum_QC_Report.csv"
        print(strQCReport)
        if os.path.exists(strQCReport):
            CMD = "rm -f " + strQCReport
            os.system(CMD)
        
        file = open(strQCReport, 'a')
        
        # Update Fields
        # Add two additioal column (CGR ID and Topoff status) 
        vFields = self.strFields.split('Sample,')        
        stFields = "Flowcell ID,Sample,CGR ID,Topoff Status," + vFields[1]
        #print(stFields) 
        file.write(stFields + "\n")        
         
        # Add QC Info for each sample
        for sample in self.vSample:
            vItem = sample.strQCInfo.split(',')
            strInfo = (sample.strFlowcellID + "," +
                       vItem[0] + "," + 
                       sample.strCGRID  + "," +
                       str(sample.bTopoff)  + "," + 
                       ','.join(vItem[1:]))            
            #print(strInfo)
            file.write(strInfo + "\n")
        file.close()
        print("All Set!")
        
def GetCurBatchSample(vCurBatchSample, strBatchDir):
    vCurBatchSample.clear()
    
    CMD = "find " + strBatchDir + " -maxdepth 5 -type d -iname 'Sample_*'"
    strSampleList = subprocess.getoutput(CMD)
    vSampleDir = strSampleList.split('\n')    
    
    for samnpleDir in vSampleDir:        
        objSample = ClsSample()
        objSample.InitBySampleDir(samnpleDir)
        vCurBatchSample.append(objSample)
              
    print("vCurBatchrSample:", len(vCurBatchSample))
    print("GetCurBatchSample is All Set", '\n')
    
def GetKTSample(vKTSample, strKeyTable):
    vKTSample.clear()
    
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
            vKTSample.append(objSample)
            iIndex += 1            
        
    file.close()      
    print("vKTSample:", len(vKTSample))
    print("GetKTSample is All Set", '\n')
    
def GetQCReportSample(objQCReport, vQCReportSample, strQCReportDir):
    # 1: Get All csv file
    CMD = "find " + strQCReportDir + " -maxdepth 2 -type f -iname '*.csv'"
    strCSVList = subprocess.getoutput(CMD)
    vCSVFile = strCSVList.split('\n')
    
    # 2: Get Sample and field
    vQCReportSample.clear()
    bSetFiled = False
    strFields = ""
    for csvFile in vCSVFile:
        # Get Flowcell Name 
        strFlowcellName = csvFile.split('/')[-2]
        strFlowcellID = strFlowcellName.split('_')[-1]  
        #print(csvFile)
        file = open(csvFile, 'r')    
        iIndex = 0
        while True:        
            strline = file.readline()
            if not strline:
                break;
            if iIndex == 0:
                if not bSetFiled:
                    strFields = strline.strip()
                    bSetFiled = True
                iIndex += 1
                continue
            else:
                objSample = ClsSample()
                objSample.InitByCSVLine(strline.strip(), strFlowcellID)
                vQCReportSample.append(objSample)
                iIndex += 1                
        
    #print(strField)
    objQCReport.strFields = strFields
    print("vQCReportSample:", len(vQCReportSample))          
    print("GetQCReportSample is All Set", '\n')

def UpdateCurBatchSample(vCurBatchSample, vKTSample, vQCReportSample):
    for curSample in  vCurBatchSample:        
        # Update CGR ID and Topoff from vKTSample
        #print(len(vKTSample))
        for ktSample in vKTSample:
            if ktSample.CheckSame(curSample):                
                curSample.strCGRID = ktSample.strCGRID
                curSample.bTopoff = ktSample.bTopoff
                break

        # Update QC Info from  vQCReportSample
        #print(len(vQCReportSample))
        for qcrSample in vQCReportSample:
            if qcrSample.CheckSame(curSample):
                curSample.strQCInfo = qcrSample.strQCInfo                
                break
    for curSample in vCurBatchSample:
        curSample.Print()

    # Update Fields (2 additional)    
    print("vCurBatchSample has been updated!", '\n')

def UpdateKTSample(vKTSample, vQCReportSample):
    for ktSample in  vKTSample:        
        # Update QC Info from  vQCReportSample
        #print(len(vQCReportSample))
        for qcrSample in vQCReportSample:
            if qcrSample.CheckSame(ktSample):
                ktSample.strBarcode = qcrSample.strBarcode
                ktSample.strLaneIndex = qcrSample.strLaneIndex
                ktSample.strQCInfo = qcrSample.strQCInfo                
                break

    # for ktSample in vKTSample:
    #     ktSample.Print()
        
    # Update Fields (2 additional)    
    print("UpdateKTSample has been updated!", '\n')

def GenerateQC4CurBatch(vCuBatchrSample, strBatchDir):
    print("GenerateQC4CurBatch is All Set", '\n')

'''
Input:
1: LowInput/StdInput Folder
2: Keytable
3: Original QC Report Folder

example:
Low Input　-->
Arg1: /home/lixin/lxwg/ad-hoc/COVID/08_14_2021/LowInput/Processed
Arg2: /home/lixin/lxwg/ad-hoc/COVID/08_14_2021/Keytable/NCI-USU_COVNET_docs/Sample_Manifests-QC_Information/Batch_2_Low_Input_WGS_pilot/Keytable_and_file_inventory_from_USU/pI8.tagmentation.N96.keytable.jul.23.2021.csv
Arg3: /home/lixin/lxwg/ad-hoc/COVID/08_14_2021/Summary/QCReport

Std Input  -->
Arg1: /home/lixin/lxwg/ad-hoc/COVID/08_14_2021/StdInput/Processed
Arg2: /home/lixin/lxwg/ad-hoc/COVID/08_14_2021/Keytable/NCI-USU_COVNET_docs/Sample_Manifests-QC_Information/Batch_1_Standard_WGS/Keytable_and_file_inventory_from_USU/pI3.COVNET.N361.keytable.june.3.2021+pI3.COVNET.N18.15topoff.keytable.june.30.2021.csv
Arg3: /home/lixin/lxwg/ad-hoc/COVID/08_14_2021/Summary/QCReport
'''
def main():
    strBatchDir = sys.argv[1]
    strKeyTable = sys.argv[2]
    strQCReportDir = sys.argv[3]
    
    print("strBatchDir   :", strBatchDir)
    print("strKeyTable   :", strKeyTable)
    print("strQCReportDir:", strQCReportDir, "\n")
    
    # Step 1: Get all samples in current batch -> Done
    vCurBatchSample = []
    GetCurBatchSample(vCurBatchSample, strBatchDir)
    
    # Step 2: Get all sample info from Keytable -> Done
    vKTSample = []
    GetKTSample(vKTSample, strKeyTable)
    
    # Step 3: Get all sample from QC report -> Go
    objQCReport = ClsQCReport()    
    vQCReportSample = []    
    GetQCReportSample(objQCReport, vQCReportSample, strQCReportDir)
        
    # # (Step 4: Update vCurBatchSample)
    # UpdateCurBatchSample(vCurBatchSample, vKTSample, vQCReportSample)
    # Step 4: Update vKTSample
    UpdateKTSample(vKTSample, vQCReportSample)
    
    # Step ４: Generate QC Report for current Batch
    objQCReport.vSample = vKTSample
    objQCReport.strRoot = strBatchDir
    strType = os.path.basename(strKeyTable).split('.csv')[0] + ".std_input"
    print(strType)
    objQCReport.GenerateQC4KT(strType)

if __name__ == "__main__":
    main()