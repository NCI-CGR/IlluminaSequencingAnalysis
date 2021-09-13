'''
This is for merge the same sample from different flowcells

How to do it:
1: Input: 
    CSV file (contains the topoff sample list)
2: Output:
    The merged subject
3: Logic
(1) Collect group info from CSV file 
(2) Get the path for each sample 
(3) Run secondary pipeline to merge them
(4) Output the merged result (check log file)
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
        self.strBAM = ""
    
    def InitByKTLine(self, strline):
        vItem = strline.split(',')
        #print(vItem)
        self.strFlowcellID = vItem[0].split('_')[-1][1:]
        self.strUSUID = vItem[1]
        self.strCGRID = vItem[2].split('-')[-1]
        if len(vItem) == 4 and vItem[3] == "topoff":
            self.bTopoff = True
    
    def GetBAMFile(self):
        # file BAM file -> need to grab data from S3 storage space
        print(self.strCGRID)
        if self.strUSUID == "SC571825":
            self.strBAM = "/home/lixin/lxwg/project/MergeSample/2ndPipeline/BAM/SC571825_CTTCACCA-AAGAGCCA_L001.bam"
        
        if self.strUSUID == "SC571826":
            self.strBAM = "/home/lixin/lxwg/project/MergeSample/2ndPipeline/BAM/SC571826_TACCAGGA-GTACTCTC_L001.bam"     

    def PrepareManifest(self, vContent):
        if not os.path.exists(self.strBAM):
            return
        # Prepare line for current sample
        CMD = "samtools view -H " + self.strBAM + " | grep -iP '@RG\t'"
        strRGList = subprocess.getoutput(CMD)
        if strRGList == "":
            return 
                
        strRGLine = strRGList.split('\n')[0]
        vRG = strRGLine.split('\t')
        print(vRG)
        # --> Prepare info for current sample
        strInstrument = "E0180-03"
        strSeqDate = ""
        strLane = "1"
        strIndex = vRG[-2].split('.')[-1]
        strFullFlowcellID = vRG[-2].split('.')[0].split(':')[1]
        print(strFullFlowcellID)
        print(strIndex)
        strGroup = "ACS"
        strLMSIndividual = "I-0000510733"
        strEXPECTEDGENDER = "M"
        strIDENTIFILERGENDER = "M"
        strUCSCAVGCOV = "99.99"
        strASSAYID = "EZ_Exome+UTR_PE"
        strANALYSISID = "NA12878_GDNA_HS4K_KHPL_1"
        
        strLine = (strInstrument + "," + 
                   strSeqDate + "," +
                   strFullFlowcellID + "," +
                   strLane  + "," +
                   strIndex  + "," +
                   self.strUSUID  + "," +
                   strGroup  + "," +
                   strLMSIndividual  + "," +
                   strEXPECTEDGENDER  + "," +
                   strIDENTIFILERGENDER  + "," +
                   strUCSCAVGCOV  + "," +
                   strASSAYID  + "," +
                   strANALYSISID  + "," +
                   self.strCGRID + "\n")
        vContent.append(strLine)
        print(strLine)
        # <--
        
    
class ClsSubject:
    def __init__(self):
        self.strCGRID = ""
        self.vSample = []
    
    def GetBAMFile(self):
        for sample in self.vSample:
            sample.GetBAMFile()
    
    def PrepareManifest(self, vContent):
        for sample in self.vSample:
            sample.PrepareManifest(vContent)
            
class ClsBuild:
    def __init__(self):
        self.vSubject = []
        self.strManifestFile = ""
    
    # it may contain both 1,2,and multiple files for each subject (1 BAM for a single subject is allowed)
    def GetGroupInfo(self, strKeyTable):
        if not os.path.exists(strKeyTable):
            return        
        # 1: Get group info        
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
                # Check if the subject has been existed
                bFind = False
                for subject in self.vSubject:
                    if subject.strCGRID == objSample.strCGRID:
                        subject.vSample.append(objSample)
                        bFind = True
                        break
                if not bFind:
                    objSubject = ClsSubject()
                    objSubject.strCGRID = objSample.strCGRID
                    objSubject.vSample.append(objSample)
                    self.vSubject.append(objSubject)
                iIndex += 1
        file.close()
    
    def GetBAMFile(self):
        for subject in self.vSubject:
            subject.GetBAMFile()
    
    def BuidManifestFile(self, strDir):
        strFile = strDir + "/Manifest_Mimic.csv"
        if os.path.exists(strFile):
            CMD = "rm " + strFile
            os.system(CMD)
        
        vContent = []
        # Add head        
        strHead = "INSTRUMENT,SEQDATE,FLOWCELL,LANE,INDEX,CGF ID,GROUP,LIMS INDVIDUALID,EXPECTED GENDER,IDENTIFILER GENDER,UCSCAVGCOV,ASSAYID,ANALYSIS ID,SR SUBJECT ID\n"
        vContent.append(strHead)
        
        # Prepare content for each sample
        for subject in self.vSubject:
            subject.PrepareManifest(vContent)
        
        # Write to file            
        f = open(strFile, "w")
        f.writelines(vContent)        
        f.close()

def main():
    strKeyTable = sys.argv[1]
    
    objBuild = ClsBuild()
    
    # 1: Prepare group info  
    objBuild.GetGroupInfo(strKeyTable)
    
    # 2: Get BAM file for each sample
    objBuild.GetBAMFile()
    print(len(objBuild.vSubject))
    
    # 3: Build Manifest file -> GO!
    objBuild.BuidManifestFile(os.path.dirname(strKeyTable))    
    
    print("All Set!")

if __name__ == "__main__":
    main()

