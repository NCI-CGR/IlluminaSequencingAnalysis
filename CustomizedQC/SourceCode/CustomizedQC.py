'''
Do Alignment and then generate QC

Processed Data folder Structure

../ProcessedData/
├── BAM
├── Flag
├── Log
├── reports
└── Sample
'''

#1: Build contains a group of samples

import os
import sys
import subprocess
TMPFolder = "/scratch/lix33/TmpSequencing/"
REFSeq = "/scratch/lix33/DCEG/CGF/Bioinformatics/Production/data/hg38/Homo_sapiens_assembly38.fasta"
REFSamIndex = "/scratch/lix33/DCEG/CGF/Bioinformatics/Production/data/hg38/Homo_sapiens_assembly38.fasta.fai"

EMAILSender = "xin.li4@nih.gov"
EMAILReceiver = "xin.li4@nih.gov"
#EMAILReceiver = "xin.li4@nih.gov,mingyi.wang@nih.gov"

FLAGAlignmentWorking = "flag.alignment.working"
FLAGAlignmentDone = "flag.alignment.done"
FLAGQCReportWorking = "flag.qc.report.working"
FLAGQCReportDone = "flag.qc.report.done"
FLAGAllDone = "flag.all.done"

class ClsSample:
    def __init__(self):
        self.strName = ""
        self.strBarcode = ""
        self.strLane = ""
        self.strFastq1 = ""
        self.strFastq2 = ""
        self.strFlagDir = ""
        self.strLogDir = ""
        self.strRootDir = ""
        self.strRawDataDir = ""
        self.strAligner = ""
        self.strRef = ""
    
    def GetUniqueName(self):
        strName = self.strName + "_" + self.strBarcode + "_" + self.strLane
        return strName           
        
    def CreatFolder(self, strDir):
        if not os.path.exists(strDir):
            CMD = "mkdir -p " + strDir
            os.system(CMD)        
    
    def Init(self, strSampleDir, strRootSampleDir, strAlinger, strRef):
        #Get fastq info
        CMD = "find " + strSampleDir + "/ -maxdepth 1 -type f -name \"*.fastq.gz\""
        strFastqList = subprocess.getoutput(CMD)
        vFastqList = strFastqList.split('\n')
        if len(vFastqList) != 2:
            return 1
        else:
            vFastqList.sort()            
            self.strFastq1 = vFastqList[0]
            self.strFastq2 = vFastqList[1]
        
        self.strRawDataDir = strSampleDir
        strDirName = os.path.basename(strSampleDir)
        vItems = strDirName.split("_")
        #Get init info of sample
        self.strLane = vItems[-1]
        self.strBarcode = vItems[-2]
        self.strName = vItems[-3]    
    
        self.strAligner = strAlinger
        self.strRef = strRef
        
        # Create a folder in Sample
        strBackupSampleDir = strRootSampleDir + "/" + strDirName
        self.strRootDir = strBackupSampleDir 
        self.CreatFolder(strBackupSampleDir)
        self.strFlagDir = strBackupSampleDir + "/Flag/" + self.strAligner + "/" + self.strRef         
        self.CreatFolder(self.strFlagDir)
        
        # return value 
        return 0
    
    def Print(self):
        print("Name   :", self.strName)
        print("Fastq 1:", self.strFastq1)
        print("Fastq 2:", self.strFastq2)
    
    def CheckFlags(self, strPhase):    
        if strPhase == "Alignment":
            workingFlag = self.strFlagDir + "/" + FLAGAlignmentWorking
            doneFlag = self.strFlagDir + "/" + FLAGAlignmentDone
            if os.path.exists(doneFlag):
                return 0
            elif os.path.exists(workingFlag):
                return 1
            else:
                CMD = "touch " + workingFlag
                os.system(CMD) 
                return 2        
        if strPhase == "QCReport":
            workingFlag = self.strFlagDir + "/" + FLAGQCReportWorking
            doneFlag = self.strFlagDir + "/" + FLAGQCReportDone
            if os.path.exists(doneFlag):
                return 0
            elif os.path.exists(workingFlag):
                return 1
            else:
                CMD = "touch " + workingFlag
                os.system(CMD) 
                return 2 
            
    # submit jobs for each simple sample -> Go!!
    def SubmitAlignmentJob(self, strProcessedDir, strRootDir, strRootLogDir, strRootQCReportDir, strRootBAMDir):
        # Check working flag 
        iReturn = self.CheckFlags("Alignment")
        if iReturn == 0 or iReturn == 1:  # everything is done or still working
            return 0
        
        # Submit jobs
        self.CreatFolder(strRootLogDir + "/Alignment")
        
        #1: Set std output and std error
        strStdOut = strRootLogDir + "/Alignment/reads_mapping_" + self.GetUniqueName() + ".log.std"
        strStdErr = strRootLogDir + "/Alignment/reads_mapping_" + self.GetUniqueName() + ".log.err"
        if os.path.exists(strStdOut):
            CMD = "rm " + strStdOut
            os.system(CMD)  
        if os.path.exists(strStdErr):
            CMD = "rm " + strStdErr
            os.system(CMD)
                        
        #2: Set other job info
        strJobName = "RMP." + self.GetUniqueName()
        strNumCore = "16"
        
        #３: Get shell script file        
        strCurDirPath = os.path.dirname(os.path.realpath(__file__))
        strShellScript = strCurDirPath + "/reference_mapping_single.sh"
        
        #4: Prepare parameters required by shell script        
        strBamFileBase = self.GetUniqueName() + "_HQ_paired"
        strBamBase = strRootBAMDir + "/" + strBamFileBase
        strRunID = "ProcessedData"
        sampleType = "DNAWholeGenome"
        tmpFolder = TMPFolder + "/" + strRunID
        reference = REFSeq
        refSamIndex = REFSamIndex
        refSamBed = ""
        IN_DIR = strRootDir
        alignerTool = "BWA"
        refVersion = "v38"
        sizeRam = "10g"
        doneAlignment="no"
        doneDedup="no"
        donePhix="no"
        doneUMI="no"
        allLaneDemultiplexOnly = "no"
                
        CMD = ("sbatch " + " --ntasks=" + strNumCore + " --nodes=1 " + 
               "--job-name=" + strJobName + " --output=" + strStdOut + " --error=" + strStdErr + " " + 
               "--export=ALL,ILLUMINA_PROCESSED_DATA_ROOT_DIR=" + strProcessedDir  + " " +
               " " + strShellScript + " " + 
                    "\"" + self.strFastq1 + "\"" + " " + 
                    "\"" + self.strFastq2 + "\"" + " " + 
                    "\"" + strBamBase + "\"" + " " +
                    "\"" + strRunID + "\"" + " " +
                    "\"" + sampleType + "\"" + " " +
                    "\"" + tmpFolder + "\"" + " " +
                    "\"" + strJobName + "\"" + " " +
                    "\"" + reference + "\"" + " " +
                    "\"" + refSamIndex + "\"" + " " +
                    "\"" + refSamBed + "\"" + " " +
                    "\"" + IN_DIR + "\"" + " " +
                    "\"" + alignerTool + "\"" + " " +
                    "\"" + refVersion + "\"" + " " +
                    "\"" + sizeRam + "\"" + " " +
                    "\"" + doneAlignment + "\"" + " " +
                    "\"" + doneDedup + "\"" + " " +
                    "\"" + donePhix + "\"" + " " +
                    "\"" + doneUMI + "\"" + " " +                    
                    "\"" + allLaneDemultiplexOnly + "\"" + " " +
                    "\"" + strNumCore + "\"" + " " +
                    "\"" + self.strFlagDir + "/" + FLAGAlignmentWorking + "\""  + " " +
                    "\"" + self.strFlagDir + "/" + FLAGAlignmentDone + "\"")
        print(CMD)
        os.system(CMD)
    
    def SubmitQCReportJob(self, strProcessedDir, strRootDir, strRootLogDir, strRootQCReportDir, strRootBAMDir):
        # Check working flag 
        iReturn = self.CheckFlags("QCReport")
        if iReturn == 0 or iReturn == 1:  # everything is done or still working
            return 0
        
        # Submit jobs
        self.CreatFolder(strRootLogDir + "/QCReport")
        
        #1: Set std output and std error
        strStdOut = strRootLogDir + "/QCReport/qc_report_" + self.GetUniqueName() + ".log.std"
        strStdErr = strRootLogDir + "/QCReport/qc_report_" + self.GetUniqueName() + ".log.err"
        if os.path.exists(strStdOut):
            CMD = "rm " + strStdOut
            os.system(CMD)  
        if os.path.exists(strStdErr):
            CMD = "rm " + strStdErr
            os.system(CMD)
        
        #2: Set other job info
        strJobName = "QCR." + self.GetUniqueName()
        strNumCore = "1"
        
        #３: Get shell script file        
        strCurDirPath = os.path.dirname(os.path.realpath(__file__))
        strShellScript = strCurDirPath + "/generate_report_single_wgs.sh"
        
        #4: Prepare parameters required by shell script -> Go!!!
        strRunID = "ProcessedData"
        strProjectID = "WeDoNotNeedProjectID"
        strSampleName = self.GetUniqueName()
        strTmpSummaryReportFile = "WeDoNotNeedTmpSummaryReport"
        strBaseSampleDir = os.path.basename(self.strRootDir)
        strLaneNumRunLevel = "1"
        strAlignerTool = "BWA"
        strRefVersion = "v38"
        strLogDir = strRootLogDir + "/QCReport"
        strBAMDir = strRootBAMDir
        strReportDir = strRootQCReportDir
        # if strAlignerTool != "" or strRefVersion != "":
        #     strReportDir += "/" + strAlignerTool + "/" + strRefVersion  
        strSampleFlagDir = self.strFlagDir
        strContainUMI = "no"
        
        CMD = ("sbatch " + " --ntasks=" + strNumCore + " --nodes=1 " + 
               "--job-name=" + strJobName + " --output=" + strStdOut + " --error=" + strStdErr + " " + 
               "--export=ALL,ILLUMINA_PROCESSED_DATA_ROOT_DIR=" + strProcessedDir  + " " +
               " " + strShellScript + " " + 
                    "\"" + strRunID + "\"" + " " + 
                    "\"" + strProjectID+ "\"" + " " + 
                    "\"" + strSampleName + "\"" + " " +
                    "\"" + strTmpSummaryReportFile + "\"" + " " +
                    "\"" + strBaseSampleDir + "\"" + " " +
                    "\"" + strLaneNumRunLevel + "\"" + " " +
                    "\"" + strJobName + "\"" + " " +
                    "\"" + strAlignerTool + "\"" + " " +
                    "\"" + strRefVersion + "\"" + " " +
                    "\"" + strLogDir + "\"" + " " +
                    "\"" + strBAMDir + "\"" + " " +
                    "\"" + strReportDir + "\"" + " " +
                    "\"" + strSampleFlagDir + "\"" + " " +
                    "\"" + strContainUMI + "\"" + " " +
                    "\"" + self.strRawDataDir + "\"" + " " +
                    "\"" + self.strFlagDir + "/" + FLAGQCReportWorking + "\""  + " " +
                    "\"" + self.strFlagDir + "/" + FLAGQCReportDone + "\"")
        print(CMD)
        os.system(CMD)
        
class ClsBuild:
    def __init__(self):
        self.vSample = []
        self.strProcessedDir = ""       
        self.strRootDir = ""
        self.strRootFlagDir = ""
        self.strRootLogDir = ""
        self.strRootQCReportDir = ""
        self.strRootSampleDir = ""
        self.strRootBAMDir = ""
        self.strAligner = "BWA"
        self.strRef = "v38"
                    
    def Init(self, strDir, strOutput):
        self.strProcessedDir = strOutput
         
        #1: Set Root Dir
        self.strRootDir = strOutput + "/" + "ProcessedData"
        self.CreatFolder(self.strRootDir)
        
        self.strRootFlagDir = self.strRootDir + "/Flag/" + self.strAligner + "/" + self.strRef
        self.CreatFolder(self.strRootFlagDir)
         
        self.strRootLogDir = self.strRootDir + "/Log/" + self.strAligner + "/" + self.strRef
        self.CreatFolder(self.strRootLogDir)
        
        self.strRootQCReportDir = self.strRootDir + "/reports/" + self.strAligner + "/" + self.strRef
        self.CreatFolder(self.strRootQCReportDir)
        
        self.strRootSampleDir = self.strRootDir + "/Sample"
        self.CreatFolder(self.strRootSampleDir)
        
        self.strRootBAMDir = self.strRootDir + "/BAM/" + self.strAligner + "/" + self.strRef
        self.CreatFolder(self.strRootBAMDir)
        
        #2: Init Samples
        self.vSample.clear()
        CMD = "find " + strDir + "/ -maxdepth 1 -type d -name \"Sample_*\""
        strList = subprocess.getoutput(CMD)
        vList = strList.split('\n')        
        print(vList)
        for strSampleDir in vList:
            objSample = ClsSample()
            iReturn = objSample.Init(strSampleDir, self.strRootSampleDir, self.strAligner, self.strRef)
            if iReturn == 0:
                self.vSample.append(objSample)
        self.Print()
    
    def CreatFolder(self, strDir):
        if not os.path.exists(strDir):
            CMD = "mkdir -p " + strDir
            os.system(CMD)
        
    def Print(self):
        print("Sample Num:", len(self.vSample))
        iIndex = 1
        for sample in self.vSample:
            print("******", iIndex, "******")
            sample.Print()
            iIndex += 1
            
    def CheckFlags(self, strPhase):
        if strPhase == "GeneralStatus":            
            vDir = os.listdir(self.strRootFlagDir)
            if len(vDir) == 0:
                # send the email to notify people the a new analysis will be launched!
                #it is new -> send start notification
                strInputFile = os.path.dirname(self.strRootDir)
                strPlatform = "COVID 19"     
                strMsg = "============ " + strPlatform + " ============\n"
                strMsg += ("New Customized QC pipeline has been Launched in Biowulf: " + strInputFile) 
                strSubject = strInputFile + " has been Launched (Alignment and QC)"
                CMD = "echo -e \"" + strMsg + "\" | mail -r " + EMAILSender + " -s \"" + strSubject + "\" " + EMAILReceiver            
                print(strMsg)               
                os.system(CMD)
                return 1
            elif FLAGAllDone in vDir:
                print("Current Build is all set! No more action required!")
                return 0
            else:
                #check all done flag ->
                if len(vDir) == 2 and FLAGAlignmentDone in vDir and FLAGQCReportDone in vDir:
                    strInputFile = os.path.dirname(self.strRootDir)
                    # need to send all done email and create all done flag
                    # 1: send email
                    strPlatform = "COVID 19"
                    strFileName = "QCReport" + "-" + self.strAligner + "-" + self.strRef + ".csv"
                    strQCReportFile = self.strRootDir + "/" + strFileName     
                    strMsg = "============ " + strPlatform + " ============\n\n"
                    strMsg += ("New Customized QC pipeline is all set: " + strInputFile + "\n\n" + 
                               "Total Number of Sample: " + str(len(self.vSample)) + "\n\n" + 
                               "QC Report: " + strQCReportFile)
                                          
                    strSubject = strInputFile + " is all set."

                    CMD = ("echo -e \"" + strMsg + "\" | mail -r " + EMAILSender + " -a " + strQCReportFile +
                           " -s \"" + strSubject + "\" " + EMAILReceiver)
                    print(strMsg)               
                    os.system(CMD)
                    # 2: Set all done flag
                    allDoneFlag = self.strRootFlagDir + "/" + FLAGAllDone
                    if not os.path.exists(allDoneFlag):
                        CMD = "touch " + allDoneFlag
                        os.system(CMD)
                    print("Current Build is all set!")
                    return 0
                else:
                    return 2
            
        if strPhase == "Alignment":
            workingFlag = self.strRootFlagDir + "/" + FLAGAlignmentWorking
            doneFlag = self.strRootFlagDir + "/" + FLAGAlignmentDone
            if os.path.exists(doneFlag):
                return 0
            elif os.path.exists(workingFlag):
                # Check if everything is all set --> if it is then done -> Go!!
                # 1: Check the number of sample folders
                CMD = "find " + self.strRootSampleDir + " -maxdepth 1 -type d -name \"Sample_*\" | wc -l"
                iNumSampleFolder = int(subprocess.getoutput(CMD))
                
                strFindPath = self.strRootSampleDir + "/*/Flag" + "/" + self.strAligner + "/" + self.strRef 
                # 2: Check the number of working flags                
                CMD = "find " + strFindPath + " -maxdepth 1 -type f -name \"" + FLAGAlignmentWorking + "\" | wc -l"
                iNumSampleWorking = int(subprocess.getoutput(CMD))
                # 3: Check the number of done flags
                CMD = "find " + strFindPath + " -maxdepth 1 -type f -name \"" + FLAGAlignmentDone + "\" | wc -l"
                iNumSampleDone = int(subprocess.getoutput(CMD))
                # <--
                if iNumSampleWorking == 0 and iNumSampleFolder == iNumSampleDone:
                    # change main level working to done and return 0
                    CMD = "rm " + workingFlag
                    os.system(CMD)
                    CMD = "touch " + doneFlag
                    os.system(CMD)
                    print("The phase of Alignment:", "All set!")
                    print("--> Main level working flag has been updated into done flag!")
                    return 0
                else:
                    return 1
            else:
                CMD = "touch " + workingFlag
                os.system(CMD)
                return 2 # new one             
        if strPhase == "QCReport":
            workingFlag = self.strRootFlagDir + "/" + FLAGQCReportWorking
            doneFlag = self.strRootFlagDir + "/" + FLAGQCReportDone
            if os.path.exists(doneFlag):
                return 0
            elif os.path.exists(workingFlag):
                # Check if everything is all set --> if it is then done -> Go!!
                # 1: Check the number of sample folders
                CMD = "find " + self.strRootSampleDir + " -maxdepth 1 -type d -name \"Sample_*\" | wc -l"
                iNumSampleFolder = int(subprocess.getoutput(CMD))
                
                strFindPath = self.strRootSampleDir + "/*/Flag" + "/" + self.strAligner + "/" + self.strRef 
                # 2: Check the number of working flags
                CMD = "find " + strFindPath + " -maxdepth 1 -type f -name \"" + FLAGQCReportWorking + "\" | wc -l"
                iNumSampleWorking = int(subprocess.getoutput(CMD))
                # 3: Check the number of done flags
                CMD = "find " + strFindPath + " -maxdepth 1 -type f -name \"" + FLAGQCReportDone + "\" | wc -l"
                iNumSampleDone = int(subprocess.getoutput(CMD))
                # <--
                if iNumSampleWorking == 0 and iNumSampleFolder == iNumSampleDone:
                    # change main level working to done and return 0                    
                    # --> Generate final QC report
                    iReturn = self.GenerateQCSummary()
                    if iReturn == 0:
                        CMD = "rm " + workingFlag
                        os.system(CMD)
                        CMD = "touch " + doneFlag
                        os.system(CMD)
                        print("The phase of QC Report:", "All set!")
                        print("--> Main level working flag has been updated into done flag!")
                        return 0
                    else:
                        return 1
                else:
                    return 1
            else:
                CMD = "touch " + workingFlag
                os.system(CMD)
                return 2 # new one       
            
    def CustomizedAlignment(self):
        iReturn = self.CheckFlags("Alignment")
        if iReturn == 0 or iReturn == 1:  # everything is done or still working
            return iReturn
        
        # This is the new case that flowcell need a new analysis--> Go!!
        for sample in self.vSample:
            sample.SubmitAlignmentJob(self.strProcessedDir, self.strRootDir, self.strRootLogDir, self.strRootQCReportDir, self.strRootBAMDir)
        return 2
    
    def CustomizedQCReport(self):
        iReturn = self.CheckFlags("Alignment")
        if iReturn != 0:
            return 3
        # Only do the logic below if alignment has been finished !
        iReturn = self.CheckFlags("QCReport")
        if iReturn == 0 or iReturn == 1:  # everything is done or still working
            return iReturn
         
        # This is the new case that flowcell need a new analysis (QC report) --> Go!!
        for sample in self.vSample:
            sample.SubmitQCReportJob(self.strProcessedDir, self.strRootDir, self.strRootLogDir, self.strRootQCReportDir, self.strRootBAMDir)
        return 2
                
    def GenerateQCSummary(self):
        print("I am the last step: GenerateQCSummary!")
        strFileName = "QCReport" + "-" + self.strAligner + "-" + self.strRef + ".csv"
        strQCReport = self.strRootDir + "/" + strFileName
        
        # remove the old one
        if os.path.exists(strQCReport):
            CMD = "rm " + strQCReport
            os.system(CMD)
        
        #Create new file  
        CMD = ("echo -e \"" + 
       "Sample,Barcode,Lane,Project,\
Sequenced Reads,HQ Reads,Mapped HQ Reads,Mapped HQ & Dedup Reads,\
PCR Dup (Paired),% PCR Dup,Optical Dup (Paired),% Optical Dup,% All Dups,\
HQ Mapped (Dedup and phiX Removed),\
On-target Reads,Capture Kit Efficiency (% On-target/HQ Mapped),\
Overall Sequencing Efficiency (% On-target/All Sequenced),Median InsertSize,Mean InsertSize,\
Capture Kit Bases Covered,% Capture Kit Bases Covered,Avg. Coverage,Capture Kit Bases >= ${MIN_COVERAGE_THRESHOLD1}X,\
% Capture Kit Bases >=${MIN_COVERAGE_THRESHOLD1}X,Capture Kit Bases >= ${MIN_COVERAGE_THRESHOLD2}X,\
% Capture Kit Bases >=${MIN_COVERAGE_THRESHOLD2}X,UCSC CDS Bases Covered,% UCSC CDS Bases Covered,\
Avg. Coverage,UCSC CDS Bases >= ${MIN_COVERAGE_THRESHOLD1}X,% UCSC CDS Bases >=${MIN_COVERAGE_THRESHOLD1}X,\
UCSC CDS Bases >= ${MIN_COVERAGE_THRESHOLD2}X,% UCSC CDS Bases >=${MIN_COVERAGE_THRESHOLD2}X,\
Attention,Percentage of trimmed bases R1,Percentage of trimmed bases R2,\
Percentage of inserts <36 bp R1,Percentage of inserts <36 bp R2,% of the lane,% binned indices,\
Dimer Reads,% Dimer\"" + 
        " >> " + strQCReport)
        
        print(CMD)
        os.system(CMD)
        # find all report log files 
        strQCLogDir = self.strRootLogDir + "/QCReport"
        CMD = "find " + strQCLogDir + " -maxdepth 1 -type f -iname " + "\"*.log.std\""
        print(CMD)
        strLogList = subprocess.getoutput(CMD)
        vLogList = strLogList.split('\n')
        for strLogFile in vLogList:
            CMD = "grep -A 1 \"Contents to be streamed to the flowcell-level report file\" " + strLogFile + " | tail -n 1 >> " + strQCReport   
            os.system(CMD)
        return 0
        
def main():
    # Folder Path
    strDir = sys.argv[1]
    strOutputDir = sys.argv[2]
    
    # Init build -> Done
    objBuild = ClsBuild()
    objBuild.Init(strDir, strOutput)
        
    # Sent email start notification
    iReturn = objBuild.CheckFlags("GeneralStatus")
    if iReturn == 0:
        return 0 
    
    # Spawn jobs do customized alignment
    objBuild.CustomizedAlignment()    
    objBuild.CustomizedQCReport()
    return 0
            
if __name__ == "__main__":
    main() 
        