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
from datetime import datetime
TMPFolder = "/scratch/lix33/TmpSequencing/"
REFSeq = "/data/COVID_WGS/lix33/DCEG/CGF/Bioinformatics/Production/data/ref38/Homo_sapiens_assembly38.fasta"
REFSamIndex = "/data/COVID_WGS/lix33/DCEG/CGF/Bioinformatics/Production/data/ref38/Homo_sapiens_assembly38.fasta.fai"

EMAILSender = "xin.li4@nih.gov"
EMAILReceiver = "xin.li4@nih.gov"
#EMAILReceiver = "xin.li4@nih.gov,mingyi.wang@nih.gov"

FLAGAlignmentWorking = "flag.alignment.working"
FLAGAlignmentDone = "flag.alignment.done"
FLAGQCReportWorking = "flag.qc.report.working"
FLAGQCReportDone = "flag.qc.report.done"
FLAGMergeSampleWorking = "flag.merge.sample.working"
FLAGMergeSampleDone = "flag.merge.sample.done"
FLAGMultiSample = "flag.multi.sample"
FLAGAllDone = "flag.all.done"

FLAGMultiSample = "flag.multi.sample"

class ClsSample:
    def __init__(self):
        self.strName = ""
        self.strBarcode = ""
        self.strLane = ""
        self.strFastq1 = ""
        self.strFastq2 = ""
        self.strFlagDir = ""
        self.strLogDir = ""
        self.strSampleDir = ""
        self.strFlowcellDir = ""
        self.strAligner = ""
        self.strRef = ""
        self.iMergeSample = 0
    
    def GetUniqueName(self, strSplit="-"):
        strName = self.strName + strSplit + self.strBarcode + strSplit + self.strLane
        return strName           
        
    def CreatFolder(self, strDir):
        if not os.path.exists(strDir):
            CMD = "mkdir -p " + strDir
            os.system(CMD)        
    
    def Init(self, strSampleDir, strAlinger, strRef):
        #Init Sample Dir 
        self.strSampleDir = strSampleDir 
        #print("strSampleDir", strSampleDir)        
                
        # #Get fastq info
        # CMD = "find " + strSampleDir + " -maxdepth 1 -type l -name \"*.fastq.gz\""
        # #print(CMD)
        # strFastqList = subprocess.getoutput(CMD)
        # vFastqList = strFastqList.split('\n')
        # #print("vFastqList:", vFastqList)
        # if len(vFastqList) != 2:            
        #     return 1
        # else:
        #     vFastqList.sort()            
        #     self.strFastq1 = vFastqList[0]
        #     self.strFastq2 = vFastqList[1]
                
        strDirName = os.path.basename(strSampleDir)
        vItems = strDirName.split("_")
        #Get init info of sample (example: Sample_Name_Barcode_Lane)
        self.strLane = vItems[-1]
        self.strBarcode = vItems[-2]
        self.strName = vItems[-3]    
    
        self.strAligner = strAlinger
        self.strRef = strRef
        
        # Create a folder in Sample        
        self.strFlagDir = strSampleDir + "/Flag/" + self.strAligner + "/" + self.strRef         
        self.CreatFolder(self.strFlagDir)
        
        # return value 
        return 0

    def UpdateFastqInfo(self):
        # Set Fastq Info after the step of merging sample been finished! ->
        # Check if we have the special tag ->  
        strSpecialTag = self.strFlagDir + "/" + FLAGMultiSample
        CMD = ""        
        if os.path.exists(strSpecialTag): 
            # This is the case of merged sample (Check file)
            CMD = "find " + self.strSampleDir + " -maxdepth 1 -type f -name \"*.fastq.gz\""
            self.iMergeSample = 1 # True
        else:
            # Get fastq info (Check soft link)
            CMD = "find " + self.strSampleDir + " -maxdepth 1 -type l -name \"*.fastq.gz\""
            self.iMergeSample = 0 # False
        #print(CMD)
        strFastqList = subprocess.getoutput(CMD)
        vFastqList = []
        if strFastqList != "":
            vFastqList = strFastqList.split('\n')
        #print("vFastqList:", vFastqList)
        if len(vFastqList) != 2:            
            return 1, vFastqList
        else:
            vFastqList.sort()            
            self.strFastq1 = vFastqList[0]
            self.strFastq2 = vFastqList[1]
            return 0, vFastqList
        # <- 
    
    def Print(self):
        print("Name   :", self.strName)
        print("Fastq 1:", self.strFastq1)
        print("Fastq 2:", self.strFastq2)
    
    def CheckFlags(self, strPhase, bTouch=True):
        if strPhase == "MergeSample":
            workingFlag = self.strFlagDir + "/" + FLAGMergeSampleWorking
            doneFlag = self.strFlagDir + "/" + FLAGMergeSampleDone
            if os.path.exists(doneFlag):
                return 0
            elif os.path.exists(workingFlag):
                return 1
            else:
                CMD = "touch " + workingFlag
                os.system(CMD) 
                return 2 
            
        if strPhase == "Alignment":
            workingFlag = self.strFlagDir + "/" + FLAGAlignmentWorking
            doneFlag = self.strFlagDir + "/" + FLAGAlignmentDone
            if os.path.exists(doneFlag):
                return 0
            elif os.path.exists(workingFlag):
                return 1
            else:
                if bTouch:
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
            
    def SubmitMergeSampleJob(self, strProcessedDir, strFlowcellDir, strRootLogDir, strRootQCReportDir, strRootBAMDir):
        iReturn = self.CheckFlags("MergeSample")
        if iReturn == 0 or iReturn == 1:  # everything is done or still working
            return 0
        
        # Submit jobs
        self.CreatFolder(strRootLogDir + "/MergeSample")
        
        #1: Set std output and std error
        strStdOut = strRootLogDir + "/MergeSample/merge_sample_" + self.GetUniqueName() + ".log.std"
        strStdErr = strRootLogDir + "/MergeSample/merge_sample_" + self.GetUniqueName() + ".log.err"
        if os.path.exists(strStdOut):
            CMD = "rm " + strStdOut
            os.system(CMD)  
        if os.path.exists(strStdErr):
            CMD = "rm " + strStdErr
            os.system(CMD)
            
        # Get Script file
        strCurDirPath = os.path.dirname(os.path.realpath(__file__))
        strPythonScript = strCurDirPath + "/MergeSampleSingle.py"
        strRunningTime = "10-00:00:00"
        strNumCore = "1"
        strJobName = "MSP." + self.GetUniqueName() # Merge Sample
        
        #submit jobs: we only need sample dir and flag dir -> Go!
        CMD = ("sbatch " + " --ntasks=" + strNumCore + " --nodes=1 " + 
               "--job-name=" + strJobName + " --time=" + strRunningTime + " --output=" + strStdOut + " --error=" + strStdErr + 
               " --wrap=\"" + "python3" + " " + strPythonScript +  
               " '" + self.strSampleDir + "'" +  
               " '" + self.strFlagDir + "/" + FLAGMergeSampleWorking + "'" + 
               " '" + self.strFlagDir + "/" + FLAGMergeSampleDone + "'" +  "\"")
        
        # CMD = ("python3" + " " + strPythonScript +  
        #        " '" + self.strSampleDir + "'" +  
        #        " '" + self.strFlagDir + "/" + FLAGMergeSampleWorking + "'" + 
        #        " '" + self.strFlagDir + "/" + FLAGMergeSampleDone + "'")
        print(CMD)
        os.system(CMD)
            
    # submit jobs for each simple sample -> Go!!
    def SubmitAlignmentJob(self, strProcessedDir, strFlowcellDir, strRootLogDir, strRootQCReportDir, strRootBAMDir):
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
        strNumCore = "32"
        
        #３: Get shell script file        
        strCurDirPath = os.path.dirname(os.path.realpath(__file__))
        strShellScript = strCurDirPath + "/reference_mapping_single.sh"
        
        #4: Prepare parameters required by shell script        
        strBamFileBase = self.GetUniqueName("_") + "_HQ_paired"
        strBamBase = strRootBAMDir + "/" + strBamFileBase
        strRunID = os.path.basename(strFlowcellDir)
        sampleType = "DNAWholeGenome"
        tmpFolder = TMPFolder + "/" + strRunID
        reference = REFSeq
        refSamIndex = REFSamIndex
        refSamBed = ""
        IN_DIR = strFlowcellDir
        alignerTool = "BWA"
        refVersion = "v38"
        sizeRam = "10g"
        doneAlignment="no"
        doneDedup="no"
        donePhix="no"
        doneUMI="no"
        allLaneDemultiplexOnly = "no"
        strRunningTime = "10-00:00:00"
                
        CMD = ("sbatch " + " --ntasks=" + strNumCore + " --nodes=1 " + 
               "--job-name=" + strJobName + " --time=" + strRunningTime + " --output=" + strStdOut + " --error=" + strStdErr + " " + 
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
                    "\"" + self.strFlagDir + "/" + FLAGAlignmentDone + "\""  + " " +
                    "\"" + str(self.iMergeSample) + "\"")
        print(CMD)
        os.system(CMD)
    
    def SubmitQCReportJob(self, strProcessedDir, strFlowcellDir, strRootLogDir, strRootQCReportDir, strRootBAMDir):
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
        strRunID = os.path.basename(strFlowcellDir)
        strProjectID = os.path.basename(os.path.dirname(self.strSampleDir))
        strSampleName = self.GetUniqueName("_")
        strTmpSummaryReportFile = "WeDoNotNeedTmpSummaryReport"
        strBaseSampleDir = os.path.basename(self.strSampleDir)
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
        strRunningTime = "10-00:00:00"
        
        CMD = ("sbatch " + " --ntasks=" + strNumCore + " --nodes=1 " + 
               "--job-name=" + strJobName + " --time=" + strRunningTime + " --output=" + strStdOut + " --error=" + strStdErr + " " + 
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
                    "\"" + self.strSampleDir + "\"" + " " +
                    "\"" + self.strFlagDir + "/" + FLAGQCReportWorking + "\""  + " " +
                    "\"" + self.strFlagDir + "/" + FLAGQCReportDone + "\"")
        print(CMD)
        os.system(CMD)
        
# class ClsSubject:
#     def __init__(self):
#         self.sample = ClsSample()
#         self.strName = ""
#         self.strBarcode = ""
#         self.
#
        
class ClsBuild:
    def __init__(self):
        self.vSample = []
        self.strProcessedDir = ""       
        self.strFlowcellDir = "" # It is flowcell Dir
        self.strRootFlagDir = ""
        self.strRootLogDir = ""
        self.strRootQCReportDir = ""
        #self.strRootSampleDir = ""
        self.strRootBAMDir = ""
        self.strAligner = "BWA"
        self.strRef = "v38"
        self.strFlowcellName = ""
        self.bAllSet = False
                    
    def Init(self, strFlowcellDir):
        # Get main processed folder name and file basename
        self.strProcessedDir = os.path.dirname(strFlowcellDir)
        self.strFlowcellName = os.path.basename(strFlowcellDir)
         
        #1: Set Root Dir
        self.strFlowcellDir = strFlowcellDir
        self.CreatFolder(self.strFlowcellDir)
        
        self.strRootFlagDir = self.strFlowcellDir + "/Flag/" + self.strAligner + "/" + self.strRef
        self.CreatFolder(self.strRootFlagDir)
         
        self.strRootLogDir = self.strFlowcellDir + "/Log/" + self.strAligner + "/" + self.strRef
        self.CreatFolder(self.strRootLogDir)
        
        self.strRootQCReportDir = self.strFlowcellDir + "/reports/" + self.strAligner + "/" + self.strRef
        self.CreatFolder(self.strRootQCReportDir)
        
        # self.strRootSampleDir = self.strRootDir + "/Sample"
        # self.CreatFolder(self.strRootSampleDir)
        
        self.strRootBAMDir = self.strFlowcellDir + "/BAM/" + self.strAligner + "/" + self.strRef
        self.CreatFolder(self.strRootBAMDir)
        
        #2: Init Samples
        self.vSample.clear()
        # Get Sample Dir
        
        CMD = "find " + strFlowcellDir + "/CASAVA -maxdepth 3 -type d -name \"Sample_*\""
        strList = subprocess.getoutput(CMD)
        vList = strList.split('\n')        
        #print(vList)
        for strSampleDir in vList:
            objSample = ClsSample()
            iReturn = objSample.Init(strSampleDir, self.strAligner, self.strRef)
            if iReturn == 0:
                self.vSample.append(objSample)
        self.Print()
    
    def CreatFolder(self, strDir):
        if not os.path.exists(strDir):
            CMD = "mkdir -p " + strDir
            os.system(CMD)

    def Print(self):
        print("Sample Num:", len(self.vSample))
        # iIndex = 1
        # for sample in self.vSample:
        #     print("******", iIndex, "******")
        #     sample.Print()
        #     iIndex += 1
            
    def CheckFlags(self, strPhase):
        if strPhase == "GeneralStatus":            
            vDir = os.listdir(self.strRootFlagDir)
            if len(vDir) == 0:
                # send the email to notify people the a new analysis will be launched!
                #it is new -> send start notification
                strInputFile = os.path.dirname(self.strFlowcellDir)
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
                if len(vDir) == 3 and FLAGAlignmentDone in vDir and FLAGQCReportDone in vDir and FLAGMergeSampleDone in vDir:                    
                    # 1: Send successful notification email
                    strInputFile = os.path.dirname(self.strFlowcellDir)
                    # need to send all done email and create all done flag
                    #  send email
                    strPlatform = "COVID 19"
                    strFileName = "QCReport" + "-" + self.strAligner + "-" + self.strRef + ".csv"
                    strQCReportFile = self.strFlowcellDir + "/" + strFileName     
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
        
        if strPhase == "MergeSample":
            workingFlag = self.strRootFlagDir + "/" + FLAGMergeSampleWorking
            doneFlag = self.strRootFlagDir + "/" + FLAGMergeSampleDone
            if os.path.exists(doneFlag):
                return 0
            elif os.path.exists(workingFlag):
                # Check if everything is all set --> if it is then done -> Go!!
                # 1: Check the number of sample folders
                CMD = "find " + self.strFlowcellDir + "/CASAVA -maxdepth 3 -type d -name \"Sample_*\" | wc -l"
                iNumSampleFolder = int(subprocess.getoutput(CMD))
                
                strFindPath = self.strFlowcellDir + "/CASAVA/*/*/*/Flag/" + self.strAligner + "/" + self.strRef 
                # 2: Check the number of working flags                
                CMD = "find " + strFindPath + " -maxdepth 1 -type f -name \"" + FLAGMergeSampleWorking + "\" | wc -l"
                iNumSampleWorking = int(subprocess.getoutput(CMD))
                # 3: Check the number of done flags
                CMD = "find " + strFindPath + " -maxdepth 1 -type f -name \"" + FLAGMergeSampleDone + "\" | wc -l"
                iNumSampleDone = int(subprocess.getoutput(CMD))
                # <--
                if iNumSampleWorking == 0 and iNumSampleFolder == iNumSampleDone:
                    # change main level working to done and return 0
                    CMD = "rm " + workingFlag
                    os.system(CMD)
                    CMD = "touch " + doneFlag
                    os.system(CMD)
                    print("The phase of Merge Sample:", "All set!")
                    print("--> Main level working flag (Merge Sample) has been updated into done flag!")
                    return 0
                else:
                    return 1
            else:
                CMD = "touch " + workingFlag
                os.system(CMD)
                return 2 # new one
            
        if strPhase == "Alignment":
            #Confirm if MergeSample has done            
            doneFlagMergeSample = self.strRootFlagDir + "/" + FLAGMergeSampleDone
            if not os.path.exists(doneFlagMergeSample):
                return 3
                        
            workingFlag = self.strRootFlagDir + "/" + FLAGAlignmentWorking
            doneFlag = self.strRootFlagDir + "/" + FLAGAlignmentDone
            if os.path.exists(doneFlag):
                return 0
            elif os.path.exists(workingFlag):
                # Check if everything is all set --> if it is then done -> Go!!
                # 1: Check the number of sample folders
                CMD = "find " + self.strFlowcellDir + "/CASAVA -maxdepth 3 -type d -name \"Sample_*\" | wc -l"
                iNumSampleFolder = int(subprocess.getoutput(CMD))
                
                strFindPath = self.strFlowcellDir + "/CASAVA/*/*/*/Flag/" + self.strAligner + "/" + self.strRef 
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
            #Confirm if MergeSample and Alignment has done            
            doneFlagMergeSample = self.strRootFlagDir + "/" + FLAGMergeSampleDone
            doneFlagAlignment = self.strRootFlagDir + "/" + FLAGAlignmentDone
            if not os.path.exists(doneFlagMergeSample) or not os.path.exists(doneFlagAlignment):
                return 3
            
            workingFlag = self.strRootFlagDir + "/" + FLAGQCReportWorking
            doneFlag = self.strRootFlagDir + "/" + FLAGQCReportDone
            if os.path.exists(doneFlag):
                return 0
            elif os.path.exists(workingFlag):
                # Check if everything is all set --> if it is then done -> Go!!
                # 1: Check the number of sample folders
                CMD = "find " + self.strFlowcellDir + "/CASAVA -maxdepth 3 -type d -name \"Sample_*\" | wc -l"
                iNumSampleFolder = int(subprocess.getoutput(CMD))
                
                strFindPath = self.strFlowcellDir + "/CASAVA/*/*/*/Flag/" + self.strAligner + "/" + self.strRef 
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
    
    def CustomizedMergeSample(self):
        iReturn = self.CheckFlags("MergeSample")
        if iReturn == 0 or iReturn == 1:  # everything is done or still working
            return iReturn
        
        # This is the new case that flowcell need a new analysis--> Go!!
        for sample in self.vSample:
            sample.SubmitMergeSampleJob(self.strProcessedDir, self.strFlowcellDir, self.strRootLogDir, self.strRootQCReportDir, self.strRootBAMDir)
        return 2
            
    def CustomizedAlignment(self):        
        iReturn = self.CheckFlags("Alignment")
        #print(iReturn, self.strFlowcellDir)
        if iReturn == 0 or iReturn == 1 or iReturn == 3:  # everything is done or still working            
            return iReturn
        
        # This is the new case that flowcell need a new analysis--> Go!!
        for sample in self.vSample:
            # Check if the alignment has been done (the sum-fastq file will be removed automatically once aligment is all set)
            iReturn = sample.CheckFlags("Alignment", False)
            #print(sample.strName, sample.strBarcode, iReturn)
            if iReturn == 0 or iReturn == 1:  # everything is done or still working
                continue
            
            iResult, vFastqList = sample.UpdateFastqInfo()
            #print(vFastqList)
            if iResult != 0:
                print("Error: Bad number of fastq files!")
                print("vFastqList Len:", len(vFastqList))
                print("vFastqList    :", vFastqList)                
                return 4
            sample.SubmitAlignmentJob(self.strProcessedDir, self.strFlowcellDir, self.strRootLogDir, self.strRootQCReportDir, self.strRootBAMDir)
        return 2
    
    def CustomizedQCReport(self):
        # Only do the logic below if alignment has been finished !
        iReturn = self.CheckFlags("QCReport")
        if iReturn == 0 or iReturn == 1 or iReturn == 3:  # everything is done or still working
            return iReturn
         
        # This is the new case that flowcell need a new analysis (QC report) --> Go!!
        for sample in self.vSample:
            sample.SubmitQCReportJob(self.strProcessedDir, self.strFlowcellDir, self.strRootLogDir, self.strRootQCReportDir, self.strRootBAMDir)
        return 2
                
    def GenerateQCSummary(self):
        print("I am the last step: GenerateQCSummary!")
        strFileName = "QCReport" + "-" + self.strAligner + "-" + self.strRef + ".csv"
        strQCReport = self.strFlowcellDir + "/" + strFileName
        
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
Capture Kit Bases Covered,% Capture Kit Bases Covered,Avg. Coverage,Capture Kit Bases >= 10X,\
% Capture Kit Bases >=10X,Capture Kit Bases >= 15X,\
% Capture Kit Bases >=15X,UCSC CDS Bases Covered,% UCSC CDS Bases Covered,\
Avg. Coverage,UCSC CDS Bases >= 10X,% UCSC CDS Bases >=10X,\
UCSC CDS Bases >= 15X,% UCSC CDS Bases >=15X,\
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
    print('\n', "==============")
    today = datetime.today()
    print("Date:", today.strftime("%B %d, %Y"))
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print("Time:", dt_string)
    print("==============", '\n')
    # Folder Path
    strPrcessedDataDir = sys.argv[1]    
    
    # Collect All flowcells in current Folder (each flowcell is a build) -->
    vFlowcellDir = []
    for subDir in os.listdir(strPrcessedDataDir):
        strSubDirPath = strPrcessedDataDir + "/" + subDir
        if os.path.isdir(strSubDirPath):
            vFlowcellDir.append(strSubDirPath)
    # make sure the order are always the same for each run 
    vFlowcellDir.sort()
    # <-- 
    #print(vFlowcellDir)
        
    #Create Builds
    vBuild = []
    iNum = 0
    iMaxSimuFlowcell = 4
    for strFlowcellDir in vFlowcellDir:
        print()
        objBuild = ClsBuild()
        objBuild.Init(strFlowcellDir)
            
        # Sent email start notification
        iReturn = objBuild.CheckFlags("GeneralStatus")
        if iReturn == 0:
            objBuild.bAllSet = True
            vBuild.append(objBuild)
            continue
        
        # Spawn jobs do customized alignment
        objBuild.CustomizedMergeSample()
        objBuild.CustomizedAlignment()    
        objBuild.CustomizedQCReport()
        vBuild.append(objBuild)
        
        # Only do one (the first) unfinished flowcell for each run
        # We can do 3 flowcells at the same time to save the total running time
        iNum += 1
        if iNum >= iMaxSimuFlowcell:             
            break
    
    #Print status of each build
    print("\n", ">>>>> Summary (All Finished Flowcells + the first working one) <<<<<")
    for build in vBuild:
        print("Flowcell:", build.strFlowcellName, "-> All Set Status:", build.bAllSet)
    
    return 0
            
if __name__ == "__main__":
    main() 
        