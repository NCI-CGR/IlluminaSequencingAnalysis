'''
Input: keytable
e.g.
/data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Keytable/low_input/sub_keytable/low_input_01_50_keytable.csv

Output:
the running results after 2nd pipeline be finished!
'''

import sys
import os
import subprocess

DIRCustomizedQC = "/home/lix33/lxwg/Git/IlluminaSequencingAnalysis/CustomizedQC/SourceCode"
DIR2ndPipeline = "/home/lix33/lxwg/Git/IlluminaSequencingAnalysis/COVID_2nd_Pipeline/SourceCode"

DIRBAMRoot = "/data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Data/BAM/Batch"

DIR2ndPipelineTMP = "/data/COVID_WGS/lix33/Test/2ndpipeline/Build/tmp"

DIRBuildProcess = "/data/COVID_WGS/lix33/Test/2ndpipeline/Build/processed"
DIRBuffer = "/data/COVID_WGS/lix33/Test/2ndpipeline/Data/secondary_buf"
DIRBAMReformat = "/data/COVID_WGS/lix33/Test/2ndpipeline/Data/BAM_reformatted"

# Email notification List
# EMAILReceiverEnd = "xin.li4@nih.gov,kristine.jones@nih.gov,dongjing.wu@nih.gov,wen.luo@nih.gov,hicksbel@mail.nih.gov"
EMAILReceiverEnd = "xin.li4@nih.gov"
EMAILSender = "xin.li4@nih.gov"

DONEFlagEmail = "flag_email_sent_done"
DONEFlagWarningEmail = "flag_email_warning_sent_done"


# --> Defind the data structure to save the subject info from Keytable
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
        # print(vItem)
        self.strFlowcellID = vItem[0].split('_')[-1][1:]
        self.strUSUID = vItem[1]
        self.strCGRID = vItem[2].split('-')[-1]
        if len(vItem) == 4 and vItem[3] == "topoff":
            self.bTopoff = True


class ClsSubject:

    def __init__(self):
        self.strCGRID = ""
        self.strAnalysisID = ""
        self.vSample = []
        self.bHit = False
        
    def GetAnalysisID(self):
        if len(self.vSample) == 0:
            return ""
        strAnalysisID = self.vSample[0].strCGRID
        for sample in self.vSample:
            strAnalysisID += "_" + sample.strUSUID
        self.strAnalysisID = strAnalysisID 
        return  strAnalysisID


class ClsBuild:

    def __init__(self):
        self.vSubject = []
        self.strRootDir = ""
        self.strManifestFile = ""
        self.strCurBAMDir = ""
    
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
        
        # Update Analysis ID
        for subject in self.vSubject:
            subject.GetAnalysisID()
            
        # Update strBAMDir
        self.strCurBAMDir = DIRBAMRoot + "/" + os.path.basename(strKeyTable).split('.')[0]
# <--


def MergeBAM(strKeytable, iSubjectNum, strKTName, objBuild): 
    print("iSubjectNum      :", iSubjectNum)
    # Pre check if everything is all set! -->
    strFlagDir = DIRBAMRoot + "/" + strKTName + "/Flag/MergeSubject"
    # bCallScript = True    
    if os.path.exists(strFlagDir):
        CMD = "find " + strFlagDir + " -iname '*.done' | wc -l"
        print("CMD:", CMD)
        iDoneNum = int(subprocess.getoutput(CMD))
        print("iDoneNum     :", iDoneNum)
        
        CMD = "find " + strFlagDir + " -iname '*.flag.*' | wc -l"
        print("CMD:", CMD)
        iTotalFlagNum = int(subprocess.getoutput(CMD))
        print("iTotalFlagNum:", iTotalFlagNum)
        
        if iDoneNum == iSubjectNum:
            print("MergeBAM is All Set!")
            return 0
        
        if iTotalFlagNum == iSubjectNum:
            print("All jobs in the phase of merge BAM are still running!")
            return 1
        else:
            print("****** Merge Project ********")
            print("Rerun some subjects with missing flags!")
            print("Re-run Subject Num:", iSubjectNum - iTotalFlagNum)
            print("****** ************* ********")
            # Output the file that missing flags -> Go!
            vFlag = subprocess.getoutput("find " + strFlagDir + " -iname '*.flag.*'").split('\n')
            print("Size objBuild.vSubject:", len(objBuild.vSubject))
            for objSubject in objBuild.vSubject:
                for strFlag in vFlag:
                    bFind = False                
                    if(objSubject.vSample[0].strUSUID in strFlag and 
                       objSubject.vSample[0].strCGRID in strFlag):
                        objSubject.bHit = True
                        bFind = True
                        break                        
                if not bFind:
                    print("Warning: Flag does not find ->", objSubject.vSample[0].strUSUID, "-->", objSubject.vSample[0].strCGRID)  
            
    # if not, go ahead to call 
    print("\n", "Run: MergeSubject.py -->", "\n")
    strScript = DIRCustomizedQC + "/Tools/MergeSubject/MergeSubject.py"
    CMD = "python " + strScript + " " + strKeytable
    os.system(CMD)
    return 1
    # <-- 


def ReconstructBAMDir(strKTName):
    # Create working flag and done flag
    strFlagDir = DIRBAMRoot + "/" + strKTName + "/Flag/ReconstructBAMDir"
    strLogDir = DIRBAMRoot + "/" + strKTName + "/Log/ReconstructBAMDir"
    
    if not os.path.exists(strFlagDir):
        CMD = "mkdir -p " + strFlagDir
        os.system(CMD)
        
    if not os.path.exists(strLogDir):
        CMD = "mkdir -p " + strLogDir
        os.system(CMD)
        
    strFlagWorking = strFlagDir + "/reconstruct_BAM_Dir.flag.working"
    strFlagDone = strFlagDir + "/reconstruct_BAM_Dir.flag.done"
    
    if os.path.exists(strFlagWorking):
        print("ReconstructBAMDir is still working!")
        return 1
    
    if os.path.exists(strFlagDone):
        print("ReconstructBAMDir has been finished!")
        return 0

    # touch working flag 
    CMD = "touch " + strFlagWorking
    os.system(CMD)
    
    strScriptBash = "bash " + DIR2ndPipeline + "/step7b_take_incoming_bams.sh" + " " + strFlagWorking + " " + strFlagDone
    
    # Submit jobs ->
    strNumCore = "8"
    strNodeNum = "1"
    strJobName = "RBD-" + strKTName
    strRunningTime = "10-00:00:00"
    strStdOut = strLogDir + "/reconstruct_BAM_Dir.std.out"
    strStdErr = strLogDir + "/reconstruct_BAM_Dir.std.err"
    CMD = ("sbatch " + "--ntasks=" + strNumCore + " " + 
                       "--nodes=" + strNodeNum + " " + 
                       "--job-name=" + strJobName + " " + 
                       "--time=" + strRunningTime + " " + 
                       "--output=" + strStdOut + " " + 
                       "--error=" + strStdErr + " " + 
                       "--wrap=\"" + strScriptBash + "\"")
    print(CMD)
    os.system(CMD)
    return 1                                        


# Step 3
def RecalibrateBAMFile(iSubjectNum, strKTName, objBuild):
    print("iSubjectNum      :", iSubjectNum)
    strManifestFile = DIR2ndPipelineTMP + "/" + strKTName + "/Manifest_Mimic.csv"
    # Pre check if everything is all set! -->
    strFlagDir = DIRBAMRoot + "/" + strKTName + "/Flag/RecalibrateBAM"    
    if os.path.exists(strFlagDir):
        CMD = "find " + strFlagDir + " -iname '*.done' | wc -l"
        iDoneNum = int(subprocess.getoutput(CMD))
        print("iDoneNum     :", iDoneNum)
        
        CMD = "find " + strFlagDir + " -iname '*.flag.*' | wc -l"
        iTotalFlagNum = int(subprocess.getoutput(CMD))
        print("iTotalFlagNum:", iTotalFlagNum)
        
        if iDoneNum == iSubjectNum:
            print("RecalibrateBAM is All Set!")
            return 0
        
        if iTotalFlagNum == iSubjectNum:
            print("All jobs in the phase of merge BAM are still running!")
            return 1
        else:
            print("****** Recalibrate BAM File ********")
            print("Rerun some subjects with missing flags!")
            print("Re-run Subject Num:", iSubjectNum - iTotalFlagNum)
            print("****** ************* ********")
            # Output the file that missing flags -> Go!
            vFlag = subprocess.getoutput("find " + strFlagDir + " -iname '*.flag.*'").split('\n')
            print("Size objBuild.vSubject:", len(objBuild.vSubject))
            for objSubject in objBuild.vSubject:
                for strFlag in vFlag:
                    bFind = False                
                    if(objSubject.vSample[0].strUSUID in strFlag and 
                       objSubject.vSample[0].strCGRID in strFlag):
                        objSubject.bHit = True
                        bFind = True
                        break                        
                if not bFind:
                    print("Warning: Flag does not find ->", objSubject.vSample[0].strUSUID, "-->", objSubject.vSample[0].strCGRID)            
    
    # if not, go ahead to call 
    print("\n", "Run: step8_sync_and_recalibrate_bam.sh -->", "\n")
    strScript = DIR2ndPipeline + "/step8_sync_and_recalibrate_bam.sh"
    CMD = "bash " + strScript + " " + strManifestFile + " " + "GERMLINE" + " " + strKTName
    os.system(CMD)
    return 1


# step 4
def ConstructBAMRecaliPerManifest(iSubjectNum, strKTName):
    # 1: Get Manifest file
    strManifestFile = DIR2ndPipelineTMP + "/" + strKTName + "/Manifest_Mimic.csv"
    if not os.path.exists(strManifestFile):
        print("Error: Manifest File is missing! ->", strManifestFile)
        return 1
    print("Got Manifest file:", strManifestFile)
    
    # 2: Create flag folder
    strFlagDir = DIRBAMRoot + "/" + strKTName + "/Flag/ConstructBAMRecaliPerManifest"
    strLogDir = DIRBAMRoot + "/" + strKTName + "/Log/ConstructBAMRecaliPerManifest"
    
    if not os.path.exists(strFlagDir):
        CMD = "mkdir -p " + strFlagDir
        os.system(CMD)
        
    if not os.path.exists(strLogDir):
        CMD = "mkdir -p " + strLogDir
        os.system(CMD)
        
    strFlagWorking = strFlagDir + "/construct_BAM_recali_per_manifest.flag.working"
    strFlagDone = strFlagDir + "/construct_BAM_recali_per_manifest.flag.done"
    
    if os.path.exists(strFlagWorking):
        print("ConstructBAMRecaliPerManifest is still working!")
        return 1
    
    if os.path.exists(strFlagDone):
        print("ConstructBAMRecaliPerManifest has been finished!")
        return 0

    # touch working flag 
    CMD = "touch " + strFlagWorking
    os.system(CMD)
    
    strBuildName = "build_CGR_" + strKTName
    strScriptBash = ("bash " + DIR2ndPipeline + "/step9_construct_BAM_recaliberated_per_manifest.sh" + " " + 
                                                strManifestFile + " " + 
                                                strBuildName + " " + 
                                                "GERMLINE" + " " + 
                                                strFlagWorking + " " + 
                                                strFlagDone)
    
    # Submit jobs ->
    strNumCore = "1"
    strNodeNum = "1"
    strJobName = "RBD-" + strKTName
    strRunningTime = "2-00:00:00"
    strStdOut = strLogDir + "/construct_BAM_recali_per_manifest.std.out"
    strStdErr = strLogDir + "/construct_BAM_recali_per_manifest.std.err"
    CMD = ("sbatch " + "--ntasks=" + strNumCore + " " + 
                       "--nodes=" + strNodeNum + " " + 
                       "--job-name=" + strJobName + " " + 
                       "--time=" + strRunningTime + " " + 
                       "--output=" + strStdOut + " " + 
                       "--error=" + strStdErr + " " + 
                       "--wrap=\"" + strScriptBash + "\"")
    print(CMD)
    os.system(CMD)
    return 1                                    


# step 5
def CoverageReport(iSubjectNum, strKTName):
    # 1: Get Manifest file
    strManifestFile = DIR2ndPipelineTMP + "/" + strKTName + "/Manifest_Mimic.csv"
    if not os.path.exists(strManifestFile):
        print("Error: Manifest File is missing! ->", strManifestFile)
        return 1
    print("Got Manifest file:", strManifestFile)
    
    # 2: Remaining jobs
    print("iSubjectNum      :", iSubjectNum)
    # Pre check if everything is all set! -->
    strFlagDir = DIRBAMRoot + "/" + strKTName + "/Flag/CoverageReport"
    
    # For step5_2 -> generate summarized report
    strFlagSumDir = DIRBAMRoot + "/" + strKTName + "/Flag/CoverageReportSummary"
    if not os.path.exists(strFlagSumDir):
        CMD = "mkdir -p " + strFlagSumDir
        os.system(CMD)
         
    if os.path.exists(strFlagDir):
        CMD = "find " + strFlagDir + " -iname '*.done' | wc -l"
        iDoneNum = int(subprocess.getoutput(CMD))
        print("iDoneNum     :", iDoneNum)
        
        CMD = "find " + strFlagDir + " -iname '*.flag.*' | wc -l"
        iTotalFlagNum = int(subprocess.getoutput(CMD))
        print("iTotalFlagNum:", iTotalFlagNum)
        
        if iDoneNum == iSubjectNum:
            print("CoverageReport is All Set!")
            strFlagDone = strFlagSumDir + "/coverage_report_sum.flag.done"
            if not os.path.exists(strFlagDone):
                # Collect Report
                strScript = DIR2ndPipeline + "/step5_2_generate_coverage_report_batch.sh " + strKTName
                CMD = "bash " + strScript 
                os.system(CMD)
                CMD = "touch " + strFlagDone
                os.system(CMD)
                print("Collect Coverage report is All Set!")
            else:
                print("Collect Coverage report has been done before!")
            print("Log file is located in            :", "/data/COVID_WGS/lix33/Test/2ndpipeline/Data/secondary_buf/cluster_job_logs/...")
            print("Coverage report file is located in:", "/data/COVID_WGS/lix33/Test/2ndpipeline/Data/secondary_buf/coverage_report/...")
            return 0
        
        if iTotalFlagNum == iSubjectNum:
            print("All jobs in the phase of CoverageReport are still running!")
            return 1
    
    # if not, go ahead to call 
    strBuildName = "build_CGR_" + strKTName
    strDirBuildBAM = DIRBuildProcess + "/" + strBuildName + "/bam_location/WGS" 
    
    print("\n", "Run: step5_generate_coverage_report_batch.sh -->", "\n")
    strScript = DIR2ndPipeline + "/step5_generate_coverage_report_batch.sh"
    CMD = "bash " + strScript + " " + strManifestFile + " " + strDirBuildBAM + " " + strKTName
    print(CMD)
    os.system(CMD)
    return 1


# Step 6
def PreCallingQCReport(iSubjectNum, strKTName):
    # 1: Get Manifest file
    strManifestFile = DIR2ndPipelineTMP + "/" + strKTName + "/Manifest_Mimic.csv"
    if not os.path.exists(strManifestFile):
        print("Error: Manifest File is missing! ->", strManifestFile)
        return 1
    print("Got Manifest file:", strManifestFile)
    
    # 2: Remaining jobs
    print("iSubjectNum      :", iSubjectNum)
    # Pre check if everything is all set! -->
    strFlagDir = DIRBAMRoot + "/" + strKTName + "/Flag/PreCallingQCReport"
    
    # For step6_2 -> generate summarized report
    strFlagSumDir = DIRBAMRoot + "/" + strKTName + "/Flag/PreCallingQCReportSummary"
    if not os.path.exists(strFlagSumDir):
        CMD = "mkdir -p " + strFlagSumDir
        os.system(CMD)
         
    if os.path.exists(strFlagDir):
        CMD = "find " + strFlagDir + " -iname '*.done' | wc -l"
        iDoneNum = int(subprocess.getoutput(CMD))
        print("iDoneNum     :", iDoneNum)
        
        CMD = "find " + strFlagDir + " -iname '*.flag.*' | wc -l"
        iTotalFlagNum = int(subprocess.getoutput(CMD))
        print("iTotalFlagNum:", iTotalFlagNum)
        
        if iDoneNum == iSubjectNum:
            print("CoverageReport is All Set!")
            strFlagDone = strFlagSumDir + "/pre_calling_qc_report_sum.flag.done"
            if not os.path.exists(strFlagDone):
                # Find pre calling qc report
                CMD = ("find " + DIRBuffer + "/coverage_report" + " -maxdepth 1 -type f -iname 'pre_calling_qc_report*" + strKTName + "*' " + 
                       "-printf '%T@ %Tc %p\n' | sort -r | head -n 1 | awk '{print $NF}'")
                print("CMD", CMD)
                strFileList = subprocess.getoutput(CMD)
                if strFileList == "":
                    print("Error: bad pre calling QC list:", strFileList)
                    return 1
                print("strFileList:", strFileList)
                strQCReport = strFileList.split('\n')[0] 
                                
                # Collect Report
                strScript = DIR2ndPipeline + "/step6_2_generating_pre_calling_qc_report_batch.sh"
                CMD = "bash " + strScript + " " + strQCReport + " " + strKTName
                # print("CMD:", CMD)
                os.system(CMD)
                CMD = "touch " + strFlagDone
                os.system(CMD)                
                print("Generating pre-calling QC report is All Set!")
                print("Pre-calling QC report is located in:", strQCReport)
            else:
                print("Generating pre-calling QC report has been done before!")
                print("Pre-calling QC report is located in:", "/data/COVID_WGS/lix33/Test/2ndpipeline/Data/secondary_buf/coverage_report/...")
            
            # print("Log file is located in            :", "/data/COVID_WGS/lix33/Test/2ndpipeline/Data/secondary_buf/cluster_job_logs/...")
            # print("Coverage report file is located in:", "/data/COVID_WGS/lix33/Test/2ndpipeline/Data/secondary_buf/coverage_report/...")
            return 0
        
        if iTotalFlagNum == iSubjectNum:
            print("All jobs in the phase of CoverageReport are still running!")
            return 1
    
    # if not, go ahead to call 
    strBuildName = "build_CGR_" + strKTName
    strDirBuildBAM = DIRBuildProcess + "/" + strBuildName + "/bam_location/WGS" 
    
    print("\n", "Run: step6_generate_pre_calling_qc_report_batch.sh -->", "\n")
    strScript = DIR2ndPipeline + "/step6_generate_pre_calling_qc_report_batch.sh"
    CMD = "bash " + strScript + " " + strManifestFile + " " + strDirBuildBAM + " " + strKTName
    print(CMD)
    os.system(CMD)
    return 1


def BAMContaminationCheck(iSubjectNum, strKTName):
    # Check if all build dir contains the BAM contamination report
    CMD = "find " + DIRBuildProcess + " -type f -iname '*SumContaminationReport.csv'"
    strReportList = subprocess.getoutput(CMD)
    
    CMD = "find " + DIRBuildProcess + " -maxdepth 1 -type d -iname 'build_*'"
    strBuildList = subprocess.getoutput(CMD)
        
    if strReportList != "" and strBuildList != "":
        vReportList = strReportList.split('\n')
        vBuildList = strBuildList.split('\n')
        if len(vReportList) == len(vBuildList):
            print("Bam Contamination Check for all Build has been all set!")
            print("vReportList:", vReportList)
            print("vBuildList :", vBuildList)            
            return 0
        
    # If not, run python script
    strScript = DIRCustomizedQC + "/Tools/BAMContaminationCheck/ContaminationCheck.py " + strKTName
    CMD = "python3 " + strScript
    os.system(CMD)
    return 1


def SendEmailNotification(iSubjectNum, strKTName):
    strCurBuildFolder = DIRBAMRoot + "/" + strKTName
    strFlag = strCurBuildFolder + "/Flag/" + DONEFlagEmail
    
    strReportDir = strCurBuildFolder + "/Report"
    # if True:
    if not os.path.exists(strFlag): 
        # Send email & add new email done flag
        strMsg = "============ " + strKTName + " ============ \\n\\n"
        strMsg += ("Good News              : The keytable " + strKTName + " IS ALL SET (COVID 2nd Analysis Pipeline) \\n\\n" + 
                   "Total Number of Subject: " + str(iSubjectNum) + "\\n\\n" + 
                   "\\n\\n" + 
                   "======== Please Check the Results Below (Biowulf) ========  \\n\\n" + 
                   "Re-calibrated BAM File Path        : " + DIRBAMReformat + "/BAM_recalibrated/WGS" + "\\n\\n" + 
                   "Coverage Report Path               : " + DIRBuffer + "/coverage_report/{contain keywords coverage}" + "\\n\\n" + 
                   "Pre-calling QC Report Path         : " + DIRBuffer + "/coverage_report/{contain keywords pre-calling QC}" + "\\n\\n" + 
                   "BAM Contamination Check Report Path: " + DIRBuildProcess + "/build_CGR_" + strKTName + "/Report" + "\\n\\n" + 
                   "\\n\\n" + 
                   "======== Some Useful Directories are also listed below (Biowulf) ========  \\n\\n" + 
                   "Build Dir                     : " + DIRBuildProcess + "/build_CGR_" + strKTName + "\\n\\n" + 
                   "Original BAM Dir (after Merge): " + DIRBAMReformat + "/ BAM_original/WGS" + "\\n\\n" + 
                   "Retrieved Flowcell (From S3)  : " + DIRBAMRoot + "/" + strKTName + "\\n\\n" + 
                   "Flag Dir                      : " + DIRBAMReformat + "/BAM_recalibrated/Flag" + "\\n\\n" + 
                   "      >>>>>>>>>>>>>>>>>>>>>>> : " + DIRBAMRoot + "/" + strKTName + "/Flag" + "\\n\\n" + 
                   "      >>>>>>>>>>>>>>>>>>>>>>> : " + DIRBuildProcess + "/build_CGR_" + strKTName + "/Flag" + "\\n\\n" + 
                   "Log Dir:                      : " + DIRBuffer + "/cluster_job_logs" + "\\n\\n" + 
                   "      >>>>>>>>>>>>>>>>>>>>>>> : " + DIR2ndPipelineTMP + "/" + strKTName + "/Log" + "\\n\\n" + 
                   "      >>>>>>>>>>>>>>>>>>>>>>> : " + DIRBAMRoot + "/" + strKTName + "/Flag" + "\\n\\n" + 
                   "      >>>>>>>>>>>>>>>>>>>>>>> : " + DIR2ndPipelineTMP + "/" + strKTName + "/Script" + "\\n\\n" + 
                   "      >>>>>>>>>>>>>>>>>>>>>>> : " + DIRBuildProcess + "/build_CGR_" + strKTName + "/Log" + "\\n\\n" + 
                   "Report Dir                    : " + DIRBuffer + "/coverage_report" + "\\n\\n" + 
                   "Pre-QC Tmp Result Dir         : " + DIRBuffer + "/PRE_QC" + "\\n\\n" + 
                   "BAM Contamination Check Dir   : " + DIRBuildProcess + "/build_CGR_" + strKTName + "/Report,Log,Flag" + "\\n\\n" + 
                   "Mimic Manifest Dir            : " + DIR2ndPipelineTMP + "/" + strKTName + "/Manifest_Mimic.csv") 
                                                                            
        strSubject = strKTName + " is all Set (COVID 2nd Analysis Pipeline)"
        
        # ---> Prepare Attachment
        CMD = "find " + strReportDir + " -maxdepth 1 -type f"
        vFileList = subprocess.getoutput(CMD).split('\n')
        print("vFileList:", vFileList)
        strAttach = ""
        for file in vFileList:
            strAttach += " -a " + file
        # <---
                        
        CMD = "echo -e \"" + strMsg + "\" | mail -r " + EMAILSender + strAttach + " -s \"" + strSubject + "\" " + EMAILReceiverEnd            
        print(CMD)               
        os.system(CMD)                 
        # 2: Set working flag to done and                
        CMD = "touch " + strFlag
        os.system(CMD)
        print("Email has been sent successfully! -->", strKTName)
    else:
        print("No action needed! Analysis has been finished! -->", strKTName)
    
            
def main(): 
    strKeytable = sys.argv[1]
    if not os.path.exists(strKeytable):
        print("Error: Keytable does not exist! -->", strKeytable)
        return
    
    # Print time stamp -->
    print()
    print("====== COVID Auto Framework ======", flush=True)
    print("strKeytable:", strKeytable)
    os.system("date")
    print()
    # <--
    
    # Get the total number of subject in current keytable -> Go after lunch ->
    CMD = "awk 'NR > 1' " + strKeytable + " | grep -v 'topoff' | wc -l"
    iSubjectNum = int(subprocess.getoutput(CMD))
    
    # Get Keytable Name
    strKTName = os.path.basename(strKeytable).split('.')[0]
    
    # --> Get subject info for current build (keytable)
    objBuild = ClsBuild()    
    # 0: Prepare group info 
    print("\n", "==> GetGroupInfo") 
    objBuild.GetGroupInfo(strKeytable)
    # <--
    
    # SendEmailNotification(iSubjectNum, strKTName)
         
    # Step 1: Merge BAM file
    print("\n-----\n", "Step 1: Merge BAM file", "\n-----\n")
    if MergeBAM(strKeytable, iSubjectNum, strKTName, objBuild) != 0:
        print("The phase of MergeBAM is still running!")
        return 1
    
    # Step 2:
    print("\n-----\n", "Step 2: Reconstruct BAM Dir", "\n-----\n")
    if ReconstructBAMDir(strKTName) != 0:
        print("The phase of ReconstructBAMDir is still running!")
        return 1
    
    # Step 3:
    print("\n-----\n", "Step 3: Recalibrate BAM file", "\n-----\n")
    if RecalibrateBAMFile(iSubjectNum, strKTName, objBuild) != 0:
        print("The phase of RecalibrateBAMFile is still running!")
        return 1
    
    # Step 4:
    print("\n-----\n", "Step 4: Construct BAM recaliberated per manifest", "\n-----\n")
    if ConstructBAMRecaliPerManifest(iSubjectNum, strKTName) != 0:
        print("The phase of ConstructBAMRecaliPerManifest is still running!")
        return 1
    
    # Step 5:     
    print("\n-----\n", "Step 5: Generate Coverage Report", "\n-----\n")
    if CoverageReport(iSubjectNum, strKTName) != 0:
        print("The phase of CoverageReport is still running!")
        return 1
    
    # Step 6:
    print("\n-----\n", "Step 6: generate pre-calling qc report", "\n-----\n")
    iReturnPreQC = PreCallingQCReport(iSubjectNum, strKTName)
    
    # Step 7:
    print("\n-----\n", "Step 7: BAM contamination check", "\n-----\n")
    iReturnBAMContamCheck = BAMContaminationCheck(iSubjectNum, strKTName)
    
    # Step 6 and step 7 can be run simutaniously
    if iReturnPreQC != 0:
        print("The phase of PreCallingQCReport is still running!")
        return 1
    if iReturnBAMContamCheck != 0:
        print("The phase of BAMContaminationCheck is still running!")
        return 1
    
    # Step 8: 
    SendEmailNotification(iSubjectNum, strKTName)
    
    print("Enjoy the day!")
    return 0
    

if __name__ == "__main__":
    main()
