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


class ClsPhase:
    def __init__(self):
        self.strCmd = ""

def MergeBAM(strKeytable, iSubjectNum, strKTName):    
    print("iSubjectNum      :", iSubjectNum)
    # Pre check if everything is all set! -->
    strFlagDir = DIRBAMRoot + "/" + strKTName + "/Flag/MergeSubject"    
    if os.path.exists(strFlagDir):
        CMD = "find " + strFlagDir + " -iname '*.done' | wc -l"
        iDoneNum = int(subprocess.getoutput(CMD))
        print("iDoneNum     :", iDoneNum)
        
        CMD = "find " + strFlagDir + " -iname '*.flag.*' | wc -l"
        iTotalFlagNum = int(subprocess.getoutput(CMD))
        print("iTotalFlagNum:", iTotalFlagNum)
        
        if iDoneNum == iSubjectNum:
            print("MergeBAM is All Set!")
            return 0
        
        if iTotalFlagNum == iSubjectNum:
            print("All jobs in the phase of merge BAM are still running!")
            return 1
    
    # if not, go ahead to call 
    print("\n", "Run: MergeSubject.py -->", "\n")
    strScript = DIRCustomizedQC + "/Tools/MergeSubject/MergeSubject.py"
    CMD = "python " + strScript + " " + strKeytable
    os.system(CMD)
    return 1
    #<-- 

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
    
    #Submit jobs ->
    strNumCore = "1"
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
def RecalibrateBAMFile(iSubjectNum, strKTName):
    print("iSubjectNum      :", iSubjectNum)
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
    
    # if not, go ahead to call 
    print("\n", "Run: step8_sync_and_recalibrate_bam.sh -->", "\n")
    strScript = DIR2ndPipeline + "/step8_sync_and_recalibrate_bam.sh"
    CMD = "bash " + strScript + " " + strKTName
    os.system(CMD)
    return 1

# step 4
def ConstructBAMRecaliPerManifest(iSubjectNum, strKTName):
    #1: Get Manifest file
    strManifestFile = DIR2ndPipelineTMP + "/" + strKTName + "/Manifest_Mimic.csv"
    if not os.path.exists(strManifestFile):
        print("Error: Manifest File is missing! ->", strManifestFile)
        return 1
    print("Got Manifest file:", strManifestFile)
    
    #2: Create flag folder
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
    
    #Submit jobs ->
    strNumCore = "1"
    strNodeNum = "1"
    strJobName = "RBD-" + strKTName
    strRunningTime = "10-00:00:00"
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
                strScript = DIR2ndPipeline + "/step5_2_generate_coverage_report_batch.sh"
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
                CMD = "find " + DIRBuffer + "/coverage_report" + " -iname 'pre_calling_qc_report*'"
                strFileList = subprocess.getoutput(CMD)
                if strFileList == "":
                    print("Error: bad pre calling QC list:", strFileList)
                    return 1
                print("strFileList:", strFileList)
                strQCReport = strFileList.split('\n')[0] 
                                
                # Collect Report
                strScript = DIR2ndPipeline + "/step6_2_generating_pre_calling_qc_report_batch.sh"
                CMD = "bash " + strScript + " " + strQCReport 
                os.system(CMD)
                CMD = "touch " + strFlagDone
                os.system(CMD)
                print("Generating pre-calling QC report is All Set!")
            else:
                print("Generating pre-calling QC report has been done before!")
            
            print("Pre-calling QC report is located in:", "/data/COVID_WGS/lix33/Test/2ndpipeline/Data/secondary_buf/coverage_report/...")    
            #print("Log file is located in            :", "/data/COVID_WGS/lix33/Test/2ndpipeline/Data/secondary_buf/cluster_job_logs/...")
            #print("Coverage report file is located in:", "/data/COVID_WGS/lix33/Test/2ndpipeline/Data/secondary_buf/coverage_report/...")
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
        
    #If not, run python script
    strScript = DIRCustomizedQC + "/Tools/BAMContaminationCheck/ContaminationCheck.py"
    CMD = "python3 " + strScript
    os.system(CMD)
    return 1
            
def main():        
    strKeytable = sys.argv[1]
    if not os.path.exists(strKeytable):
        print("Error: Keytable does not exist! -->", strKeytable)
        return
    
    #Print time stamp -->
    print()
    print("====== COVID Auto Framework ======", flush=True)
    print("strKeytable:", strKeytable)
    os.system("date")
    print()
    #<--
    
    # Get the total number of subject in current keytable -> Go after lunch ->
    CMD = "awk 'NR > 1' " + strKeytable + " | grep -v 'topoff' | wc -l"
    iSubjectNum = int(subprocess.getoutput(CMD))
    
    #Get Keytable Name
    strKTName = os.path.basename(strKeytable).split('.')[0]    
         
    # Step 1: Merge BAM file
    print("\n-----\n", "Step 1: Merge BAM file", "\n-----\n")
    if MergeBAM(strKeytable, iSubjectNum, strKTName) != 0:
        print("The phase of MergeBAM is still running!")
        return 1
    
    # Step 2:
    print("\n-----\n", "Step 2: Reconstruct BAM Dir", "\n-----\n")
    if ReconstructBAMDir(strKTName) != 0:
        print("The phase of ReconstructBAMDir is still running!")
        return 1
    
    # Step 3:
    print("\n-----\n", "Step 3: Recalibrate BAM file", "\n-----\n")
    if RecalibrateBAMFile(iSubjectNum, strKTName) != 0:
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
    if PreCallingQCReport(iSubjectNum, strKTName) != 0:
        print("The phase of PreCallingQCReport is still running!")
        return 1
    
    # Step 7:
    print("\n-----\n", "Step 7: BAM contamination check", "\n-----\n")
    if BAMContaminationCheck(iSubjectNum, strKTName) != 0:
        print("The phase of BAMContaminationCheck is still running!")
        return 1
    
    print("Enjoy the day!")
    return 0
    

if __name__ == "__main__"