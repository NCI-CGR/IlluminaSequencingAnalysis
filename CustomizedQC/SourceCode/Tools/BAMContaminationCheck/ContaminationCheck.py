'''
1: Need to have a flag dir
'''

import os
import sys
import subprocess

DIRRootBuild = "/data/COVID_WGS/lix33/Test/2ndpipeline/Build/processed"
VCFRef = "/data/COVID_WGS/lix33/DCEG/CGF/Bioinformatics/Production/data/refVariant/Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.sites_with_chr.vcf"

class ClsSubject:
    def __init__(self):
        self.strBAM = ""
        self.strFlagDir = ""
        self.strReportDir = ""
        self.strLogDir = ""
    
    def Init(self, strBAMFile, strFlagDir, strReportDir, strLogDir):        
        self.strBAM = strBAMFile
        self.strFlagDir = strFlagDir
        self.strReportDir = strReportDir
        self.strLogDir = strLogDir
    
    def Run(self):
        strBAMName = os.path.basename(self.strBAM).split('.')[0]        
        strFlagWorking = self.strFlagDir + "/" + strBAMName + ".flag.working"  
        strFlagDone = self.strFlagDir + "/" + strBAMName + ".flag.done"
        
        if os.path.exists(strFlagDone):
            return 1
        
        if os.path.exists(strFlagWorking):
            return 0
        
        # Nothing -> try to submit jobs
        CMD = "touch " + strFlagWorking
        os.system(CMD)        
        # submit jobs    
        dirCurFile = os.path.dirname(__file__)
        print(dirCurFile)
        scriptRun = dirCurFile + "/ContaminationCheckSingle.sh"
        #strVCF = ""
        outputPrefix = os.path.basename(self.strBAM).split('.')[0]
        strReportFile = self.strReportDir + "/" + outputPrefix
        
        scriptBash = ("bash " + scriptRun + " \"" + VCFRef + "\"" + 
                               " \"" + self.strBAM + "\"" +
                               " \"" + strReportFile + "\"" +
                               " \"" + strFlagWorking + "\"" +
                               " \"" + strFlagDone + "\"")        
        #Submit jobs ->
        strNumCore = "1"
        strNodeNum = "1"
        strJobName = "BCC-" + outputPrefix
        strRunningTime = "10-00:00:00"
        strStdOut = self.strLogDir + "/" + outputPrefix + ".std.out"
        strStdErr = self.strLogDir + "/" + outputPrefix + ".std.err"
        CMD = ("sbatch " + "--ntasks=" + strNumCore + " " +  
                           "--nodes=" + strNodeNum + " " +
                           "--job-name=" + strJobName + " " + 
                           "--time=" + strRunningTime + " " + 
                           "--output=" + strStdOut + " " +
                           "--error=" + strStdErr + " " + 
                           "--wrap=\"" + scriptBash + "\"")
        print(CMD)
        os.system(CMD)
        #<- 

class ClsBuild:
    def __init__(self):
        self.vSubject = []
        self.strBuildName = ""
        self.strRootDir = ""
        self.strFlagDir = ""
        self.strReportDir = ""
        self.strLogDir = ""
    
    def InitByDir(self, strBuildDir):        
        self.strBuildName = os.path.basename(strBuildDir)
        
        self.strRootDir = strBuildDir
        # Notice: do not use "-type f ", since the "bam" file may be a softlink
        CMD = "find " + strBuildDir + " -iname '*.bam'" 
        #print(CMD)
        strBAMList = subprocess.getoutput(CMD)
        #print("strBAMList:", strBAMList)
        if strBAMList == "":
            return 
        
        #Init Flag Dir      
        self.strFlagDir = strBuildDir + "/Flag/ContaminationCheck"
        if not os.path.exists(self.strFlagDir):
            CMD = "mkdir -p " + self.strFlagDir
            os.system(CMD)
        
        #Init Report Dir
        self.strReportDir = strBuildDir + "/Report/ContaminationCheck"
        if not os.path.exists(self.strReportDir):
            CMD = "mkdir -p " + self.strReportDir
            os.system(CMD)
        
        #Init Log Dir
        self.strLogDir = strBuildDir + "/Log/ContaminationCheck"
        if not os.path.exists(self.strLogDir):
            CMD = "mkdir -p " + self.strLogDir
            os.system(CMD)
        
        # Init subject
        vBAMFileList = strBAMList.split('\n')
        self.vSubject.clear()
        for strBAMFile in vBAMFileList:
            objSubject = ClsSubject()
            objSubject.Init(strBAMFile, self.strFlagDir, self.strReportDir, self.strLogDir)
            self.vSubject.append(objSubject)
    
    def Run(self):
        for subject in self.vSubject:
            subject.Run()
    
    def BuildReport(self):
        #Check the number of "done" flag
        CMD = "find " + self.strFlagDir + " -iname '*done' | wc -l"
        iDoneFlagNum = int(subprocess.getoutput(CMD))        
         
        #Check the number of "selfSM" file
        CMD = "find " + self.strReportDir + " -iname '*selfSM' | wc -l"
        iReportNum = int(subprocess.getoutput(CMD))
        
        if iDoneFlagNum == 0 or iReportNum == 0 or iDoneFlagNum != iReportNum:
            print("self.strFlagDir  :", self.strFlagDir)
            print("self.strReportDir:", self.strReportDir)
            print("iDoneFlagNum     :", iDoneFlagNum)
            print("iReportNum       :", iReportNum)
            print("Waiting for the run be finished!")
            return
        
        CMD = "find " + self.strReportDir + " -iname '*selfSM'"
        vReport = subprocess.getoutput(CMD).split('\n')
        
        #Build VCF File to show the bam contamination checking result
        strSumReport = self.strReportDir + "/" + self.strBuildName + "_SumContaminationReport.csv"
        if os.path.exists(strSumReport):
            print("Report was existing: ", strSumReport)
            print("No further action needed!")
            return
        
        # Add header in to the CSV File
        fs = open(strSumReport, 'w')
        strHead = "SubjectName,FREEMIX\n"
        fs.write(strHead)
        for reportFile in vReport:
            strSubjectName = os.path.basename(reportFile).split('.')[0]
            CMD = "awk '{print $7}' " + reportFile + " | tail -n 1"  
            strCCRatio = subprocess.getoutput(CMD)
            strLine = strSubjectName + "," + strCCRatio + "\n"
            fs.write(strLine)
        fs.close()        
        print("Report location:", strSumReport)
        print("Report All Set!")

def GetBuildInfo(vBuild):
    vBuild.clear()
    CMD = "find " + DIRRootBuild + " -mindepth 1 -maxdepth 1 -type d"
    strBuildList = subprocess.getoutput(CMD)
    if strBuildList == "":
        return 
    vBuildDirList = strBuildList.split('\n')
    #print(vBuildDirList)
    for strBuildDir in vBuildDirList:
        objBuild = ClsBuild()
        objBuild.InitByDir(strBuildDir)
        vBuild.append(objBuild)
        
def RunContaminationCheck(vBuild):
    for build in vBuild:
        build.Run()

def BuildReport(vBuild):
    for build in vBuild:
        build.BuildReport()    

def main():    
    vBuild = []
    
    # Check each sub-sub folder
    GetBuildInfo(vBuild)
    
    # Submit jobs
    RunContaminationCheck(vBuild)
    
    # Build BAM Contamination Report
    BuildReport(vBuild)  
    
    print("Main Finished!")
    
if __name__ == "__main__":
    main()