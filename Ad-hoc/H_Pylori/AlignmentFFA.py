import sys
import os
import subprocess
import gzip

LSTRefName  = ["26695-1MET", "26695", "ELS37", "F57"]
SUFFIXRef   = "fsa"
SUFFIXReads = "ffa.gz"

DIRRef              = "/scratch/lix33/lxwg/Project/H_pylori/Data/to-Xin/ref"
DIROutputBAMRoot    = "/scratch/lix33/lxwg/Project/H_pylori/Processed/BAM"
DIROutputLogRoot    = "/scratch/lix33/lxwg/Project/H_pylori/Processed/Log"
DIROutputFlagRoot   = "/scratch/lix33/lxwg/Project/H_pylori/Processed/Flag"
DIRRawSample        = "/scratch/lix33/lxwg/Project/H_pylori/Data/to-Xin/1012set"
DIRUpdatedRawSample = "/scratch/lix33/lxwg/Project/H_pylori/Processed/UpdatedRawSample"

# DIRRawSample        = "/home/lixin/lxwg/Data/H-Pylori/RawReads"
# DIRUpdatedRawSample = "/home/lixin/lxwg/Data/H-Pylori/UpdatedRawSample"

CORE = "2"
BASHScript = "/scratch/lix33/lxwg/SourceCode/H_Pylori/Alignment.sh"

class ClsReads:
    def __init__(self):
        self.strName = ""
        self.strSeq = ""
        self.iIndex = 0
    
    def UpdateReadsName(self):
        if self.strName == "":
            return
        vItem = self.strName.split(' ')
        vItem[0] += "|" + str(self.iIndex)
        self.strName = ' '.join(vItem)

class ClsSample :
    def __init__(self):
        self.strName = ""
        self.strReads = ""
        self.strUpdatedReads = ""
    
    def Init(self, strReadsPath):
        self.strName = strReadsPath.split('/')[-1].split('.')[0]
        self.strReads = strReadsPath
    
    def UpdateReads(self):
        # 1: Get File Name
        strFileName = os.path.basename(self.strReads)
        # Update reads Name (by adding index in the end of the reads)
        vReads = []
        objReads = ClsReads()
        print(self.strReads)
        ifs = gzip.GzipFile(self.strReads)
        strSeq = ""
        iCount = 0
        for line in ifs: 
            strLine = line.decode("utf-8")
            if strLine[0] == '>':
                iCount += 1
                if strSeq != "":
                    objReads.strSeq = strSeq
                    vReads.append(objReads)
                    strSeq = ""
                objReads = ClsReads()
                objReads.strName = strLine
                objReads.iIndex = iCount
                objReads.UpdateReadsName()
            else:
                strSeq += strLine
        # Add the last reads
        objReads.strSeq = strSeq
        vReads.append(objReads) 
        ifs.close()
        
        # Save updated reads to file
        strUpdatedReads = DIRUpdatedRawSample + "/" + strFileName
        with gzip.open(strUpdatedReads, 'wb') as f:
            for reads in vReads:
                f.write(reads.strName.encode())
                f.write(reads.strSeq.encode())
        self.strUpdatedReads = strUpdatedReads

    def SubmitJob(self, strRefName, strRefFile, strOutputBAMDir, strOutputLogDir, strOutputFlagDir):
        
        flagWorking = strOutputFlagDir + "/" + self.strName + ".mapping.working"
        flagDone = strOutputFlagDir + "/" + self.strName + ".mapping.done"
        if os.path.exists(flagWorking) or os.path.exists(flagDone):
            print("Current Sample is working or done!")
            return
        
        strOutStd = strOutputLogDir + "/" + self.strName + ".out.std"
        strOutErr = strOutputLogDir + "/" + self.strName + ".out.err"
        strName = "RMP." + strRefName + "." + self.strName
        CMD = "qsub -cwd -q long.q -pe by_node " + CORE + " " + \
              "-o " + strOutStd + " -e " + strOutErr + " " + \
              "-N " + strName + " " + \
              "-S /bin/bash " + BASHScript + " " + strRefFile + " " + \
               self.strUpdatedReads + " " + self.strName + " " + \
               strOutputBAMDir + " " + strOutputFlagDir + " " + CORE
        print(CMD)
        os.system(CMD)

def InitSample(vSample):
    # Notice: ignore hidden files
    CMD = "find " + DIRRawSample + " -maxdepth 1 -mindepth 1 -not -path '*/.*' -type f -iname '*" + SUFFIXReads + "'"
    vReadsList = subprocess.getoutput(CMD).split('\n')
    #print(len(vReadsList), vReadsList[0])
    
    # Init Sample
    vSample.clear()
    for reads in vReadsList:
        objSample = ClsSample()
        objSample.Init(reads)
        vSample.append(objSample)
    
    #print(vSample[0].strName)
    

def MergeBAM(strRefName, strOutputBAMDir, strOutputFlagDir):
    # Get the number of done flag
    CMD = "find " + strOutputBAMDir + " -mindepth 1 -maxdepth 1 -type f -iname '*.bam'"
    vBAM = subprocess.getoutput(CMD).split('\n')
    
    CMD = "find " + strOutputFlagDir + " -mindepth 1 -maxdepth 1 -type f -iname '*.done'"
    vDoneFlag = subprocess.getoutput(CMD).split('\n')
    
    if len(vBAM) != len(vDoneFlag):
        print("vBAM     :", len(vBAM))
        print("vDoneFlag:", len(vDoneFlag))
        print("Files are not equal! In correct!")
        return
    
    print("vBAM     :", len(vBAM))
    print("vDoneFlag:", len(vDoneFlag))

    # merge BAM
    strBAMList = " ".join(vBAM)
    strMergedBAM = strOutputBAMDir + "/merged.all." + strRefName + ".bam"
    CMD = "samtools merge -@ 4 -f " + strMergedBAM + " " + strBAMList
    os.system(CMD)
    
    strMergedSortedBAM = strOutputBAMDir + "/merged.all." + strRefName + ".sorted.bam"
    CMD = "samtools sort -@ 4 -o " + strMergedSortedBAM + " " + strMergedBAM
    os.system(CMD)
    
    SortedIndex = strMergedSortedBAM + ".bai"
    CMD="samtools index -@ 4 " + strMergedSortedBAM + " " + SortedIndex
    os.system(CMD)
    
    CMD = "rm " + strMergedBAM
    os.system(CMD)
    
def main():
    # Prepare Sample
    vSample = []
    InitSample(vSample)
    
    for sample in vSample:
        sample.UpdateReads()
    
    # Prepare
    for strRefName in LSTRefName:
        # Ref File 
        strRefFile = DIRRef + "/" + strRefName + "." + SUFFIXRef
        # Output BAM Dir 
        strOutputBAMDir = DIROutputBAMRoot + "/" + strRefName
        CMD = "mkdir -p " + strOutputBAMDir
        os.system(CMD) 
        # Output Log Dir
        strOutputLogDir = DIROutputLogRoot + "/" + strRefName
        CMD = "mkdir -p " + strOutputLogDir
        os.system(CMD)
        # Output Log Dir
        strOutputFlagDir = DIROutputFlagRoot + "/" + strRefName
        CMD = "mkdir -p " + strOutputFlagDir
        os.system(CMD)
        
        for sample in vSample:
            sample.SubmitJob(strRefName, strRefFile, strOutputBAMDir, strOutputLogDir, strOutputFlagDir)
    
        # MergeBAM
        #MergeBAM(strRefName, strOutputBAMDir, strOutputFlagDir)
            
if __name__ == "__main__":
    main()