import sys
import os
import operator

class ClsRegion:
    def __init__(self):
        #self.strChrom = ""
        self.iStart = ""
        self.iEnd = ""
        self.strComment = ""
    
    def Write2File(self, strChrom, f):
        strLine = strChrom + "\t" + str(self.iStart) + "\t" + str(self.iEnd) + "\t" + self.strComment
        f.write(strLine)
    
    def Write2FileBrief(self, strChrom, f):
        strLine = strChrom + "\t" + str(self.iStart) + "\t" + str(self.iEnd) + "\n"
        f.write(strLine)
    
    def Print(self):
        print(self.iStart, self.iEnd)
    
class ClsChromRegion:
    def __init__(self):
        self.strChrom = ""
        self.vRegion = []
    
    def GetMergedRegion(self, vOrgRegion):
        self.vRegion.clear()
        i = 0
        while i < len(vOrgRegion):
            if i + 1 < len(vOrgRegion):
                iStart = vOrgRegion[i].iStart
                iEnd = vOrgRegion[i].iEnd
                strComment = vOrgRegion[i].strComment
                while i + 1 < len(vOrgRegion) and iEnd >= vOrgRegion[i+1].iStart:
                    #print(vOrgRegion[i].iEnd, vOrgRegion[i+1].iEnd)
                    if iEnd < vOrgRegion[i+1].iEnd:
                        iEnd = vOrgRegion[i+1].iEnd
                    i += 1
                    
                objRegion = ClsRegion()
                objRegion.iStart = iStart
                objRegion.iEnd = iEnd
                objRegion.strComment = strComment
                self.vRegion.append(objRegion)
                i += 1
                continue
                    
            if i == len(vOrgRegion) - 1:
                objRegion = ClsRegion()
                objRegion.iStart = vOrgRegion[i].iStart
                objRegion.iEnd = vOrgRegion[i].iEnd
                objRegion.strComment = vOrgRegion[i].strComment
                self.vRegion.append(objRegion)
                i += 1
                continue
    
    def Write2File(self, f):
        for region in self.vRegion:
            region.Write2File(self.strChrom, f)
    
    def Write2FileBrief(self, f):
        for region in self.vRegion:
            region.Write2FileBrief(self.strChrom, f)
    
    def Print(self):
        print(self.strChrom)
        for objRegion in self.vRegion:
            objRegion.Print()

def ParseBed(vChromRegion, strBedFile):
    file = open(strBedFile,'r')
    iNum = 100
    iIndex = 0
    while True:
        strLine = file.readline()    
        if not strLine:
            break;                        
        vItem = strLine.split('\t')
        strChrom = vItem[0]
        bFind = False        
        for chromRegion in vChromRegion:
            if chromRegion.strChrom == strChrom:
                objRegion = ClsRegion()
                objRegion.iStart = int(vItem[1])
                objRegion.iEnd = int(vItem[2])
                if len(vItem) > 3:
                    objRegion.strComment = "\t".join(vItem[3:])
                chromRegion.vRegion.append(objRegion)
                bFind = True
                break
        if not bFind:
            objChromRegion = ClsChromRegion()
            objChromRegion.strChrom = strChrom
            
            objRegion = ClsRegion()
            objRegion.iStart = int(vItem[1])
            objRegion.iEnd = int(vItem[2])
            objRegion.strComment = "\t".join(vItem[3:])
            
            objChromRegion.vRegion.append(objRegion)
            vChromRegion.append(objChromRegion)
        # iIndex += 1
        # if iIndex > iNum:
        #     break    
    file.close()
    
    # for chromRegion in vChromRegion:
    #     chromRegion.Print()
    # print()

def MergeOverlap(vChromRegion, vMergedChromRegion):
    for chromRegion in vChromRegion:
        objChromRegion = ClsChromRegion()
        objChromRegion.strChrom = chromRegion.strChrom
        # Sort 
        chromRegion.vRegion.sort(key=operator.attrgetter('iStart'))
        # Merge
        objChromRegion.GetMergedRegion(chromRegion.vRegion)
        vMergedChromRegion.append(objChromRegion)
    
    # for chromRegion in vChromRegion:
    #     chromRegion.Print()
    # print()
    #
    # for chromRegion in vMergedChromRegion:
    #     chromRegion.Print()
    # print()

def SaveUpdatedBAM(vMergedChromRegion, strOrgBedFile):
    # Regular: full info
    strFileName = os.path.basename(strOrgBedFile).split(".")[0] + ".MergedOverlap.bed"
    strFilePath = os.path.dirname(strOrgBedFile) + "/" + strFileName
    if os.path.exists(strFilePath):
        CMD = "rm " + strFilePath
        os.system(CMD)
        
    f = open(strFilePath, "a")
    for chromRegion in vMergedChromRegion:
        chromRegion.Write2File(f)
    f.close()
    
    # Brief version: the first 3 columns
    strFileName = os.path.basename(strOrgBedFile).split(".")[0] + ".MergedOverlap.Brief.bed"
    strFilePath = os.path.dirname(strOrgBedFile) + "/" + strFileName
    if os.path.exists(strFilePath):
        CMD = "rm " + strFilePath
        os.system(CMD)
        
    f = open(strFilePath, "a")
    for chromRegion in vMergedChromRegion:
        chromRegion.Write2FileBrief(f)
    f.close()
    
    
    print("All Set!")

def main():    
    strBedFile = sys.argv[1]
    vChromRegion = []
    ParseBed(vChromRegion, strBedFile)
    
    vMergedChromRegion = []
    MergeOverlap(vChromRegion, vMergedChromRegion)
    
    SaveUpdatedBAM(vMergedChromRegion, strBedFile)

if __name__ == "__main__":        
    main()