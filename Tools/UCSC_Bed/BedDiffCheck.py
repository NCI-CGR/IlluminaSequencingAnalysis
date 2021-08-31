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
        self.bHit = False
    
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

class ClsCmpUnit:
    def __init__(self):
        self.strChrom = ""
        self.iBaseNum1st = 0
        self.iBaseNum2nd = 0
        self.iBaseNumOverlap = 0
    
    def UnitCompare(self, chromRegion1st, chromRegion2nd):
        chromRegion1st.bHit = True
        chromRegion2nd.bHit = True
        self.strChrom = chromRegion1st.strChrom
        
        # Get iBaseNum1st
        self.iBaseNum1st = 0
        for region in chromRegion1st.vRegion:
            self.iBaseNum1st += region.iEnd - region.iStart
        
        # Get iBaseNum2nd
        self.iBaseNum2nd = 0
        for region in chromRegion2nd.vRegion:            
            self.iBaseNum2nd += region.iEnd - region.iStart
            
        #Get overlap base
        i = 0
        iIndex2nd = 0
        iOverlapBase = 0
        while i < len(chromRegion1st.vRegion):
            if iIndex2nd >= len(chromRegion2nd.vRegion):
                break
            
            if chromRegion1st.vRegion[i].iEnd <= chromRegion2nd.vRegion[iIndex2nd].iStart:
                #print(chromRegion1st.vRegion[i].iStart, chromRegion1st.vRegion[i].iEnd, "-->", chromRegion2nd.vRegion[iIndex2nd].iStart, chromRegion2nd.vRegion[iIndex2nd].iEnd)
                i += 1                
                continue
            # print("i:", i)
            # print(chromRegion1st.vRegion[i].iStart, chromRegion1st.vRegion[i].iEnd, "-->", chromRegion2nd.vRegion[iIndex2nd].iStart, chromRegion2nd.vRegion[iIndex2nd].iEnd)
            # print
            # if i >= 100:
            #     print("will break:", i)
            #     break;
            while iIndex2nd < len(chromRegion2nd.vRegion):
                if chromRegion1st.vRegion[i].iEnd > chromRegion2nd.vRegion[iIndex2nd].iStart and chromRegion1st.vRegion[i].iStart < chromRegion2nd.vRegion[iIndex2nd].iEnd:
                    #Get overlap
                    iStart = chromRegion1st.vRegion[i].iStart if chromRegion1st.vRegion[i].iStart > chromRegion2nd.vRegion[iIndex2nd].iStart else chromRegion2nd.vRegion[iIndex2nd].iStart
                    iEnd = chromRegion1st.vRegion[i].iEnd if chromRegion1st.vRegion[i].iEnd < chromRegion2nd.vRegion[iIndex2nd].iEnd else chromRegion2nd.vRegion[iIndex2nd].iEnd
                    # print(i, iIndex2nd)
                    # print(iStart, iEnd)
                    # print("1st:", chromRegion1st.vRegion[i].iStart, chromRegion1st.vRegion[i].iEnd)
                    # print("2nd:", chromRegion2nd.vRegion[iIndex2nd].iStart, chromRegion2nd.vRegion[iIndex2nd].iEnd)                    
                    iOverlapBase += iEnd - iStart
                    
                if chromRegion1st.vRegion[i].iEnd >= chromRegion2nd.vRegion[iIndex2nd].iEnd or chromRegion1st.vRegion[i].iStart > chromRegion2nd.vRegion[iIndex2nd].iEnd:
                    iIndex2nd += 1
                else:
                    break
            # print("New:", i, iIndex2nd)
            # print()
            i += 1            
        
        self.iBaseNumOverlap = iOverlapBase
        #print("iOverlapBase(" + self.strChrom + "):", iOverlapBase)
        
    def PrintStat(self):
        print(self.strChrom + " (OverlapBase):", self.iBaseNumOverlap)
        print("\t", "hg19 Base Num:", self.iBaseNum1st, " --> Overlap/Base_hg19:",  "{:.1%}".format(self.iBaseNumOverlap/self.iBaseNum1st) if self.iBaseNum1st !=0 else 0)
        print("\t", "v38 Base Num :", self.iBaseNum2nd, " --> Overlap/Base_v38 :",  "{:.1%}".format(self.iBaseNumOverlap/self.iBaseNum2nd) if self.iBaseNum2nd !=0 else 0)
        print("---")
        

class ClsBedComparison:
    def __init__(self):
        self.vCmpUnit = []
    
    def PrintStat(self):
        iTotalBaseOverlap = 0
        iTotalBase1st = 0
        iTotalBase2nd = 0
            
        for unit in self.vCmpUnit:
            iTotalBaseOverlap += unit.iBaseNumOverlap
            iTotalBase1st += unit.iBaseNum1st
            iTotalBase2nd += unit.iBaseNum2nd
        
        print("************** Summary *****************")
        print("Total Overlap Base:", iTotalBaseOverlap)
        print("hg19 Base Num     :", iTotalBase1st, " --> Overlap/Total_hg19:",  "{:.1%}".format(iTotalBaseOverlap/iTotalBase1st))
        print("v38 Base Num      :", iTotalBase2nd, " --> Overlap/Total_v38 :",  "{:.1%}".format(iTotalBaseOverlap/iTotalBase2nd))
        print("************** ******* *****************")
        
        print()
        print("======== Details ========")
        
        for unit in self.vCmpUnit:
            unit.PrintStat()

def BedComparison(vChromRegion1st, vChromRegion2nd):
    objBedComparison = ClsBedComparison()
     
    for chromRegion2nd in vChromRegion2nd:
        for chromRegion1st in vChromRegion1st:        
            if chromRegion1st.strChrom == chromRegion2nd.strChrom:
                #Sort it first
                chromRegion1st.vRegion.sort(key=operator.attrgetter('iStart'))
                chromRegion2nd.vRegion.sort(key=operator.attrgetter('iStart'))
                
                objCmpUnit = ClsCmpUnit()
                objCmpUnit.UnitCompare(chromRegion1st, chromRegion2nd)
                objBedComparison.vCmpUnit.append(objCmpUnit)
                break; 
        #break; 
    
    for chromRegion1st in vChromRegion1st:
        if not chromRegion1st.bHit:
            objCmpUnit = ClsCmpUnit()
            objCmpUnit.strChrom = chromRegion1st.strChrom
            for region in chromRegion1st.vRegion:
                objCmpUnit.iBaseNum1st += region.iEnd - region.iStart
            objBedComparison.vCmpUnit.append(objCmpUnit) 
            print(chromRegion1st.strChrom)
        
    for chromRegion2nd in vChromRegion2nd:
        if not chromRegion2nd.bHit:
            objCmpUnit = ClsCmpUnit()
            objCmpUnit.strChrom = chromRegion2nd.strChrom
            for region in chromRegion2nd.vRegion:
                objCmpUnit.iBaseNum2nd += region.iEnd - region.iStart
            objBedComparison.vCmpUnit.append(objCmpUnit) 
            print(chromRegion2nd.strChrom)
    
    objBedComparison.PrintStat()
         
    print("Good to know!", len(objBedComparison.vCmpUnit))

def ParseBed(vChromRegion, strBedFile):
    file = open(strBedFile,'r')
    iNum = 100
    iIndex = 0
    iBaseNum = 0
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
                iBaseNum += objRegion.iEnd - objRegion.iStart
                bFind = True
                break
        if not bFind:
            objChromRegion = ClsChromRegion()
            objChromRegion.strChrom = strChrom
            
            objRegion = ClsRegion()
            objRegion.iStart = int(vItem[1])
            objRegion.iEnd = int(vItem[2])
            objRegion.strComment = "\t".join(vItem[3:])
            iBaseNum += objRegion.iEnd - objRegion.iStart  
            
            objChromRegion.vRegion.append(objRegion)
            vChromRegion.append(objChromRegion)
        # iIndex += 1
        # if iIndex > iNum:
        #     break    
    file.close()
    # print("BaseNum :-->", iBaseNum)
    # for chromRegion in vChromRegion:
    #     chromRegion.Print()
    # print()
    

def main():
    strBedFile = sys.argv[1]
    vChromRegion1st = []
    ParseBed(vChromRegion1st, strBedFile)
    # for chromregion in vChromRegion1st:
    #     print(chromregion.strChrom)
    
    strBedFile = sys.argv[2]
    vChromRegion2nd = []
    ParseBed(vChromRegion2nd, strBedFile)
    
    BedComparison(vChromRegion1st, vChromRegion2nd)
    
    print("Everything is all set!")
    
    
if __name__ == "__main__":        
    main()