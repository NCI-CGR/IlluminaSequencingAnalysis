import pysam
import csv

DICCIGAR = {0: "M",
            1: "I",
            2: "D",
            3: "N",
            4: "S",
            5: "H",
            6: "P",
            7: "=",
            8: "X",
            9: "B"}

DICCIGAREDRTAIL = {0: "BAM_CMATCH",
                   1: "BAM_CINS",
                   2: "BAM_CDEL",
                   3: "BAM_CREF_SKIP",
                   4: "BAM_CSOFT_CLIP",
                   5: "BAM_CHARD_CLIP",
                   6: "BAM_CPAD",
                   7: "BAM_CEQUAL",
                   8: "BAM_CDIFF",
                   9: "BAM_CBACK"}


class ClsVariant():
    def __init__(self):
        self.strReadsName = ""
        self.iPos = 0
        self.strType = ""
        self.iLen = 0
    
    def Init(self, strName, iPos, tuple):
        self.strReadsName = strName
        self.iPos = iPos
        self.strType = DICCIGAREDRTAIL[tuple[0]]
        self.iLen = int(tuple[1])
    
    def GetInfo(self, listRow):
        listRow.clear()
        listRow.append(self.strReadsName)
        listRow.append(self.iPos)
        listRow.append(self.strType)
        listRow.append(self.iLen)

class ClsReads:
    def __init__(self):
        self.strName = ""
        self.iPos = ""
        self.strCIGAR = ""
        self.listCIGAR = []
        self.iAlignedLen = ""
        self.iReadsLength = ""
        self.strReadsSeq = ""
            
    def Init(self, reads):
        self.strName = reads.query_name
        self.iPos = reads.reference_start
        self.strCIGAR = reads.cigarstring
        self.listCIGAR = reads.cigartuples
        self.iReadsLength = reads.query_length
        self.iAlignedLen = reads.reference_length
        self.strReadsSeq = reads.get_forward_sequence()
            
    def GetVariant(self, vVariant, vNoneReads):
        if type(self.listCIGAR) == type(None):
            vNoneReads.append(self)
            return
        iPos = int(self.iPos)
        for tuple in self.listCIGAR:
            #skip the cases of "M" and "=" 
            if DICCIGAR[tuple[0]] == "M" or DICCIGAR[tuple[0]] == "=":
                iPos += int(tuple[1])
                continue
            objVariant = ClsVariant()
            objVariant.Init(self.strName, iPos, tuple)
            vVariant.append(objVariant)
            if DICCIGAR[tuple[0]] == "D" or DICCIGAR[tuple[0]] == "N" or DICCIGAR[tuple[0]] == "X": 
                iPos += int(tuple[1])

def OutputVairant(vVariant, strVariantFile):
    #Sort vVariant
    vVariant.sort(key=lambda x: x.iPos, reverse=False)
    
    f = open(strVariantFile, 'w')
    writer = csv.writer(f)
    listRow = ["GeneName", "Position", "AlignmentType", "Length"]
    writer.writerow(listRow)
    for variant in vVariant:
        variant.GetInfo(listRow)
        writer.writerow(listRow)
    f.close()
    

def OutputNoneReads(vNoneReads, strUnalignedReads):
    f = open(strUnalignedReads, 'w')
    for reads in vNoneReads:
        f.write(">" + reads.strName + "\n")
        f.write(reads.strReadsSeq + "\n")
    f.close()

def AlignAnalyze(dicAlignment):
    # for key in dicAlignment.keys():
    #     if len(dicAlignment[key]) > 1:
    #         print(key)
    #         print(len(dicAlignment[key]))
    #         for read in dicAlignment[key]:
    #             print(read.iAlignedLen, read.strCIGAR) 
    #         print("----------")
    
    # Get best Alignment (refine the dictionary)
    for key in dicAlignment.keys():
        if len(dicAlignment[key]) > 1:
            objBestReads = dicAlignment[key][0]
            for read in dicAlignment[key]:
                if read.iAlignedLen > objBestReads.iAlignedLen:
                    objBestReads = read
            dicAlignment[key].clear() 
            dicAlignment[key].append(objBestReads)
            #print(dicAlignment[key][0].iAlignedLen)
            #print("----------")
    
    # Next Export the diff parts Go!!
    vVariant = []
    vNoneReads = []
    for key in dicAlignment.keys():
        for reads in dicAlignment[key]:
            reads.GetVariant(vVariant, vNoneReads)
    
    strVariantFile = "/home/lixin/lxwg/ad-hoc/H_Pylori/CallVariant/VariantProfile.vcf"
    OutputVairant(vVariant, strVariantFile)
    
    strUnalignedGenome = "/home/lixin/lxwg/ad-hoc/H_Pylori/CallVariant/UnalignedGenome.ffa"
    OutputNoneReads(vNoneReads, strUnalignedGenome)
       
def ReadBAM(strBAMFile):
    bamfile = pysam.AlignmentFile(strBAMFile, "rb")
    header = bamfile.header.copy()
    print(header)
    #print(type(bamfile.fetch(until_eof=True)))
    dicAlignment = {}
    for reads in bamfile.fetch(until_eof=True):
        objReads = ClsReads()
        objReads.Init(reads)
        if not objReads.strName in dicAlignment.keys():
            dicAlignment[objReads.strName] = []
            dicAlignment[objReads.strName].append(objReads)
        else:
            bDup = False
            for existRead in dicAlignment[objReads.strName]:
                if existRead.strCIGAR == objReads.strCIGAR:
                    bDup = True
                    break
            if not bDup:
                dicAlignment[objReads.strName].append(objReads)
        
    print(len(dicAlignment))
    AlignAnalyze(dicAlignment)
    
    
def main():
    fileName = "/home/lixin/lxwg/ad-hoc/H_Pylori/Data/Merged.all.sorted.bam"
    ReadBAM(fileName)

if __name__ == "__main__":
    main()