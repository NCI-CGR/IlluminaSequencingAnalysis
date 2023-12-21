'''
This is for H-Pylori Project 
Find the potential methylation
We assume the methylation comes from the disorded mapped reads (genes)
'''

import sys
import os
import subprocess
import gzip
import pysam
import csv
#import ast
from Bio import SeqIO

LSTRefName  = ["26695-1MET", "26695", "ELS37", "F57"]
SUFFIXRef   = "fsa"
SUFFIXReads = "ffa.gz"

DIRRef                   = "/scratch/lix33/lxwg/Project/H_pylori/Data/to-Xin/ref"
DIROutputBAMRoot         = "/scratch/lix33/lxwg/Project/H_pylori/Processed/BAM"
DIROutputLogRoot         = "/scratch/lix33/lxwg/Project/H_pylori/Processed/Log"
DIROutputFlagRoot        = "/scratch/lix33/lxwg/Project/H_pylori/Processed/Flag"
DIRRawSample             = "/scratch/lix33/lxwg/Project/H_pylori/Data/to-Xin/1012set"
DIRGFFV3                 = "/scratch/lix33/lxwg/Project/H_pylori/Data/to-Xin/GFFV3"
DIRUpdatedRawSample      = "/scratch/lix33/lxwg/Project/H_pylori/Processed/UpdatedRawSample"
DIROutputMethylationRoot = "/scratch/lix33/lxwg/Project/H_pylori/Processed/Methylation"



# DIRRawSample             = "/home/lixin/lxwg/Data/H-Pylori/RawReads"
# DIRUpdatedRawSample      = "/home/lixin/lxwg/Data/H-Pylori/UpdatedRawSample"
# DIROutputBAMRoot         = "/home/lixin/lxwg/Data/H-Pylori/BAM"
# DIROutputMethylationRoot = "/home/lixin/lxwg/Data/H-Pylori/Methylation"
# DIRGFFV3                 = "/home/lixin/lxwg/Data/H-Pylori/GFFV3" 
# DIRRef                   = "/home/lixin/lxwg/Data/H-Pylori/Ref"

CORE = "2"
BASHScript = "/scratch/lix33/lxwg/SourceCode/H_Pylori/Alignment.sh"

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

DICCIGAROFFSETREF = {"M": 1,
                     "I": 0,
                     "D": 1,
                     "N": 1,
                     "S": 0,
                     "H": 0,
                     "P": 0,
                     "=": 1,
                     "X": 1,
                     "B": 1}

DICCIGAROFFSETREADS = {"M": 1,
                       "I": 1,
                       "D": 0,
                       "N": 1,
                       "S": 1,
                       "H": 1,
                       "P": 1,
                       "=": 1,
                       "X": 1,
                       "B": 1}


class ClsRef:
    def __init__(self):
        self.dicInfo = {}
    
    def Init(self, strRefFile):
        fasta_sequences = SeqIO.parse(open(strRefFile),'fasta')       
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            self.dicInfo[name] = sequence

class ClsReads:
    def __init__(self):
        self.strName = ""
        self.strSeq = ""
        self.iIndex = 0
        # For Alignment  Result --> 
        self.iPos = ""
        self.strCIGAR = ""
        self.listCIGAR = []
        self.iAlignedLen = ""
        self.iReadsLength = ""
        self.strReadsSeq = ""
        self.strOrgName = ""
        self.strRefName = ""
        # <--
        self.strRawReadsName = ""
        self.iRawStartPos = -1
        self.iRawEndPos = -1
        
        self.vMathyGFFV3Pos = []
        self.vMathyPos = []
        self.vOrgGFFV3MathyPos = []
        self.vMathyMask = []
        self.vMathyStrand = []
        
        self.vMathyPosOffSetInGene = []
        self.vMathyPosRef = []
        self.vMathyCIGAR = []
        self.vMathyRefBase = []
        
        #tripleSeq for both ref and reads
        self.vTripleSeqRef = []
        self.vTripleSeqReads = []
        
        #Add two additional variable
        self.bIsReverse = False
        self.strRefSeq = ""
        self.vAdjustOffsetInRef = []
        
        self.strReadsFFAStrand = ""
    
    def UpdateReadsName(self):
        if self.strName == "":
            return
        vItem = self.strName.split(' ')
        vItem[0] += "|" + str(self.iIndex)
        self.strName = ' '.join(vItem)
    
    def Init(self, reads):
        self.strName = reads.query_name 
        self.strSeq = reads.get_forward_sequence() #reads.query_sequence
        self.strRefSeq = reads.get_reference_sequence()
        self.iPos = reads.reference_start
        self.strCIGAR = reads.cigarstring
        self.listCIGAR = reads.cigartuples
        self.iReadsLength = reads.query_length
        self.iAlignedLen = reads.reference_length
        self.strOrgName = '|'.join(self.strName.split('|')[:-1])
        self.strRefName = reads.reference_name
        self.bIsReverse = reads.is_reverse
        
        #print("Referene Name:", self.strRefName)
        
    def GetReadsMathyPos(self, stTmpStrand, strGFFV3Strand, iRawStartPos, iRawEndPos, iGFFV3MathyPos):
        if strGFFV3Strand == "+":
            if stTmpStrand == '+':
                return iGFFV3MathyPos
            else:
                return iRawEndPos - (iGFFV3MathyPos - iRawStartPos)
        else: # for the case in the - strand
            if stTmpStrand == '-':
                return iGFFV3MathyPos
            else:
                return iRawEndPos - (iGFFV3MathyPos - iRawStartPos)
    
    def GetRawReadsInfo(self, strRawReads, vMathyPos, vMathyMask, vMathyStrand, objRef):
        print(strRawReads)
        strPatten = self.strOrgName
        strPatten = strPatten.replace("[", "\[" )
        strPatten = strPatten.replace("]", "\]" )
        CMD = "zcat " + strRawReads + " | grep \"" + strPatten + "\"" 
        #print(CMD)
        self.strRawReadsName = subprocess.getoutput(CMD)
        #print(self.strRawReadsName)
        strTmp = ""
        stTmpStrand = ""
        if "Reverse" in self.strRawReadsName:
            strTmp = self.strRawReadsName.split(' Reverse')[0].split('|')[-1]
            stTmpStrand = "-"
        elif "Forward" in self.strRawReadsName:
            strTmp = self.strRawReadsName.split(' Forward')[0].split('|')[-1]
            stTmpStrand = "+"
        print(strTmp)
        self.iRawStartPos = int(strTmp.split(':')[0])
        self.iRawEndPos = int(strTmp.split(':')[1])
        self.strReadsFFAStrand = stTmpStrand
        # Get the related Mathylation
        iIndex = 0
        for strPos in vMathyPos:
            if (int(strPos) >= self.iRawStartPos and
                int(strPos) <= self.iRawEndPos):
                self.vMathyGFFV3Pos.append(int(strPos))
                iAdjustMethyPos = self.GetReadsMathyPos(stTmpStrand, vMathyStrand[iIndex], self.iRawStartPos, self.iRawEndPos, int(strPos))
                self.vMathyPos.append(iAdjustMethyPos)
                self.vOrgGFFV3MathyPos.append(int(strPos))
                self.vMathyMask.append(vMathyMask[iIndex])
                self.vMathyStrand.append(vMathyStrand[iIndex])
                iOffset = int(iAdjustMethyPos) - self.iRawStartPos
                self.vMathyPosOffSetInGene.append(iOffset)
                
                strReadsTripleSeq = ""
                if iOffset-1 < 0:
                    if iOffset + 2 <= len(self.strSeq):
                        strReadsTripleSeq = "#" + self.strSeq[iOffset:iOffset+2]
                    else:
                        strReadsTripleSeq = "#" + self.strSeq[iOffset] + "#"
                else:
                    if iOffset + 2 <= len(self.strSeq):
                        strReadsTripleSeq = self.strSeq[iOffset-1:iOffset+2]
                    else:
                        strReadsTripleSeq = self.strSeq[iOffset-1:iOffset+1] + "#" 
                self.vTripleSeqReads.append(strReadsTripleSeq)
                
                # #print(self.vMathyPos)
                # #print(self.vMathyPosOffSetInGene)
                #
                # # print("vMathyPos:", len(vMathyPos))
                # # print("vMathyMask:", len(vMathyMask))
                # # print("vMathyStrand:", len(vMathyStrand))
                # # print("iIndex", iIndex)
                # # print(strRawReads)
                # # print()
                # print("strOrgName", self.strOrgName.strip('\"'))
                # print("iPos                 :", self.iPos)
                # print("iRawStartPos         :", self.iRawStartPos)
                # print("iRawEndPos           :", self.iRawEndPos)
                # print("vMathyPos            :", self.vMathyPos)
                # print("vMathyPosOffSetInGene:", self.vMathyPosOffSetInGene)
                # print("self.bIsReverse      :", self.bIsReverse)
                # print("self.strSeq   :", self.strSeq)
                # print()
                # print("self.strRefSeq:", self.strRefSeq)
                # print()
                # print("============================")
                # print()
                # #break
                
            # Index only need to be increased here -->     
            if int(strPos) > self.iRawEndPos:
                break
            else:
                iIndex += 1
            # <--
        
        # Get Mathylation Ref Pos based on the case of RC
        self.GetMathylationRefPos(objRef)
        
        # print("self.vMathyPos            :", self.vMathyPos)
        # print("self.vMathyPosOffSetInGene:", self.vMathyPosOffSetInGene)
        # print("self.vTripleSeqReads      :", self.vTripleSeqReads)
        # print("ReadsLength               :", len(self.strSeq))
        # print()
        
        print("strOrgName", self.strOrgName.strip('\"'))
        print("iPos                 :", self.iPos)
        print("iRawStartPos         :", self.iRawStartPos)
        print("iRawEndPos           :", self.iRawEndPos)
        print("vOrgGFFV3MathyPos    :", self.vOrgGFFV3MathyPos)
        print("vMathyPos            :", self.vMathyPos)
        print("vMathyPosOffSetInGene:", self.vMathyPosOffSetInGene)
        print("vAdjustOffsetInRef   :", self.vAdjustOffsetInRef)
        print("self.bIsReverse      :", self.bIsReverse)
        print("vMathyCIGAR    :", self.vMathyCIGAR)
        print("vMathyRefBase  :", self.vMathyRefBase)
        print("vTripleSeqReads:", self.vTripleSeqReads)
        print("vTripleSeqRef  :", self.vTripleSeqRef)
        print()
        print("self.strSeq    :", self.strSeq)
        print()
        print("self.strRefSeq :", self.strRefSeq)
        print()
        print("============================")
        print()
        #break
        
    def GetMathylationRefPos(self, objRef):
        #print(self.strCIGAR)
        # print(self.listCIGAR)
        # for cigar in self.listCIGAR:
        #     print(cigar[0])
        #     print(cigar[1])
        self.vAdjustOffsetInRef.clear()
        for iOffSet in self.vMathyPosOffSetInGene:
            iAdjustOffSet = iOffSet
            # Update the Offset for the case of RC (reverse complementary)　
            if self.bIsReverse:
                iAdjustOffSet = len(self.strSeq) - 1 - iOffSet
            self.vAdjustOffsetInRef.append(iAdjustOffSet)
            
        
            iPosLeft, strRefBaseLeft, strCIGARLeft = self.GetSingleRefSeqInfo(objRef, iAdjustOffSet-1)
            # This is info for current Methylation -->
            iPosMiddle, strRefBaseMiddle, strCIGARMiddle = self.GetSingleRefSeqInfo(objRef, iAdjustOffSet)
            # <--
            iPosRight, strRefBaseRight, strCIGARRight = self.GetSingleRefSeqInfo(objRef, iAdjustOffSet+1)
            
            self.vMathyPosRef.append(iPosMiddle)
            self.vMathyCIGAR.append(strCIGARMiddle)
            self.vMathyRefBase.append(strRefBaseMiddle)
            # if self.bIsReverse:
            #     self.vTripleSeqRef.append(strRefBaseRight + strRefBaseMiddle + strRefBaseLeft)
            # else: # this is the normal string
            #     self.vTripleSeqRef.append(strRefBaseLeft + strRefBaseMiddle + strRefBaseRight)
            self.vTripleSeqRef.append(strRefBaseLeft + strRefBaseMiddle + strRefBaseRight)
              
    def GetSingleRefSeqInfo(self, objRef, iOffSet):
        iCount = 0
        iPos = self.iPos
        strCIGAR = ""
        for cigar in self.listCIGAR:
            if iCount + cigar[1] * DICCIGAROFFSETREADS[DICCIGAR[cigar[0]]] < iOffSet + 1:
                iPos += cigar[1] * DICCIGAROFFSETREF[DICCIGAR[cigar[0]]]
                iCount += cigar[1] * DICCIGAROFFSETREADS[DICCIGAR[cigar[0]]]
                strCIGAR += str(cigar[1]) + DICCIGAR[cigar[0]]
            else:
                iPos += (iOffSet - iCount + 1) * DICCIGAROFFSETREF[DICCIGAR[cigar[0]]]
                strCIGAR += str(iOffSet - iCount + 1) + DICCIGAR[cigar[0]]
                break
        strRefBase = ""
        if DICCIGAROFFSETREF[strCIGAR[-1]] == 0:
            strRefBase = "#"
            iPos = str(iPos) + "(#)"
        else:
            strRefBase = objRef.dicInfo[self.strRefName][iPos-1]
        return iPos, strRefBase, strCIGAR
    
    def GetInfo(self, listRow, iIndex):
        listRow.clear()
        listRow.append(self.strOrgName.strip('\"'))
        listRow.append(str(self.iPos))
        listRow.append(str(self.iRawStartPos))
        listRow.append(str(self.iRawEndPos))
        listRow.append(str(self.vMathyPos[iIndex]))
        listRow.append(str(self.vMathyPosOffSetInGene[iIndex]))
        # Two new important info --> 
        listRow.append(str(self.vAdjustOffsetInRef[iIndex]))
        listRow.append(self.bIsReverse)
        # <--
        listRow.append(str(self.vMathyPosRef[iIndex]))
        #-> new info
        listRow.append(str(self.vMathyMask[iIndex]))
        listRow.append(str(self.vMathyRefBase[iIndex]))
        listRow.append(str(self.vMathyStrand[iIndex]))
        listRow.append(str(self.strReadsFFAStrand))
        listRow.append(str(self.vMathyGFFV3Pos[iIndex]))
        #<-
        #-> Obtain triple sequence
        listRow.append(self.vTripleSeqRef[iIndex])
        listRow.append(self.vTripleSeqReads[iIndex])
        #<-
        
        listRow.append(str(self.vMathyCIGAR[iIndex]))

class ClsSample :
    def __init__(self):
        self.strName = ""
        self.strReads = "" # this is reads file
        self.strUpdatedReads = ""
        self.vReads = []
    
    def Init(self, strReadsPath):
        self.strName = strReadsPath.split('/')[-1].split('.')[0]
        self.strReads = strReadsPath
    
    def UpdateReads(self):
        # 1: Get File Name
        strFileName = os.path.basename(self.strReads)
        strUpdatedReads = DIRUpdatedRawSample + "/" + strFileName
        if not os.path.exists(strUpdatedReads):
            # Update reads Name (by adding index in the end of the reads)
            vReads = []
            objReads = ClsReads()
            #print(self.strReads)
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
    
    def GetOutofOrderReads(self, strRefName):
        strMappedBAMFile = DIROutputBAMRoot + "/" + strRefName + "/Mapped/" + self.strName + ".mapped.sorted.bam"
        strGeneListFile = DIROutputBAMRoot + "/" + strRefName + "/Mapped/" + self.strName + ".mapped.sorted.bam.gene.list"
        if not os.path.exists(strGeneListFile):
            print("GeneList does not exist! Return!")
            #print(strGeneListFile)
            return
        # Do the remaining jobs -> Go!
        # 0: Get the gene list
        CMD = "cat " + strGeneListFile
        vGeneList = subprocess.getoutput(CMD).split('\n')
        
        # 1: Get the org index list
        CMD = "awk -F '|' '{print $NF}' " + strGeneListFile
        vOrgIndex = subprocess.getoutput(CMD).split('\n')
        print(vOrgIndex)
        print()
        
        # 2: Get the sorted list (not continuous, since some reads may be unmapped!)
        CMD = "awk -F '|' '{print $NF}' " + strGeneListFile + " | sort -n"
        vSortedIndex = subprocess.getoutput(CMD).split('\n')
        print(vSortedIndex)
        print()
        
        # 3: Get the sequenced trunk one by one for vOrgIndex
        # Notice: this trunk may be 
        # (1) in order
        # (2) in reverse order
        vTrunk = []
        vItem = []
        for item in vOrgIndex:
            if len(vItem) == 0:
                vItem.append(item)
                continue
            if (int(vSortedIndex[vSortedIndex.index(vItem[-1]) + 1]) == int(item) or 
                int(vSortedIndex[vSortedIndex.index(vItem[-1]) - 1]) == int(item)):
                vItem.append(item)
                continue
            else:
                vTrunk.append(vItem)
                vItem = []
                vItem.append(item)
        if len(vItem) > 0:
            vTrunk.append(vItem)
        
        print(vTrunk)
        print()
        
        for item in vTrunk:
            print(item)
        
        print(strGeneListFile)
        
        #pick up the most obvious disordered reads
        vDisorderResult = []
        iCount = 0
        for i in range(0, len(vTrunk)):
            if len(vTrunk[i]) > 1:
                continue
            iCount += 1
            # skip the boundary trunks 
            if i - 1 < 0 or i + 1 >= len(vTrunk):
                continue
            if (int(vSortedIndex[vSortedIndex.index(vTrunk[i-1][-1]) + 1]) == int(vTrunk[i+1][0]) or 
                int(vSortedIndex[vSortedIndex.index(vTrunk[i-1][-1]) + 1]) == int(vTrunk[i+1][-1]) or 
                int(vSortedIndex[vSortedIndex.index(vTrunk[i-1][0]) + 1]) == int(vTrunk[i+1][0]) or
                int(vSortedIndex[vSortedIndex.index(vTrunk[i-1][0]) + 1]) == int(vTrunk[i+1][-1])):
                vDisorderResult.append(vTrunk[i][0])
            #break
        print(len(vTrunk))
        print(iCount)
        print(len(vDisorderResult))
        print(vDisorderResult)
        
        # Get the GeneName
        vDisorderedGeneList = []
        for item in vDisorderResult:
            vDisorderedGeneList.append(vGeneList[vOrgIndex.index(item)])
        for item in vDisorderedGeneList:
            print(item)
        
        strDisorderGeneListFile =  strGeneListFile + ".disorder"
        ofs = open(strDisorderGeneListFile, 'w')
        for item in vDisorderedGeneList:
            ofs.write(item + "\n")
        ofs.close()
    
    
    def GetMethylationCoordinate(self, strRefName, objRef):
        strMappedBAMFile = DIROutputBAMRoot + "/" + strRefName + "/Mapped/" + self.strName + ".mapped.sorted.bam"
        strGeneListFile = DIROutputBAMRoot + "/" + strRefName + "/Mapped/" + self.strName + ".mapped.sorted.bam.gene.list"
        strGFFV3File = DIRGFFV3 + "/" + self.strName + ".gff3.gz"
        strRawReads = DIRRawSample + "/" + self.strName + "." + SUFFIXReads
        #print(strGFFV3File)
        if (not os.path.exists(strGeneListFile) or 
            not os.path.exists(strGFFV3File) or 
            not os.path.exists(strRawReads)):
            #print("GeneList does not exist! Return!")
            #print(strGeneListFile)
            return
        print(strGFFV3File)
        # Know check each reads
        # Check aligned results one by one
        # 1: Scan aligned reads 
        # 2: Get the reads raw region 
        # 3: Get it's associated Mathylation
        # 4: Get the correct coordinate for each Mathylation
        # 5: Save the result to csv file
        # Go -> Parse current BAM file
        self.ParseBAM(strMappedBAMFile, strRawReads, strGFFV3File, objRef)
    
    def ParseBAM(self, strMappedBAMFile, strRawReads, strGFFV3File, objRef):
        # 1: Get the Methylation Pos
        CMD = "zcat " + strGFFV3File + " | tail -n +4 | awk '{print $4}'" #| sort -n"
        vMathyPos = subprocess.getoutput(CMD).split('\n')
        
        CMD = "zcat " + strGFFV3File + " | tail -n +4 | awk '{print $3}'" # | sort -n"
        vMathyMask = subprocess.getoutput(CMD).split('\n')
        
        CMD = "zcat " + strGFFV3File + " | tail -n +4 | awk '{print $7}'" # | sort -n"
        vMathyStrand = subprocess.getoutput(CMD).split('\n')
        #print(vMathyPos)
        #return 
        
        bamfile = pysam.AlignmentFile(strMappedBAMFile, "rb")
        header = bamfile.header.copy()
        #print(header)
        #print(type(bamfile.fetch(until_eof=True)))
        dicAlignment = {}
        self.vReads = []
        #iCount = 0
        for reads in bamfile.fetch(until_eof=True):
            objReads = ClsReads()
            objReads.Init(reads)
            #print(objReads.strOrgName)
            objReads.GetRawReadsInfo(strRawReads, vMathyPos, vMathyMask, vMathyStrand, objRef)
            #print(objReads.strRawReadsName)
            self.vReads.append(objReads)
            # #break
            # iCount += 1
            # if iCount == 5:
            #     break
            # break
        
        print(len(self.vReads))
    
    def PrintMethylationCoordinate(self, strRefName):
        strMappedBAMFile = DIROutputBAMRoot + "/" + strRefName + "/Mapped/" + self.strName + ".mapped.sorted.bam"
        strGeneListFile = DIROutputBAMRoot + "/" + strRefName + "/Mapped/" + self.strName + ".mapped.sorted.bam.gene.list"
        strGFFV3File = DIRGFFV3 + "/" + self.strName + ".gff3.gz"
        strRawReads = DIRRawSample + "/" + self.strName + "." + SUFFIXReads
        #print(strGFFV3File)
        if (not os.path.exists(strGeneListFile) or 
            not os.path.exists(strGFFV3File) or 
            not os.path.exists(strRawReads)):
            #print("GeneList does not exist! Return!")
            #print(strGeneListFile)
            return
        
        
        strCSVDir = DIROutputMethylationRoot + "/" + strRefName
        if not os.path.exists(strCSVDir):
            CMD = "mkdir -p " + strCSVDir
            os.system(CMD)
        strCSVFile = strCSVDir + "/" + self.strName + ".Methylation.list.csv"
        print("Preparing CSV File:", strCSVFile)
        f = open(strCSVFile, 'w')
        writer = csv.writer(f)
        listRow = ["GeneName", "AlignedRefPos(Gene)", "RawStartPosRef(Gene)", "RawEndPosRef(Gene)", "RawPosRef(Methylation)", "OffSetPosGene(Methylation)", "OffSetPosRef(Methylation)", "IsReverseComplement", "LiftoverAlignedPosRef(Methylation)", "MethyMask", "RefBase", "Strand(GFF3)", "Strand(FFA)", "MethyPos(GFF3)", "TripleSeqRef", "TripleSeqReads", "CIGAR(Calculation)"]
        writer.writerow(listRow)
        
        for reads in self.vReads:
            iMethyNum = len(reads.vMathyPosRef)
            #print(iMethyNum)
            for i in range(0, iMethyNum):
                reads.GetInfo(listRow, i)
                writer.writerow(listRow)
        f.close()
        print("Finished!")
        print()
        
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
        # update the reads name by binding "index" and grab the reads sequence
        sample.UpdateReads()
    
    # # Prepare Aligned BAM
    # for strRefName in LSTRefName:
    #     # Ref File 
    #     strRefFile = DIRRef + "/" + strRefName + "." + SUFFIXRef
    #     # Output BAM Dir 
    #     strOutputBAMDir = DIROutputBAMRoot + "/" + strRefName
    #     CMD = "mkdir -p " + strOutputBAMDir
    #     os.system(CMD) 
    #     # Output Log Dir
    #     strOutputLogDir = DIROutputLogRoot + "/" + strRefName
    #     CMD = "mkdir -p " + strOutputLogDir
    #     os.system(CMD)
    #     # Output Log Dir
    #     strOutputFlagDir = DIROutputFlagRoot + "/" + strRefName
    #     CMD = "mkdir -p " + strOutputFlagDir
    #     os.system(CMD)
    #
    #     for sample in vSample:
    #         sample.SubmitJob(strRefName, strRefFile, strOutputBAMDir, strOutputLogDir, strOutputFlagDir)
    
    # # Merge BAM for each reference respectively
    # for strRefName in LSTRefName:
    #     strOutputBAMDir = DIROutputBAMRoot + "/" + strRefName
    #     strOutputFlagDir = DIROutputFlagRoot + "/" + strRefName
    #     # MergeBAM
    #     MergeBAM(strRefName, strOutputBAMDir, strOutputFlagDir)
        
    # # Get the out of order Reads From Each BAM
    # # Prepare Aligned BAM
    # for strRefName in LSTRefName:
    #     for sample in vSample:
    #         sample.GetOutofOrderReads(strRefName)
    #         break
    
    # Get the Correct Methylation coordinate From Each BAM
    # Prepare Aligned BAM
    
    for strRefName in LSTRefName:
        # Read current reference 
        strRefFile = DIRRef + "/" + strRefName + "." + SUFFIXRef
        objRef = ClsRef()
        objRef.Init(strRefFile)
        #print(strRefFile)
        #print(objRef.dicInfo)
        for sample in vSample:
            print(sample.strName)
            sample.GetMethylationCoordinate(strRefName, objRef)
            # Print Result 
            sample.PrintMethylationCoordinate(strRefName)
            # break
        #break
        
if __name__ == "__main__":
    main()