### MergeOverlapInBedFile.py ###
1. One input
example:
/home/lixin/lxwg/Test/CGR/UCSC/BedFileForRef38_CCDS


2. Output: 
There are two output, including:
(1) Full version of Merged Bed file
-rw-rw---- 1 lix33 ncicgf_dceg_exome 12427761 Aug 27 14:56 BedFileForRef38_CCDS.MergedOverlap.bed

(2)Brief version of Merged Bed file, which only contain the first 3 columns
-rw-rw---- 1 lix33 ncicgf_dceg_exome  4595387 Aug 27 14:56 BedFileForRef38_CCDS.MergedOverlap.Brief.bed

### BedDiffCheck.py ###
1: Input
(1) The first argument: 1st BED file: hg 19
    example:
    /home/lixin/lxwg/Test/CGR/UCSC/Test_Hg19/hg19_cds.MergedOverlap.Brief.bed
(2) The second argument: 2nd BED file: v38
    example:
    /home/lixin/lxwg/Test/CGR/UCSC/LiftOver/v38_to_hg19_liftover.bed
    
2: Output
statistic comparison result

For details, Please check https://github.com/NCI-CGR/IlluminaSequencingAnalysis/issues/33

