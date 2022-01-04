### The Meaning of "-U"
To filter out specific regions from a BAM file, you could use the -U option of samtools view:

samtools view -b -L specificRegions.bed -U myFileWithoutSpecificRegions.bam myFile.bam > overlappingSpecificRegions.bam

-b Output in the BAM format
-L FILE Only output alignments overlapping the input BED FILE
-U FILE Write alignments that are not selected by the various filter options to FILE. When this option is used, all alignments intersecting the regions specified are written to either the output file or this file, but never both
