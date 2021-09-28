#Version 4 (20190411)

The new features:
1. The new @RG in the BAMs are changed as:
@RG      ID: AHTLY7BBXX.8            SM: SC310221    LB: SC310221_ATGCCTAA             PL:ILLUMINA     PU: AHTLY7BBXX20180606.8. ATGCCTAA      CN:CGR
Remove PU, LB in the body part 
The benefits are: saving storage space; more precisely for BQSR and dedup processing.

2. Skipping some unnessary steps, removing phix147 (time consuming), (most of them have been removed in the primary analysis scripts; picard reorderSam can do that even phix147 not removed); consolidate to one line @RG (will lost information in reclibrated BAMs); header check (skip validate_or_update_bam_header.sh);
3. Picard version upgrade and using VALIDATION_STRINGENCY=STRICT 

Usage:
1. For this new format, in step2, if the lane-level BAMs are old ones, the script will reformat the old lane-level BAMs and generate the new BAM suffixed as "_reformated.bam" in the flowcell folder. If the lane-level BAMs have been archieved, can grep them from log stdout files: grep lane /DCEG/Projects/Exome/SequencingData/variant_scripts/logs/GATK/patch_build_bam_${date}/*.stdout | cut -d" " -f6
2. Add a new script step7a_reformat_bams.sh $MANIFEST $RESTORE_FILE, if the large build, the merging record has not been changed but the subject-level BAMs are still in old format, the script will modify that according the information in the Manifest file.

