### New Requirement from Cenk's research group
1. Original Email 
```
Hi Lisa and Meredith,
 
I hope your summers have been going well! I'd like to give you folks an update on whats been happening on our side:
•	Expanded ImmunoTyper to genotype:
o	IGLV - waiting on validated samples but have 95% allele-presence concordance on 9 WGS replicates
o	TRAV  - 95% precision / 99% recall for allele-presence on 10 validated individuals
•	Completed COVID severity association on the n=1108 NIAID cohort using IGHV allele, IGHV SNV, IGLV allele and IGLV SNV as independent variables
o	We were able to move away from meta-analysis because the NIAID folks provided us with PCAs that we're generated all together 
o	1 significant IGHV allele after FDR correction but it is a low confidence allele call
o	1 significant risk IGHV SNV after FDR correction
o	No significant IGLV associations
•	Performed IGHV, IGLV genotyping on the Chernobyl family of 4 we currently have - Mendelian concordance for the two offspring using allele-presence: 
  (0.968, 0.975) for IGHV, (0.99, 1.0) for IGLV
o	This excludes any notion of copies since we can't do haplotype phasing
We're currently working on expanding to the other TCR loci and IGK, and apply those to COVID severity association (recognizing that this will hurt the FDR 
correction). The NIAID team is unwilling to share WGS for TCR since they feel that overlaps with a TCR repertoire study they are doing on that dataset. 
Would you all be willing to share the WGS for those loci from the COVNET dataset? It would give us somewhere to start, even if the dataset is a bit messier. 
We can provide you a BED file with the new regions to extract. 
 
We'd also like to continue with the mendelian concordance analysis on the chernobyl sample - its just such a great dataset. We're ready for new WGS data 
to be put in the shared drive (/data/DCEG_Trios/Cenk - what's in there is good to be removed from our side) whenever you have time to do the transfer. 
```

2. What we need to do
  * Prepare the sliced BAMs from TCR loci and IGK. 
  * The bed file is below
     * Trio-Bed: [ImmunoTyper_regions_of_interest_GRCh38.zip](https://github.com/NCI-CGR/IlluminaSequencingAnalysis/files/9521866/ImmunoTyper_regions_of_interest_GRCh38.zip) 
     * Updated-Bed (used in program): [ImmunoTyper_regions_of_interest_GRCh38.new.zip](https://github.com/NCI-CGR/IlluminaSequencingAnalysis/files/9521865/ImmunoTyper_regions_of_interest_GRCh38.new.zip)

### Implementation
1. The implementation is almost identical with "https://github.com/NCI-CGR/IlluminaSequencingAnalysis/tree/main/Ad-hoc/Cenk"

2. the running steps is identical with "https://github.com/NCI-CGR/IlluminaSequencingAnalysis/tree/main/Ad-hoc/Cenk"  

### Results
1. All sliced BAM files based on “ImmunoTyper_regions_of_interest_GRCh38.bed” are located in: 
   * /data/COVID_ADHOC/Sequencing/COVID_WGS/ad-hoc/Cenk/TargetBAM_TCR_IGK
 
2. The delivered Sample list are located in: 
   * /data/COVID_ADHOC/Sequencing/COVID_WGS/ad-hoc/Cenk/UpdatedExcel/TargetBAM_TCR_IGK
3. The total number of delivered sample is 475. 
   * BTW, I just did a quick check for these sliced BAMs, and most of them only contain around 14 regions defined in your bed file (not all of them).  
