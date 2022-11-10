threads=8
RG="@RG\\tID:hpylori\\tSM:WES"
ref="/home/lix33/DAATeam_Xin/ad_hoc/H_Pylori/Data/HpGP-26695-ATCC.fsa"
sampleName="HpGP-26695-ATCC"
Sample="/home/lix33/DAATeam_Xin/ad_hoc/H_Pylori/Data/${sampleName}.ffa"
OutSAM="/home/lix33/DAATeam_Xin/ad_hoc/H_Pylori/Test/Output/${sampleName}.sam"
OutBAM="/home/lix33/DAATeam_Xin/ad_hoc/H_Pylori/Test/Output/${sampleName}.bam"
OutSortBAM="/home/lix33/DAATeam_Xin/ad_hoc/H_Pylori/Test/Output/${sampleName}.sorted.bam"
SortedIndex="${OutSortBAM}.bai"

#minimap2 -t ${threads} -Y --MD -R '${RG}' -a ${ref} ${SampleControl} > ${OutSAM}

minimap2 -t ${threads} -Y --MD -a ${ref} ${Sample} > ${OutSAM}

samtools view -@ ${threads} -S -b ${OutSAM} > ${OutBAM}
samtools sort -@ ${threads} -o ${OutSortBAM} ${OutBAM}
samtools index -@ ${threads} ${OutSortBAM} ${SortedIndex}
rm ${OutSAM}
rm ${OutBAM}


OutSortBAM="/home/lix33/DAATeam_Xin/ad_hoc/H_Pylori/Test/Output/Merged.sorted.bam"
SortedIndex="/home/lix33/DAATeam_Xin/ad_hoc/H_Pylori/Test/Output/Merged.sorted.bam.bai"

OutSubSortBAM1="/home/lix33/DAATeam_Xin/ad_hoc/H_Pylori/Test/Output/HpGP-26695-ATCC.sorted.bam"
OutSubSortBAM2="/home/lix33/DAATeam_Xin/ad_hoc/H_Pylori/Test/Output/HpGP-LAT-002.sorted.bam"
OutSubSortBAM3="/home/lix33/DAATeam_Xin/ad_hoc/H_Pylori/Test/Output/HpGP-SWT-002.sorted.bam"
samtools merge -@ {threads} -f ${OutSortBAM} ${OutSubSortBAM1} ${OutSubSortBAM2} ${OutSubSortBAM3}

threads=8
OutSortBAM="/home/lix33/DAATeam_Xin/ad_hoc/H_Pylori/Test/Output/Merged.sorted.bam"
OutBAM="/home/lix33/DAATeam_Xin/ad_hoc/H_Pylori/Test/Output/Merged.all.sorted.bam"
samtools sort -@ ${threads} -o ${OutBAM} ${OutSortBAM}
samtools index -@ ${threads} ${OutBAM}
rm ${OutSortBAM}



~~~~~~~~~~~~~~~~~~~~~~~~~~


+
>>>>>>>>>>>>
. /etc/profile.d/modules.sh
module purge

module load miniconda/4.8.3 sge

#. /home/lix33/.bashrc

#module load sge
#unset module
#source /scratch/lix33/lxwg/software/conda/miniconda3/condabin/activate envCGR

source /home/lix33/.conda/envs/envCGR/etc/profile.d/conda.sh
conda activate /home/lix33/.conda/envs/envCGR

#module load python3/3.8.12 samtools/1.8 tabix/1.9 bgzip/0.2.6 gcc/7.2.0 zlib/1.2.11 bedtools/2.27.1 python3/3.7.0 
module load samtools/1.8 tabix/1.9 bgzip/0.2.6 gcc/7.2.0 zlib/1.2.11 bedtools/2.27.1 python3/3.7.0 python3/3.8.12
