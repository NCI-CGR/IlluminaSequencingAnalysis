#!/bin/sh
SCRIPT=$(readlink -f "$0")
SCRIPT_HOME=$(dirname "$SCRIPT")
. ${SCRIPT_HOME}/global_config_bash.rc
if [ ! -d "$GATK_LOG_BUILD_BAM_DIR" ]; then
mkdir $GATK_LOG_BUILD_BAM_DIR
fi
CMD="qsub -q $QUEUE -o ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA12878_GDNA_HS4K_KHPL_1.stdout -e ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA12878_GDNA_HS4K_KHPL_1.stderr -N BAMCOMB.CTRL_NA12878_GDNA_HS4K_KHPL_1 -v SCRIPT_HOME=$SCRIPT_HOME -S /bin/sh ${SCRIPT_HOME}/gatk_build_bam_for_single_name_v4.sh CTRL_NA12878_GDNA_HS4K_KHPL_1 1 /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160613_K00278_0053_AH7FTFBBXX/BAM/SB887769_1_CAGATC_L007.bam /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160613_K00278_0053_AH7FTFBBXX/BAM/SB887769_CGATGT_L007.bam "
echo $CMD
eval $CMD
CMD="qsub -q $QUEUE -o ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA12878_GDNA_HS4K_KHPL_2.stdout -e ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA12878_GDNA_HS4K_KHPL_2.stderr -N BAMCOMB.CTRL_NA12878_GDNA_HS4K_KHPL_2 -v SCRIPT_HOME=$SCRIPT_HOME -S /bin/sh ${SCRIPT_HOME}/gatk_build_bam_for_single_name_v4.sh CTRL_NA12878_GDNA_HS4K_KHPL_2 1 /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160613_K00278_0053_AH7FTFBBXX/BAM/SB887769_1_CAGATC_L008.bam /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160613_K00278_0053_AH7FTFBBXX/BAM/SB887769_CGATGT_L008.bam "
echo $CMD
eval $CMD
CMD="qsub -q $QUEUE -o ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA12878_GDNA_HSV4_KHPL_1.stdout -e ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA12878_GDNA_HSV4_KHPL_1.stderr -N BAMCOMB.CTRL_NA12878_GDNA_HSV4_KHPL_1 -v SCRIPT_HOME=$SCRIPT_HOME -S /bin/sh ${SCRIPT_HOME}/gatk_build_bam_for_single_name_v4.sh CTRL_NA12878_GDNA_HSV4_KHPL_1 1 /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160607_D00620_0070_BC6PWPANXX/BAM/SB887769_1_CAGATC_L001.bam /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160607_D00620_0070_BC6PWPANXX/BAM/SB887769_CGATGT_L001.bam "
echo $CMD
eval $CMD
CMD="qsub -q $QUEUE -o ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA12878_GDNA_HSV4_KHPL_2.stdout -e ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA12878_GDNA_HSV4_KHPL_2.stderr -N BAMCOMB.CTRL_NA12878_GDNA_HSV4_KHPL_2 -v SCRIPT_HOME=$SCRIPT_HOME -S /bin/sh ${SCRIPT_HOME}/gatk_build_bam_for_single_name_v4.sh CTRL_NA12878_GDNA_HSV4_KHPL_2 1 /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160607_D00620_0070_BC6PWPANXX/BAM/SB887769_1_CAGATC_L002.bam /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160607_D00620_0070_BC6PWPANXX/BAM/SB887769_CGATGT_L002.bam "
echo $CMD
eval $CMD
CMD="qsub -q $QUEUE -o ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24143_GDNA_HS4K_KHPL_1.stdout -e ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24143_GDNA_HS4K_KHPL_1.stderr -N BAMCOMB.CTRL_NA24143_GDNA_HS4K_KHPL_1 -v SCRIPT_HOME=$SCRIPT_HOME -S /bin/sh ${SCRIPT_HOME}/gatk_build_bam_for_single_name_v4.sh CTRL_NA24143_GDNA_HS4K_KHPL_1 1 /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160613_K00278_0053_AH7FTFBBXX/BAM/SB988949_1_ATCACG_L007.bam /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160613_K00278_0053_AH7FTFBBXX/BAM/SB988949_ACAGTG_L007.bam "
echo $CMD
eval $CMD
CMD="qsub -q $QUEUE -o ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24143_GDNA_HS4K_KHPL_2.stdout -e ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24143_GDNA_HS4K_KHPL_2.stderr -N BAMCOMB.CTRL_NA24143_GDNA_HS4K_KHPL_2 -v SCRIPT_HOME=$SCRIPT_HOME -S /bin/sh ${SCRIPT_HOME}/gatk_build_bam_for_single_name_v4.sh CTRL_NA24143_GDNA_HS4K_KHPL_2 1 /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160613_K00278_0053_AH7FTFBBXX/BAM/SB988949_1_ATCACG_L008.bam /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160613_K00278_0053_AH7FTFBBXX/BAM/SB988949_ACAGTG_L008.bam "
echo $CMD
eval $CMD
CMD="qsub -q $QUEUE -o ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24143_GDNA_HSV4_KHPL_1.stdout -e ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24143_GDNA_HSV4_KHPL_1.stderr -N BAMCOMB.CTRL_NA24143_GDNA_HSV4_KHPL_1 -v SCRIPT_HOME=$SCRIPT_HOME -S /bin/sh ${SCRIPT_HOME}/gatk_build_bam_for_single_name_v4.sh CTRL_NA24143_GDNA_HSV4_KHPL_1 1 /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160607_D00620_0070_BC6PWPANXX/BAM/SB988949_1_ATCACG_L001.bam /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160607_D00620_0070_BC6PWPANXX/BAM/SB988949_ACAGTG_L001.bam "
echo $CMD
eval $CMD
CMD="qsub -q $QUEUE -o ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24143_GDNA_HSV4_KHPL_2.stdout -e ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24143_GDNA_HSV4_KHPL_2.stderr -N BAMCOMB.CTRL_NA24143_GDNA_HSV4_KHPL_2 -v SCRIPT_HOME=$SCRIPT_HOME -S /bin/sh ${SCRIPT_HOME}/gatk_build_bam_for_single_name_v4.sh CTRL_NA24143_GDNA_HSV4_KHPL_2 1 /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160607_D00620_0070_BC6PWPANXX/BAM/SB988949_1_ATCACG_L002.bam /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160607_D00620_0070_BC6PWPANXX/BAM/SB988949_ACAGTG_L002.bam "
echo $CMD
eval $CMD
CMD="qsub -q $QUEUE -o ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24149_GDNA_HS4K_KHPL_1.stdout -e ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24149_GDNA_HS4K_KHPL_1.stderr -N BAMCOMB.CTRL_NA24149_GDNA_HS4K_KHPL_1 -v SCRIPT_HOME=$SCRIPT_HOME -S /bin/sh ${SCRIPT_HOME}/gatk_build_bam_for_single_name_v4.sh CTRL_NA24149_GDNA_HS4K_KHPL_1 1 /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160613_K00278_0053_AH7FTFBBXX/BAM/SB988955_1_TTAGGC_L007.bam /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160613_K00278_0053_AH7FTFBBXX/BAM/SB988955_GCCAAT_L007.bam "
echo $CMD
eval $CMD
CMD="qsub -q $QUEUE -o ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24149_GDNA_HS4K_KHPL_2.stdout -e ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24149_GDNA_HS4K_KHPL_2.stderr -N BAMCOMB.CTRL_NA24149_GDNA_HS4K_KHPL_2 -v SCRIPT_HOME=$SCRIPT_HOME -S /bin/sh ${SCRIPT_HOME}/gatk_build_bam_for_single_name_v4.sh CTRL_NA24149_GDNA_HS4K_KHPL_2 1 /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160613_K00278_0053_AH7FTFBBXX/BAM/SB988955_1_TTAGGC_L008.bam /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160613_K00278_0053_AH7FTFBBXX/BAM/SB988955_GCCAAT_L008.bam "
echo $CMD
eval $CMD
CMD="qsub -q $QUEUE -o ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24149_GDNA_HSV4_KHPL_1.stdout -e ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24149_GDNA_HSV4_KHPL_1.stderr -N BAMCOMB.CTRL_NA24149_GDNA_HSV4_KHPL_1 -v SCRIPT_HOME=$SCRIPT_HOME -S /bin/sh ${SCRIPT_HOME}/gatk_build_bam_for_single_name_v4.sh CTRL_NA24149_GDNA_HSV4_KHPL_1 1 /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160607_D00620_0070_BC6PWPANXX/BAM/SB988955_1_TTAGGC_L001.bam /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160607_D00620_0070_BC6PWPANXX/BAM/SB988955_GCCAAT_L001.bam "
echo $CMD
eval $CMD
CMD="qsub -q $QUEUE -o ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24149_GDNA_HSV4_KHPL_2.stdout -e ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24149_GDNA_HSV4_KHPL_2.stderr -N BAMCOMB.CTRL_NA24149_GDNA_HSV4_KHPL_2 -v SCRIPT_HOME=$SCRIPT_HOME -S /bin/sh ${SCRIPT_HOME}/gatk_build_bam_for_single_name_v4.sh CTRL_NA24149_GDNA_HSV4_KHPL_2 1 /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160607_D00620_0070_BC6PWPANXX/BAM/SB988955_1_TTAGGC_L002.bam /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160607_D00620_0070_BC6PWPANXX/BAM/SB988955_GCCAAT_L002.bam "
echo $CMD
eval $CMD
CMD="qsub -q $QUEUE -o ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24385_GDNA_HS4K_KHPL_1.stdout -e ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24385_GDNA_HS4K_KHPL_1.stderr -N BAMCOMB.CTRL_NA24385_GDNA_HS4K_KHPL_1 -v SCRIPT_HOME=$SCRIPT_HOME -S /bin/sh ${SCRIPT_HOME}/gatk_build_bam_for_single_name_v4.sh CTRL_NA24385_GDNA_HS4K_KHPL_1 1 /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160613_K00278_0053_AH7FTFBBXX/BAM/SB988943_1_CTTGTA_L007.bam /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160613_K00278_0053_AH7FTFBBXX/BAM/SB988943_TGACCA_L007.bam "
echo $CMD
eval $CMD
CMD="qsub -q $QUEUE -o ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24385_GDNA_HS4K_KHPL_2.stdout -e ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24385_GDNA_HS4K_KHPL_2.stderr -N BAMCOMB.CTRL_NA24385_GDNA_HS4K_KHPL_2 -v SCRIPT_HOME=$SCRIPT_HOME -S /bin/sh ${SCRIPT_HOME}/gatk_build_bam_for_single_name_v4.sh CTRL_NA24385_GDNA_HS4K_KHPL_2 1 /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160613_K00278_0053_AH7FTFBBXX/BAM/SB988943_1_CTTGTA_L008.bam /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160613_K00278_0053_AH7FTFBBXX/BAM/SB988943_TGACCA_L008.bam "
echo $CMD
eval $CMD
CMD="qsub -q $QUEUE -o ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24385_GDNA_HSV4_KHPL_1.stdout -e ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24385_GDNA_HSV4_KHPL_1.stderr -N BAMCOMB.CTRL_NA24385_GDNA_HSV4_KHPL_1 -v SCRIPT_HOME=$SCRIPT_HOME -S /bin/sh ${SCRIPT_HOME}/gatk_build_bam_for_single_name_v4.sh CTRL_NA24385_GDNA_HSV4_KHPL_1 1 /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160607_D00620_0070_BC6PWPANXX/BAM/SB988943_1_CTTGTA_L001.bam /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160607_D00620_0070_BC6PWPANXX/BAM/SB988943_TGACCA_L001.bam "
echo $CMD
eval $CMD
CMD="qsub -q $QUEUE -o ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24385_GDNA_HSV4_KHPL_2.stdout -e ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_CTRL_NA24385_GDNA_HSV4_KHPL_2.stderr -N BAMCOMB.CTRL_NA24385_GDNA_HSV4_KHPL_2 -v SCRIPT_HOME=$SCRIPT_HOME -S /bin/sh ${SCRIPT_HOME}/gatk_build_bam_for_single_name_v4.sh CTRL_NA24385_GDNA_HSV4_KHPL_2 1 /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160607_D00620_0070_BC6PWPANXX/BAM/SB988943_1_CTTGTA_L002.bam /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data/160607_D00620_0070_BC6PWPANXX/BAM/SB988943_TGACCA_L002.bam "
echo $CMD
eval $CMD







