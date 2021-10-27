#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mail-type=END
#SBATCH --mail-user=lix33@nih.gov
#SBATCH --output=/data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Log/BackupCOVID2ndResults/std_10_27_2021_01_200_std_output.out
#SBATCH --error=/data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Log/BackupCOVID2ndResults/std_10_27_2021_01_200_std_output.err
#SBATCH --job-name=BackupCOVID2ndResults
#SBATCH --time=10-00:00:00

module load python/3.7

python3 --version

#sourceCode="/home/lix33/lxwg/Git/IlluminaSequencingAnalysis/CustomizedQC/SourceCode/ObjectStorage/Test.py"
sourceCode="/home/lix33/lxwg/Git/IlluminaSequencingAnalysis/CustomizedQC/SourceCode/ObjectStorage/BackupCOVID2ndResults/BackupCovid2ndPipelineResults.py"

python3 ${sourceCode}

echo "sbatch (BackupCovid2ndPipelineResults.py) end!"
