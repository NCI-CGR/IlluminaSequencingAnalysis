#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mail-type=END
#SBATCH --mail-user=lix33@nih.gov
#SBATCH --output=/data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Log/BackupFastq/std.out
#SBATCH --error=/data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Log/BackupFastq/std.err
#SBATCH --job-name=S3Fastq
#SBATCH --time=10-00:00:00

#module load python/3.7

#python3 --version

#sourceCode="/home/lix33/lxwg/Git/IlluminaSequencingAnalysis/CustomizedQC/SourceCode/ObjectStorage/Test.py"
sourceCode="/home/lix33/lxwg/Git/IlluminaSequencingAnalysis/CustomizedQC/SourceCode/ObjectStorage/BackupRawUSUFastq/BackupRawUSUFastq.sh"

bash ${sourceCode}

echo "sbatch end!"