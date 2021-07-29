#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mail-type=END
#SBATCH --mail-user=lix33@nih.gov
#SBATCH --output=std.out
#SBATCH --error=std.err
#SBATCH --job-name=ObjStorage
#SBATCH --time=10-00:00:00

module load python/3.7

python3 --version

#sourceCode="/home/lix33/lxwg/Git/IlluminaSequencingAnalysis/CustomizedQC/SourceCode/ObjectStorage/Test.py"
sourceCode="/home/lix33/lxwg/Git/IlluminaSequencingAnalysis/CustomizedQC/SourceCode/ObjectStorage/Backup2S3.py"

python3 ${sourceCode}

echo "sbatch end!"
