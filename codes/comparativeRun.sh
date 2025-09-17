#!/usr/bin/env bash

#SBATCH --job-name="nalif_00"
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

source /mgpfs/apps/bioinformatics/apps/miniforge3/24.3.0-0/etc/profile.d/conda.sh
source /mgpfs/apps/bioinformatics/apps/miniforge3/24.3.0-0/etc/profile.d/mamba.sh
mamba activate kpenv

python3 comparativeAnalysis.py