#!/usr/bin/env bash
#SBATCH --job-name=variant-meta-analysis
#SBATCH --output=/well/lindgren/dpalmer/logs/%x_%j.log
#SBATCH --error=/well/lindgren/dpalmer/logs/%x_%j.err

module load BWA/0.7.17-GCCcore-12.3.0
module load SAMtools/1.18-GCC-12.3.0
module load R
export PATH="$HOME/bin:$PATH"
export BCFTOOLS_PLUGINS="$HOME/bin"

FILE=$1
OUT=$2

Rscript prepare_files_for_bcftools_metal.r --file_path $FILE --out_data_dir $OUT
