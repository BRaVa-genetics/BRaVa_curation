#!/usr/bin/env bash
#SBATCH --job-name=variant-meta-analysis
#SBATCH --output=/well/lindgren/dpalmer/logs/%x_%j.log
#SBATCH --error=/well/lindgren/dpalmer/logs/%x_%j.err

module load BWA/0.7.17-GCCcore-12.3.0
module load SAMtools/1.18-GCC-12.3.0
export PATH="$HOME/bin:$PATH"
export BCFTOOLS_PLUGINS="$HOME/bin"

FILES_VARIANT=$1
OUT=$2

REGION_FILE="/well/lindgren/dpalmer/protein_coding_regions_hg38_no_padding_no_UTR_v39.bed"

echo "bcftools +metal -Oz -o $OUT -e 'AF>0.01 & AF<0.99' -R $REGION_FILE $FILES_VARIANT"
bcftools +metal -Oz -o $OUT -e 'AF>0.01 & AF<0.99' -R $REGION_FILE $FILES_VARIANT
