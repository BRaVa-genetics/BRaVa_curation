#!/usr/bin/env bash
#SBATCH --job-name=cauchy_combine
#SBATCH --output=/well/lindgren/dpalmer/logs/%x_%j.log
#SBATCH --error=/well/lindgren/dpalmer/logs/%x_%j.err

FILE=$1
OUT=$2

Rscript cauchy_combination_of_results.r --results_file $FILE --out $OUT
