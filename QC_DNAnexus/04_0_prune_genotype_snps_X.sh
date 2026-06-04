#!/usr/bin/env bash

# UKB genotype calls
ukb_bed="/mnt/project/Bulk/Genotype Results/Genotype calls/ukb22418_cX_b0_v2.bed"
ukb_bim="/mnt/project/Bulk/Genotype Results/Genotype calls/ukb22418_cX_b0_v2.bim"
ukb_fam="/mnt/project/Bulk/Genotype Results/Genotype calls/ukb22418_cX_b0_v2.fam"

out="ukb_snp_chrX_pruned"

# Create the vcf fam
vcf_samples="UKBB_WES450K.txt"
cat /mnt/project/ukb_wes_450k_qc/data/ukb_wes_450k.qced.sample_ids.tsv | awk '{print $1, $1}' | tail -n +2 > $vcf_samples

./plink --bed "$ukb_bed" --bim "$ukb_bim" --fam "$ukb_fam" \
      --keep "$vcf_samples" --maf 0.05 --geno 0.02 \
      --indep 50 5 2 --out ${out}
