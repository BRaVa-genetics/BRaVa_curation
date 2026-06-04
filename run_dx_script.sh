#!/bin/bash

# "e" stops process if error occurs
# "u" throws error if variables are undefined
set -eu 

WD="/Users/barneyh/brava/BRaVa_curation/"

## Job parameters
# instance_type="mem1_ssd1_v2_x72" # Very large
# instance_type="mem3_ssd1_v2_x32" # Large
instance_type="mem3_ssd1_v2_x96" # DEFAULT
# instance_type="mem1_ssd1_v2_x8" # DEFAULT
# instance_type="mem1_ssd1_v2_x2" # Small


# instance_count="64" # Very large
instance_count="1"
# instance_count="16" # Large
# instance_count="8" # Medium
# instance_count="4" # Small-mid


# Duration in minutes
duration="1000" # arbitrarily large, cost limit will automatically terminate if needed

# Cost limit (in GBP)
cost_limit="500"

#pyscript="02_prefilter_variants.py"
#pyscript="03_0_initial_sample_qc.py"
pyscript="10_2_create_qc_plink.py"

dx upload QC_DNAnexus/$pyscript --destination /Barney/scripts/

for chrom in 19; do
    # Name of job
    name="${pyscript}_c${chrom}-${instance_count}x"
    dx run \
      dxjupyterlab_spark_cluster \
        --instance-type="${instance_type}" \
        --instance-count="${instance_count}" \
        -icmd="python3 /mnt/project/Barney/scripts/${pyscript} --chrom ${chrom}" \
        --name="${name}" \
        --priority="low" \
        --yes \
        -iduration="${duration}" \
        -ifeature="HAIL" \
        --brief \
        --cost-limit="${cost_limit}" \
        --extra-args='{"default_task_dx_attributes": {"runSpec": {"executionPolicy": {"restartOn": {"*": 3}}}}}'

  done
done
