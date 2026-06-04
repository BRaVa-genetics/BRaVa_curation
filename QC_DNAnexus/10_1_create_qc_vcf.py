import hail as hl
import sys
import pyspark
import dxpy
import pandas as pd
from math import ceil
import json
import os
import re
import glob
import argparse

def export_file(path, out_folder):
    dxpy.upload_local_file(
        filename=path,
        name=path.split('/')[-1],
        folder=out_folder,
        parents=True
    )


def main(chrom):

    my_database = dxpy.find_one_data_object(
        name="my_database", 
        project=dxpy.find_one_project()["id"]
    )["id"]
    database_dir = f'dnax://{my_database}'
    sc = pyspark.SparkContext()
    spark = pyspark.sql.SparkSession(sc)
    hl.init(sc=sc, tmp_dir=f'{database_dir}/tmp/')


    def get_final_filter_mt_path(chrom):
        return f'{database_dir}/04_final_filter_v2_write_to_mt/ukb_wes_450k.qced.chr{chrom}.mt'

    # Inputs:
    SUPERPOPS = "file:///mnt/project/Barney/qc2/05_estimate_superpopulation/superpopulation_labels.tsv"
    FINAL_VARIANT_LIST = 'file:///mnt/project/Barney/qc2/08_0_final_variant_qc/08_final_qc.pop.keep.variant_list'
    FINAL_SAMPLE_LIST = 'file:///mnt/project/Barney/qc2/09_0_final_sample_qc/09_final_qc.keep.BRaVa.sample_list'

    # Outputs
    ht_final_samples = hl.import_table(FINAL_SAMPLE_LIST, no_header=True, key='f0', delimiter=',', types={'f0': hl.tstr})
    ht_final_variants = hl.import_table(FINAL_VARIANT_LIST, types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
    ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)
    ht_superpops = hl.import_table(SUPERPOPS, impute=True, types={'f0': hl.tstr}).key_by("sample.ID").select("classification_strict")

    # Loop over populations
    
    for pop in ["AFR", "AMR", "EAS", "EUR", "SAS"]:

        mt = hl.read_matrix_table(get_final_filter_mt_path(chrom))
        mt = mt.annotate_cols(pops = ht_superpops[mt.s])
        mt = mt.drop('qual', 'info', 'filters')

        mt = mt.filter_cols((mt.pops.classification_strict == pop))
        # Filter
        mt = mt.filter_cols(hl.is_defined(ht_final_samples[mt.col_key]))
        ht_final_variants_tmp = ht_final_variants.filter(ht_final_variants.pop==pop)
        mt = mt.filter_rows(hl.is_defined(ht_final_variants_tmp[mt.row_key]))

        mt = mt.select_entries(GT = hl.unphased_diploid_gt_index_call(mt.GT.n_alt_alleles()))

        vcf_file=f'ukb_wes_450k.qced.barney_qc.{pop}.chr{chrom}.vcf.bgz'
        hl.export_vcf(mt, output=vcf_file)

        os.system(f'hdfs dfs -get {vcf_file} .')

        export_file(
            path=vcf_file, 
            out_folder='/Barney/qc/10_1_create_qc_vcf/'
		)

if __name__=='__main__':

    parser = argparse.ArgumentParser(description="Process some integers.")
    parser.add_argument('--chrom', required=True, help='Chromosome number')

    args = parser.parse_args()
    main(chrom=args.chrom)
