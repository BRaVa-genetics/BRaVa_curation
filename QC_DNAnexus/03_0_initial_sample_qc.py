import hail as hl
import sys
import pyspark
import dxpy
import hail as hl
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
	hl.init(sc=sc)
     
	def get_final_filter_mt(chrom):
		INITIAL_QC_SAMPLES = 'file:///mnt/project/ukb_wes_450k_qc/data/ukb_wes_450k.qced.sample_ids.tsv'
		ht_initial_samples = hl.import_table(INITIAL_QC_SAMPLES, no_header=True, key='f0')

		mt = hl.read_matrix_table(f'{database_dir}/04_final_filter_v2_write_to_mt/ukb_wes_450k.qced.chr{chrom}.mt')
		mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))

		return mt

	# Create sample QC metrics restricted and not restricted (target plus padding) 
	# to the target intervals.

	# Inputs
	INITIAL_VARIANT_LIST = f'file:///mnt/project/Barney/qc2/02_prefilter_variants/02_variant_qc_chr{chrom}_list.tsv'
	TARGET_INTERVALS = 'file:///mnt/project/resources/xgen_plus_spikein.b38.chr_prefix.bed'
	REFERENCE = 'GRCh38'

	# Outputs
	INITIAL_SAMPLE_QC_FILE  = f'03_initial_sample_qc_chr{chrom}.tsv'

	variants_to_filter = hl.import_table(INITIAL_VARIANT_LIST,
		types={'locus':hl.tlocus(reference_genome=REFERENCE), 'alleles':hl.tarray(hl.tstr)})
	variants_to_filter = variants_to_filter.key_by(locus=variants_to_filter.locus, alleles=variants_to_filter.alleles)

	mt = get_final_filter_mt(chrom)
	mt = mt.filter_rows(hl.is_defined(variants_to_filter[mt.row_key]))
	mt = mt.annotate_cols(gq = hl.agg.stats(mt.GQ), dp = hl.agg.stats(mt.DP))
	mt = hl.sample_qc(mt, name='qc_padded_target')

	# Import the target interval lists.
	target_intervals = hl.import_locus_intervals(TARGET_INTERVALS, reference_genome=REFERENCE)
	mt = mt.annotate_rows(not_in_target_intervals = ~hl.is_defined(target_intervals[mt.locus]))
	mt = mt.filter_rows(mt.not_in_target_intervals, keep=False)
	mt = hl.sample_qc(mt, name='qc_target')

	mt.cols().select('qc_padded_target', 'qc_target', 'gq', 'dp').flatten().export(output=INITIAL_SAMPLE_QC_FILE)

	os.system('hdfs dfs -get ./*.tsv .')

	export_file(
		path=INITIAL_SAMPLE_QC_FILE, 
		out_folder='/Barney/qc2/03_initial_sample/'
	)

if __name__=='__main__':

    parser = argparse.ArgumentParser(description="Process some integers.")
    parser.add_argument('--chrom', required=True, help='Chromosome number')

    args = parser.parse_args()
    main(chrom=args.chrom)
