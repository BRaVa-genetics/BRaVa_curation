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

def main(chrom, pop):

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

	# Inputs
	MT = get_final_filter_mt_path(chrom)
	IMPUTESEX_TABLE = 'file:///mnt/project/Barney/qc2/06_impute_sex/04_imputesex_'
	SUPERPOPS = "file:///mnt/project/Barney/qc2/05_estimate_superpopulation/superpopulation_labels.tsv"
	SEXCHECK_LIST = 'file:///mnt/project/Barney/qc2/06_impute_sex/06_sexcheck.remove.BRaVa.sample_list'
	RELATED_SAMPLES = 'file:///mnt/project/Barney/qc2/07_king/07_king.related.sample_list'
	INITIAL_VARIANT_LIST = f'file:///mnt/project/Barney/qc2/02_prefilter_variants/02_variant_qc_chr{chrom}_list.tsv'
	SAMPLE_LIST_INITIAL_QC = 'file:///mnt/project/Barney/qc/03_initial_sample/03_initial_qc.keep.sample_list'

	REFERENCE = 'GRCh38'

	# Outputs
	VARIANT_QC_FILE = f'08_final_qc.variants_chr{chrom}_'

	print("Inputs:")
	print('MT; input matrix table: ', MT)
	print('IMPUTESEX_TABLE; input file prefix of .tsvs file to plot imputed sex information output from 06_0_impute_sex_superpopulations.py: ', IMPUTESEX_TABLE)
	print('SUPERPOPS; list of 1000G labels to loop over:', SUPERPOPS)
	print('SEXCHECK_LIST; list of sex swaps to remove:', SEXCHECK_LIST)
	print('RELATED_SAMPLES; input set of related samples output from 07_0_ibd.py, 07_0_pc_relate.py, or 07_0_ukb_relatedness_king.py', RELATED_SAMPLES)
	print('INITIAL_VARIANT_LIST; set of variants that pass initial filtering \
		(output from 02_prefilter_variants.py): ', INITIAL_VARIANT_LIST)
	print('SAMPLE_LIST_INITIAL_QC; set of initial samples output from 03_01_initial_sample_qc_filter.r: ', SAMPLE_LIST_INITIAL_QC)

	print("Outputs:")
	print('VARIANT_QC_FILE; output .tsv file variant QC information: ', VARIANT_QC_FILE)

	ht_superpops = hl.import_table(SUPERPOPS, impute=True).key_by("sample.ID").select("classification_strict")
	ht_related_samples = hl.import_table(RELATED_SAMPLES, no_header=True, key='f0')
	ht_initial_variants = hl.import_table(INITIAL_VARIANT_LIST, types={"locus": hl.tstr, "alleles":hl.dtype('array<str>')})
	ht_initial_variants = ht_initial_variants.annotate(locus=hl.expr.functions.parse_locus(ht_initial_variants.locus, reference_genome='GRCh38'))
	ht_initial_variants = ht_initial_variants.key_by('locus', 'alleles')

	ht_sexcheck_samples = hl.import_table(SEXCHECK_LIST, no_header=True, key='f0')
	ht_initial_samples = hl.import_table(SAMPLE_LIST_INITIAL_QC, no_header=True, key='f0')

	mt = hl.read_matrix_table(MT)
	mt = mt.annotate_cols(**ht_superpops[mt.s])
	mt = mt.filter_rows(hl.is_defined(ht_initial_variants[mt.row_key]))
	mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
	mt = mt.filter_cols(~hl.is_defined(ht_sexcheck_samples[mt.col_key]))
	mt = mt.filter_cols(~hl.is_defined(ht_related_samples[mt.col_key]))

	impute_sex_annotations = hl.import_table(IMPUTESEX_TABLE + pop + '.tsv.bgz', types={'s': hl.tstr, 'impute_sex.is_female': hl.tbool, 'impute_sex.f_stat': hl.tfloat, 'impute_sex.n_called': hl.tint, 'impute_sex.expected_homs': hl.tfloat, 'impute_sex.observed_homs': hl.tint})
	impute_sex_annotations = impute_sex_annotations.key_by('s')
	
	mt = mt.annotate_cols(filter_pop = (mt.classification_strict == pop))
	mt = mt.filter_cols(mt.filter_pop == True)
	mt = mt.annotate_cols(imputesex = impute_sex_annotations[mt.col_key])
	mt = hl.variant_qc(mt, name='variant_qc')
	mt = mt.annotate_rows(
		variant_qc=mt.variant_qc.annotate(AC=mt.variant_qc.AC[1],
		AF=mt.variant_qc.AF[1],
		homozygote_count=mt.variant_qc.homozygote_count[1]))
		
	mt = mt.annotate_rows(
		variant_qc=mt.variant_qc.annotate(
			p_value_hwe=hl.case()
			.when(mt.locus.in_autosome(), mt.variant_qc.p_value_hwe)
			.default(hl.agg.filter(mt.imputesex['impute_sex.is_female'],
								hl.agg.hardy_weinberg_test(mt.GT).p_value)),
			het_freq_hwe=hl.case()
			.when(mt.locus.in_autosome(), mt.variant_qc.het_freq_hwe)
			.default(hl.agg.filter(mt.imputesex['impute_sex.is_female'],
								hl.agg.hardy_weinberg_test(mt.GT).het_freq_hwe))
		)
	)
	
	mt.rows().select('variant_qc').flatten().export(f'{VARIANT_QC_FILE}{pop}.tsv.bgz')
	
	os.system(f'hdfs dfs -get ./{VARIANT_QC_FILE}{pop}.tsv.bgz .')

	export_file(
		path=f'{VARIANT_QC_FILE}{pop}.tsv.bgz', 
		out_folder='/Barney/qc2/08_0_final_variant_qc/'
	)

if __name__=='__main__':

	parser = argparse.ArgumentParser(description="Process some integers.")
	parser.add_argument('--chrom', required=True, help='Chromosome number')
	parser.add_argument("--pop", type=str, required=True, help="Population to process (e.g., AFR, AMR, EAS, EUR, SAS)")

	args = parser.parse_args()
	main(chrom=args.chrom, pop=args.pop)
