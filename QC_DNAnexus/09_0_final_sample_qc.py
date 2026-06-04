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

    # Inputs:
    MT_HARDCALLS = get_final_filter_mt_path(chrom)
    IMPUTESEX_TABLE = 'file:///mnt/project/Barney/qc2/06_impute_sex/04_imputesex_'
    SUPERPOPS = "file:///mnt/project/Barney/qc2/05_estimate_superpopulation/superpopulation_labels.tsv"
    SEXCHECK_LIST = 'file:///mnt/project/Barney/qc2/06_impute_sex/06_sexcheck.remove.BRaVa.sample_list'
    SAMPLE_LIST_INITIAL_QC = 'file:///mnt/project/Barney/qc2/03_initial_sample/03_initial_qc.keep.sample_list'
    INITIAL_VARIANT_LIST = f'file:///mnt/project/Barney/qc2/02_prefilter_variants/02_variant_qc_chr{chrom}_list.tsv'
    FINAL_VARIANT_LIST = 'file:///mnt/project/Barney/qc2/08_0_final_variant_qc/08_final_qc.pop.keep.variant_list'
    REFERENCE = 'GRCh38'

    # Outputs:
    SAMPLE_BEFORE_QC_FILE = f'09_final_qc_chr{chrom}.before'
    SAMPLE_AFTER_QC_FILE = f'09_final_qc_chr{chrom}.after'

    ht_superpops = hl.import_table(SUPERPOPS, impute=True).key_by("sample.ID").select("classification_strict")
                                            
    ht_initial_variants = hl.import_table(INITIAL_VARIANT_LIST, types={"locus": hl.tstr, "alleles":hl.dtype('array<str>')})
    ht_initial_variants = ht_initial_variants.annotate(locus=hl.expr.functions.parse_locus(ht_initial_variants.locus, reference_genome='GRCh38'))
    ht_initial_variants = ht_initial_variants.key_by('locus', 'alleles')
                                            
    ht_initial_samples = hl.import_table(SAMPLE_LIST_INITIAL_QC, no_header=True, key='f0')
    ht_sexcheck_samples = hl.import_table(SEXCHECK_LIST, no_header=True, key='f0')

    ht_final_variants = hl.import_table(FINAL_VARIANT_LIST,
        types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
    ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)

    mt_before = hl.read_matrix_table(MT_HARDCALLS)
    mt_before = mt_before.annotate_cols(**ht_superpops[mt_before.s])
    mt_before = mt_before.filter_rows(hl.is_defined(ht_initial_variants[mt_before.row_key]))

    impute_sex_annotations = hl.import_table(IMPUTESEX_TABLE + pop + '.tsv.bgz', types={'s': hl.tstr, 'impute_sex.is_female': hl.tbool, 'impute_sex.f_stat': hl.tfloat, 'impute_sex.n_called': hl.tint, 'impute_sex.expected_homs': hl.tfloat, 'impute_sex.observed_homs': hl.tint})
    impute_sex_annotations = impute_sex_annotations.key_by('s')
    
    mt_before = mt_before.annotate_cols(filter_pop = (mt_before.classification_strict == pop))
    mt_before_pop = mt_before.filter_cols(mt_before.filter_pop == True)
    mt_before_pop = mt_before_pop.filter_cols(hl.is_defined(ht_initial_samples[mt_before_pop.col_key]))
    mt_before_pop = mt_before_pop.filter_cols(~hl.is_defined(ht_sexcheck_samples[mt_before_pop.col_key]))
    mt_before_pop = hl.variant_qc(mt_before_pop, name = 'variant_qc')
    mt_before_pop = mt_before_pop.annotate_rows(
        variant_qc = mt_before_pop.variant_qc.annotate(
            AC=mt_before_pop.variant_qc.AC[1],
            AF = mt_before_pop.variant_qc.AF[1],
            homozygote_count = mt_before_pop.variant_qc.homozygote_count[1]
            )
        )
    mt_before_pop = mt_before_pop.filter_rows(
        (mt_before_pop.variant_qc.AF > 0) & (mt_before_pop.variant_qc.AF < 1)
        )
    mt_before_pop = hl.sample_qc(mt_before_pop, name='sample_qc')
    ht_final_variants_pop = ht_final_variants.annotate(filter_pop = (ht_final_variants.pop == pop))
    ht_final_variants_pop = ht_final_variants_pop.filter(ht_final_variants_pop.filter_pop)
    mt_after_pop = mt_before_pop.filter_rows(hl.is_defined(ht_final_variants[mt_before_pop.row_key]))
    mt_after_pop = hl.sample_qc(mt_after_pop, name='sample_qc')
    mt_before_pop.cols().select("sample_qc").flatten().export(f'{SAMPLE_BEFORE_QC_FILE}.{pop}.samples.tsv.bgz')
    mt_after_pop.cols().select("sample_qc").flatten().export(f'{SAMPLE_AFTER_QC_FILE}.{pop}.samples.tsv.bgz')
    
    os.system(f'hdfs dfs -get {SAMPLE_BEFORE_QC_FILE}.{pop}.samples.tsv.bgz .')
    os.system(f'hdfs dfs -get {SAMPLE_AFTER_QC_FILE}.{pop}.samples.tsv.bgz .')

    export_file(
        path=f'{SAMPLE_BEFORE_QC_FILE}.{pop}.samples.tsv.bgz', 
        out_folder='/Barney/qc2/09_0_final_sample_qc/'
    )

    export_file(
        path=f'{SAMPLE_AFTER_QC_FILE}.{pop}.samples.tsv.bgz', 
        out_folder='/Barney/qc2/09_0_final_sample_qc/'
    )
                                            
if __name__=='__main__':

    parser = argparse.ArgumentParser(description="Process some integers.")
    parser.add_argument('--chrom', required=True, help='Chromosome number')
    parser.add_argument("--pop", type=str, required=True, help="Population to process (e.g., AFR, AMR, EAS, EUR, SAS)")

    args = parser.parse_args()
    main(chrom=args.chrom, pop=args.pop)
