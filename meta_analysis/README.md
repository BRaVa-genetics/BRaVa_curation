Summary of the figures dealt with in this readme

Figure 1 create_Figure4_count_of_assocs_plots.r, create_Figure4_meta_results_plots.r, create_Figure4_summary_information.r
Figure 2: create_Figure5_counts_of_ancestry_specific_assocs_plots.r
Figure 3: plot_meta_variant_rvas.r
Figure 4: plot_meta_variant_rvas.r, create_Figure4_meta_results_plots.r, plotting_results_not_found_in_biobanks.r, extract_burden_effects.r

Figure EDF1: NA
Figure EDF2: create_Figure2.r
Figure EDF3: NA

Figure S1: create_Figure2.r
Figure S2: create_Figure2.r
Figure S3: compare_prevalences_between_population_and_hospital_plots.r
Figure S4: create_inflation_figure.r
Figure S5: compare_effect_sizes_across_ancestries_plots.r
Figure S6: comparing_Neff_independent_plot.r
Figure S7: create_inflation_figure.r
Figure S8: create_Figure4_count_of_assocs_plots.r
Figure S9: plotting_results_not_found_in_biobanks.r
Figure S10: plotting_results_not_found_in_biobanks.r
Figure S11: plotting_results_not_found_in_biobanks.r
Figure S12: plotting_results_not_found_in_biobanks.r
Figure S13: create_Figure4_count_of_assocs_plots.r
Figure S14: create_Figure4_count_of_assocs_plots.r
Figure S15: create_Figure5_counts_of_ancestry_specific_assocs_plots.r
Figure S16: plot_meta_variant_rvas.r
Figure S17: plot_meta_variant_rvas.r
Figure S18: plot_meta_variant_rvas.r
Figure S19: plot_meta_variant_rvas.r
Figure S20: plot_meta_variant_rvas.r
Figure S21: plot_meta_variant_rvas.r
Figure S22: plot_meta_variant_rvas.r
Figure S23: meta_meta_analyses_comparisons.r (note, this is currently in the Jurgens_meta folder)
Figure S24: meta_meta_analyses_comparisons.r (note, this is currently in the Jurgens_meta folder)
Figure S25: meta_meta_analyses_comparisons.r (note, this is currently in the Jurgens_meta folder)
Figure S26: comparing_Neff_independent_plot.r
Figure S27: plot_before_vs_after.r
Figure S28: plot_before_vs_after.r
Figure S29: plot_before_vs_after.r

Run the munging for each biobank that has uploaded results
upload:
rsync -avz --progress --relative -- */cleaned/gene/* qen698@cluster2.bmrc.ox.ac.uk:/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/
Don't upload the extra_cauchy results files - these should go into a different folder

Submit the gene QQ plot jobs
Rscript run_gene_analysis_qq_gcloud_bmrc.r --analysis_results_folder /well/lindgren/dpalmer/BRaVa_meta-analysis_inputs --out_dir /well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/plots_biobank_specific

Determine genomic inflation of all of the tests
Rscript extract_genomic_inflation.r (keep defaults)
Determines metrics of inflation of test statistics before running meta-analysis
This needs to be run and then filtered on prior to meta-analysis. It will flag things such as the issue with MatHem, and VE in men from CCPM.
Note that at this point, we place the v7 version of BMI for AMR in place of the v8 version. However, v7 is never run through the pipeline due to the inflation in v8 which is used for filtering. During review, we can consider rerunning v8 or running v7 in it's place to boost power for BMI.
If we do the latter (v7 in place), we need to perform the following:
Importantly, delete the v8 version of BMI for AMR after this function is run. We then allow all-of-us BMI to not be excluded in the meta-analysis in the gene-based tests by hard-coding.

Download the results and take a look
Download the Neff weights files where possible
Combine them in the Neff weights file that is passed to the meta-analysis script
Upload the new meta-analysis script
Upload the Neff weights file

submit the meta-analysis scripts at the gene-level
Rscript run_gene_meta_analysis_gcloud_bmrc.r (note that Neff file is hard-coded as default in the slurm call - make sure it is correct)
check the meta-analysis QQ plotting functions
submit the meta-analysis QQ plotting functions
Rscript run_gene_meta_analysis_qq_gcloud_bmrc.r

Rscript run_gene_meta_analysis_gcloud_no_Neff_file_bmrc.r (to compare Neff)

Determine the genomic inflation after the meta-analysis to check
Rscript extract_genomic_inflation_meta.r
Rscript extract_burden_effects.r contains the forest plot code - this is used in Figure 4 d, e, f

Create the supplementary tables
Rscript create_supplementary_table_4_7.r
Rscript create_supplementary_table_8.r
Rscript compare_prevalences_between_population_and_hospital.r
Rscript compare_prevalences_between_population_and_hospital_plots.r

I created a short script to loop over the meta-analysed sumstats files (across all of them, and restrict to the pairs below a defined p-value cutoff)
Rscript generate_candidate_file_for_conditioning.r
Note that this reads a file that contains the significant associations in the full meta-analysis. This table is generated from a script within the manuscript_figures folder:
create_Figure4_count_of_assocs.r - so run this first

You now have an updated candidate list - upload to the bucket!

You can now create the collection of (gene, phenotype) pairs to pass to the agent.
Use genes_for_agent.r. Note that this is somewhat manual, as we're iterating, manually including gene-bass information etc, before downloading, so work through the file, reading what it says.

The supp tables containing all the gene-level results (long and wide) are all generated in this script (genes_for_agent.r).

plotting_results_not_found_in_biobanks.r
requires: manuscript_figures/biobank_results_at_sig_genes.tsv.gz, manuscript_tables/Tables/novelty_agent_results.tsv, manuscript_tables/Tables/gene_phenotype_pairs_for_agent_superset.tsv.gz, manuscript_tables/Tables/gene_phenotype_pairs_for_agent.tsv.gz, all of which are generated in genes_for_agent.r.

Next, create all the figures.

Figure EDF1: N/A

Figure EDF2
Rscript create_Figure2.r (note this is unlikely to change - it's the PC plots)
creates PCs and ternary plot: Figure2.pdf, FigureS1.pdf, FigureS2.pdf.

Figure EDF3: N/A

Figure 1
Rscript create_Figure4_count_of_assocs.r (reminder to run this first)
Rscript create_Figure4_count_of_assocs_plots.r
creates: cts_unique_hits_count.pdf = 1d, binary_unique_hits_count.pdf = 1e
Rscript create_Figure4_meta_results.r
Rscript create_Figure4_meta_results_plots.r
creates: `{meta_analysis, AFR, AMR, EAS, EUR, SAS, non_EUR}_unique_{cts, binary}_categories.{pdf,png}` and `{meta_analysis, AFR, AMR, EAS, EUR, SAS}_unique_{cts, binary}.{pdf,png}` (meta_analysis are 1f, 1g, 1h and 1i)
move 
`{AFR, AMR, EAS, EUR, SAS, non_EUR}_{cts, binary}_categories.png` to FigureS{3, 4, 5, 6, 7, 8}.png

Rscript create_Figure4_summary_information.r
creates: ancestry_summary.pdf = 1a and 1b, and phenotype_counts.pdf 1c.

These are combined manually in inkscape and stored at
`Figures_main_text/Figure4.{svg,pdf}`

Figure 2
Rscript create_Figure5_counts_of_ancestry_specific_assocs.r
Rscript create_Figure5_counts_of_ancestry_specific_assocs_plots.r
creates: 'venn' diagrams - `euler_plot_{all_ancestries, class, group, superpops, ukb_aou}.pdf`
upset_plot.pdf and upset_plot_{class, group, superpops, ukb_aou}.pdf

These are combined manually in inkscape and stored at `Figures_main_text/Figure5.{svg, pdf}`

Figure 3
(see below for details)
Rscript variant_output_results_for_plotting.r
Rscript plot_meta_variant_rvas.r

Figure 4
Rscript plot_meta_variant_rvas.r - phenotype specific at variant level are used (AFib)
Rscript create_Figure4_meta_results_plots.r - phenotype specific at gene level are used (T2D, AFib)
Rscript extract_burden_effects.r - forest plots for examples in d, e, f.
Rscript plotting_results_not_found_in_biobanks.r - a. Use the non-height version (as height makes it too big of a pic - demote to supp)

Supplementary Figures

Rscript create_inflation_figure.r
creates: `{AFR,AMR,EAS,EUR,SAS}_inflation.png` and meta_inflation.png

Rscript compare_effect_sizes_across_ancestries.r
Rscript compare_effect_sizes_across_ancestries_deming.r
Rscript compare_effect_sizes_across_ancestries_plots.r

comparing_Neff_independent_plot.r
creates comparing_Neff_independent.tsv.gz and Neff_biobank_plot.tsv.gz
The former is then used to plot results assuming independent samples for combining across cohorts (to define the weights), vs using the Neff estimated from the GRM. We also plot stouffer vs inverse_variance_burden.

Next, the investigation of sample overlap
First - null_correlation_investigation.r - change this to act directly on the files on BMRC
plot_corr_matrix.r - creates the supp fig plots.
Also, run the meta-analysis sample overlap code, based on these results.

Also, need the comparison of association results allowing for the sample overlap
plot_overlap_check.r

And allowing for the CAF to contribute
plot_caf_stouffer_check.r

Variant meta-analysis

install_bcftools_plugins.sh - code to get bcftools +metal up and running on the BMRC cluster
Rscript run_prepare_files_for_bcftools_metal.r
prepare_files_for_bcftools_metal.r - this munges the SAIGE output files ready for input into bcftools +metal as summary statistic vcf files. I've added code to estimate AF for GEL
run_variant_meta_analysis_gcloud_bmrc.r - this loops all phenotypes and performs all required variant level meta-analysis (LOBO, ancestry specific etc)
run_meta_analysis_bcftools_metal_variant_gcloud_bmrc.sh - the submission script called for each meta-analysis in the run_variant_meta_analysis_gcloud_bmrc.r script

Plotting the variant meta-analysis results

Rscript manuscript_figures/variant_output_results_for_plotting.r
Note that the above creates manhattan_plots_null_results_downsampled.tsv.gz which I use in create_variant_candidates.r
It also creates manhattan_plots.tsv.gz - which does not downsample, but only contains the 'all' meta-analysis results.
And it creates manhattan_plots_up_to_0.01.tsv.gz - which also does not downsample, only contains the 'all' meta-analysis results, but with a MAF filter of 0.01, rather than 0.001.
Finally, it creates all_assocs_from_meta_051125.tsv.gz which contains all of the gene-based results for the damaging annotations in the 'all' meta (though this is useful, it's currently not used for anything).
Then, running plot_meta_variant_rvas.r creates the 'manhattan' style plots
create_variant_candidates.r takes the downsampled (null portion is downsampled) data and filters to the significant associations and creates a table.
Rscript plot_meta_variant_rvas.r
This creates the manhattan plots for the variant level data

Rscript manuscript_figures/variant_qq_plotting.r

sign_agreement_between_pops.r has code that checks concordance in sign between gene-burden associations.

plotting the CHIP genes
CHIP_genes_plot.r

create correlation of phenotypes, with covariates regressed out
find_correlations_phenotypes.r
creates a pdf of the correlation matrix with black outlines around boxes with abs(r) > 0.05, and bonf sig.

Create confusion matrices using the columns in supplementary tables (taken from Jeremy's AI agent)
create_confusion_matrices.r

statistics_for_the_paper.r - this script determines a collection of stats that we report in the paper, for the meta and biobank level results.
This script calls count_phen_loci.r and map_variants_to_loci.r

final_meta_analysis.r contains the meta-analysis of the conditional results, plot and numbers for the paper.

Create upset plots for leave one out (LOO) analysis
create_Figure5_counts_of_LOO_assocs.r
create_Figure5_counts_of_LOO_plots.r