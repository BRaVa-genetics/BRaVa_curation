#!/bin/Rscript
library(data.table)
library(dplyr)
library(argparse)

tests <- c("Burden", "SKAT", "SKAT-O")
source("meta_analysis_utils.r")
dict <- c(
  "other_missense"     = "other_missense_or_protein_altering",
  "damaging_missense"  = "damaging_missense_or_protein_altering"
)
groups_for_cauchy <- c("pLoF", "damaging_missense_or_protein_altering",
    "pLoF;damaging_missense_or_protein_altering")
max_MAFs_for_cauchy <- c(1e-4, 1e-3)

# library(googlesheets4)
# dt_cts <- read_sheet("https://docs.google.com/spreadsheets/d/1LMtyhOIWgyPfs_UGv2CA6yh8jnXV7gCjN_5zRSi7ORw/edit?gid=70010180#gid=70010180", sheet="Table S14")
# dt_binary <- read_sheet("https://docs.google.com/spreadsheets/d/1LMtyhOIWgyPfs_UGv2CA6yh8jnXV7gCjN_5zRSi7ORw/edit?gid=70010180#gid=70010180", sheet="Table S15")
# dt <- rbind(dt_cts, dt_binary, use_names=TRUE)
# dt <- dt %>% rename(Region = `Gene ID`, phenotype = `Phenotype ID`, max_MAF = `max MAF`, Group = Mask) %>% select(Region, phenotype)
# setDT(dt)
# dt <- unique(dt)
# dt <- dt %>% filter(Region != TRUE) # Remove blank lines
# setkeyv(dt, c("Region", "phenotype"))
# # Final meta-analysis script to take the conditioned results and run Stouffers and inverse-variance weighted meta-analysis.
# # biobank_sig <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/biobank_results_at_sig_genes.tsv.gz")
# biobank_sig <- fread("manuscript_figures/biobank_results_at_sig_genes.tsv.gz")
# setkeyv(biobank_sig, c("Region", "phenotype"))
# biobank_sig_specific <- merge(dt, biobank_sig)
# biobank_sig_specific <- biobank_sig_specific %>% 
#     rename(biobank = dataset, pop = ancestry, phenotypeID = phenotype)
# setkeyv(biobank_sig_specific, c("Region", "phenotypeID", "max_MAF", "Group", "pop", "biobank"))
# fwrite(biobank_sig_specific, file="manuscript_figures/biobank_results_at_sig_genes_specific.tsv.gz")

# THIS WILL BE UPDATED!
biobank_sig_specific <- "/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/biobank_results_at_sig_genes_specific.tsv.gz"
# Read in and determine sample sizes to merge in
# Let's just get all of the sample sizes etc for everything
file_information <- rbindlist(lapply(grep("cleaned/gene",
    dir("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks", recursive=TRUE), value=TRUE),
    extract_file_info), fill=TRUE) %>% rename(N=n, N_ctrl=n_controls, N_case=n_cases, pop=ancestry, phenotypeID=phenotype) %>%
    mutate(biobank = gsub(".*gene/", "", dataset),
        N=as.integer(N), N_case=as.integer(N_case), N_ctrl=as.integer(N_ctrl)) %>% 
    select(biobank, pop, phenotypeID, sex, N, N_case, N_ctrl)
setkeyv(file_information, c("biobank", "pop", "phenotypeID"))

# All the (gene, trait, ancestry) tuples that should have had conditional analysis attempted.
AFR_cond <- fread(cmd = "gsutil cat gs://brava-meta-pilot-analysis/gene_phenotype_pairs_with_conditioning_variants_251125_AFR.tsv")
AMR_cond <- fread(cmd = "gsutil cat gs://brava-meta-pilot-analysis/gene_phenotype_pairs_with_conditioning_variants_251125_AMR.tsv")
EAS_cond <- fread(cmd = "gsutil cat gs://brava-meta-pilot-analysis/gene_phenotype_pairs_with_conditioning_variants_251125_EAS.tsv")
EUR_cond <- fread(cmd = "gsutil cat gs://brava-meta-pilot-analysis/gene_phenotype_pairs_with_conditioning_variants_251125_EUR.tsv")
SAS_cond <- fread(cmd = "gsutil cat gs://brava-meta-pilot-analysis/gene_phenotype_pairs_with_conditioning_variants_251125_SAS.tsv")

cond <- rbind(AFR_cond, AMR_cond, EAS_cond, EUR_cond, SAS_cond, use.names=TRUE)
cond <- cond %>% rename(Region = Gene, phenotypeID = Trait, pop = ancestry, max_MAF = MAF_cutoff_for_conditioning_variants)
cond <- unique(cond %>% select(Region, phenotypeID, pop))
setkeyv(cond, c("Region", "phenotypeID", "pop"))

main <- function(args)
{
	# Can then also consider appropriate reweighting and accounting for sample overlap.
	gene_dir <- args$gene_results_data_dir
	variant_dir <- args$variant_results_data_dir
	phe <- args$phenotypeID

	# First, read in the gene-level results files
	gene_results_files <- grep("analysis_results.txt$",
		dir(gene_dir, recursive=TRUE, full.names=TRUE), value=TRUE)
	variant_results_files <- grep("analysis_results.txt.singleAssoc.txt$",
		dir(variant_dir, recursive=TRUE, full.names=TRUE), value=TRUE)

	read_and_annotate <- function(file, gene_dir) {
		dt <- fread(file)
		dt$biobank <- gsub(paste0(gene_dir, "/([A-Za-z0-9-]+)/.*"), "\\1", file)
		dt$pop <- gsub(paste0(gene_dir, "/[A-Za-z0-9-]+/([A-Z]{3})/.*"), "\\1", file)
        if ("trait" %in% names(dt)) {
            dt <- dt %>% rename(Trait = trait)
        }
		return(dt)
	}

    munge_dt_unconditioned <- function(path)
    {
        dt_unconditioned <- fread(path) %>% 
            rename(Pvalue_prev = Pvalue, Pvalue_SKAT_prev = Pvalue_SKAT,
                Pvalue_Burden_prev = Pvalue_Burden, BETA_Burden_prev = BETA_Burden,
                SE_Burden_prev = SE_Burden, MAC_prev = MAC, MAC_case_prev = MAC_case,
                MAC_control_prev = MAC_control, Number_rare_prev = Number_rare,
                Number_ultra_rare_prev = Number_ultra_rare)
        setkeyv(dt_unconditioned, c("biobank", "pop", "phenotypeID"))
        dt_unconditioned <- merge(dt_unconditioned, file_information %>% 
            rename(N_old=N, N_case_old=N_case, N_ctrl_old=N_ctrl))
        dt_unconditioned <- dt_unconditioned %>% 
            mutate(N_eff_old = ifelse(!is.na(N_old), N_old, 4/(1/N_case_old + 1/N_ctrl_old)))
        setkeyv(dt_unconditioned, c("Region", "phenotypeID", "max_MAF", "Group", "pop", "biobank"))
        return(dt_unconditioned)
    }

	dt_gene <- rbindlist(lapply(gene_results_files, read_and_annotate, gene_dir=gene_dir), fill=TRUE)
	dt_variant <- rbindlist(lapply(variant_results_files, read_and_annotate, gene_dir=variant_dir), fill=TRUE)

	# Merge in the Ns
	dt_sample_size <- unique(
		dt_variant %>% select(biobank, pop, Trait, N, N_case, N_ctrl)) %>% 
		rename(phenotypeID = Trait)
	dt_gene <- dt_gene %>% rename(phenotypeID = Trait)
	dt_sample_size <- dt_sample_size %>% mutate(N_eff = ifelse(!is.na(N), N, 4/(1/N_case + 1/N_ctrl)))
	setkeyv(dt_gene, c("biobank", "pop", "phenotypeID"))

	dt_conditioned <- merge(dt_sample_size, dt_gene)
    dt_conditioned[, Group := sapply(
        strsplit(Group, ";"),
        function(tokens) {
            mapped <- ifelse(tokens %in% names(dict),
                dict[tokens], tokens)
            paste(mapped, collapse = ";")
        }
        )
    ]

	if (is.null(phe)) {
		phes <- c(
			"AMD", "Asth", "AFib", "BenCervUterNeo", "BenIntNeo",
			"BreastCanc", "CervCanc", "COPD", "CRF", "ColonRectCanc",
			"CAD", "EFRMB", "FemInf", "Gout", "HF", "HTN", "IBD",
			"IFHern", "ILDSarc", "MatHem", "NonRheuValv", "Pancreat",
			"PeptUlcer", "PAD", "Psori", "RheumHeaDis", "RheumArth",
			"Stroke", "T2Diab", "Urolith", "VaricVeins", "VTE", "ALT",
			"AlcCons", "AST", "BMI", "HDLC", "Height", "LDLC",
			"TChol", "TG", "WHRBMI", "HipRep", "CRP"
		)
	} else {
		phes <- phe
	}

	dt_conditioned <- dt_conditioned %>% filter(phenotypeID %in% phes)

    # Now read in the biobank data - ensure that the correlation of the P-values in (ancestry, biobank)
    # tuples is extremely high (> 0.99), if not, remove those from the subsequent analysis.
    # Also, at this point - merge in the (gene, trait) associations that do not have any conditioning
    # variants to be used. This will also allow us to sanity check against the numbers in Figure 2.
    thres <- 0.99
    dt_unconditioned <- munge_dt_unconditioned(biobank_sig_specific)
    dt_conditioned <- dt_conditioned %>% filter(Group %in% c("pLoF", "damaging_missense_or_protein_altering",
        "pLoF;damaging_missense_or_protein_altering"))
    setkeyv(dt_conditioned, c("Region", "phenotypeID", "max_MAF", "Group", "pop", "biobank"))

    # From this, we can then determine the collection of (gene, phenotype, ancestry) tuples that should have been
    # pulled through step 2. We can get this by creating a subset of dt_results restricted to the collection of genes
    # that have conditioning variants.
    cond_gene_trait <- unique(cond %>% select(Region, phenotypeID))
    setkeyv(cond_gene_trait, c("Region", "phenotypeID"))
    dt_unconditioned <- merge(cond_gene_trait, dt_unconditioned)
    setkeyv(dt_unconditioned, c("Region", "phenotypeID", "max_MAF", "Group", "pop", "biobank"))
    dt <- merge(dt_unconditioned, dt_conditioned, all.x = TRUE)

    # Also, include all the significant gene-trait associations that did not
    # require any variants to condition on - these are the set of (gene, trait) pairs that 
    # are not in 'cond_gene_trait', but are present in biobank_sig_specific.
    dt_unconditioned <- munge_dt_unconditioned(biobank_sig_specific)
    dt_no_conditioning_required <- unique(dt_unconditioned %>% select(Region, phenotypeID))
    dt_no_conditioning_required_key <- setdiff(dt_no_conditioning_required, unique(dt %>% select(Region, phenotypeID)))
    setkeyv(dt_no_conditioning_required_key, c("Region", "phenotypeID"))
    dt_no_conditioning_required <- merge(dt_unconditioned, dt_no_conditioning_required_key)
    dt <- rbind(dt, dt_no_conditioning_required, use.names=TRUE, fill=TRUE)

    correlation_with_previous <- data.frame(dt %>% group_by(biobank, pop) %>% 
        summarise(
            cor_SKATO = cor(Pvalue, Pvalue_prev, use="pairwise.complete.obs"),
            cor_Burden = cor(Pvalue_Burden, Pvalue_Burden_prev, use="pairwise.complete.obs"),
            cor_SKAT = cor(Pvalue_SKAT, Pvalue_SKAT_prev, use="pairwise.complete.obs"), n=n()))
    
    dt_tmp <- unique(dt %>%
        filter((!is.na(N) | !is.na(N_ctrl)) & (!is.na(N_old) | !is.na(N_ctrl_old))) %>%
        select(N, N_ctrl, N_case, N_old, N_ctrl_old, N_case_old, phenotypeID, biobank, pop))
    
    # Check that everything matches
    data.frame(dt_tmp %>% group_by(biobank, pop) %>% summarise(
        cor(N, N_old, use="pairwise.complete.obs"),
        cor(N_case, N_case_old, use="pairwise.complete.obs"),
        cor(N_ctrl, N_ctrl_old, use="pairwise.complete.obs")))

    # Here, we must remove any files that have been deemed to be inflated.
    dt_inflation <- fread(args$inflation_file)
    dt_inflation <- unique(dt_inflation %>% filter(Group == "synonymous") %>% 
        filter(max_MAF != 0.01,
            lambda_value > 1.3,
            !(lambda_type %in% c("lambda_50_Burden", "lambda_50_SKAT", "lambda_50"))) %>% 
        select(phenotype, dataset, ancestry))
    # Manual curation, adding the following (biobank, trait) tuples containing spurious 
    # associations
    dt_inflation <- rbind(dt_inflation, data.table(
        phenotype = c("ColonRectCanc", "Height"),
        dataset = c("egcut", "mgbb"),
        ancestry = c("EUR", "AMR"))
    ) %>% rename(phenotypeID = phenotype, biobank = dataset, pop = ancestry)
    dt_inflation <- setdiff(dt_inflation, data.table(
    phenotypeID = c("Height"),
    biobank = c("uk-biobank"),
    pop = c("EUR")))

    setkeyv(dt_inflation, c("phenotypeID", "biobank", "pop"))
    setkeyv(dt,  c("phenotypeID", "biobank", "pop"))
    dt <- setdiff(dt, merge(dt_inflation, dt)) # Need to fill in with the N counts here.

    dt <- dt %>% mutate(
        Pvalue = ifelse(is.na(Pvalue), Pvalue_prev, Pvalue),
        Pvalue_Burden = ifelse(is.na(Pvalue_Burden), Pvalue_Burden_prev, Pvalue_Burden),
        Pvalue_SKAT = ifelse(is.na(Pvalue_SKAT), Pvalue_SKAT_prev, Pvalue_SKAT),
        BETA_Burden = ifelse(is.na(BETA_Burden), BETA_Burden_prev, BETA_Burden),
        SE_Burden = ifelse(is.na(SE_Burden), SE_Burden_prev, SE_Burden),
        MAC = ifelse(is.na(MAC), MAC_prev, MAC),
        MAC_case = ifelse(is.na(MAC_case), MAC_case_prev, MAC_case),
        MAC_control = ifelse(is.na(MAC_control), MAC_control_prev, MAC_control),
        Number_rare = ifelse(is.na(Number_rare), Number_rare_prev, Number_rare),
        Number_ultra_rare = ifelse(is.na(Number_ultra_rare), Number_ultra_rare_prev, Number_ultra_rare))
    dt <- dt %>% mutate(Pvalue_cond = ifelse(is.na(Pvalue_cond), Pvalue, Pvalue_cond),
        Pvalue_Burden_cond = ifelse(is.na(Pvalue_Burden_cond), Pvalue, Pvalue_Burden_cond),
        Pvalue_SKAT_cond = ifelse(is.na(Pvalue_SKAT_cond), Pvalue, Pvalue_SKAT_cond),
        BETA_Burden_cond = ifelse(is.na(BETA_Burden_cond), BETA_Burden, BETA_Burden_cond),
        SE_Burden_cond = ifelse(is.na(SE_Burden_cond), SE_Burden, SE_Burden_cond))

    # If the biobank is present in dt, then include the corresponding results for that biobank where
    # it is absent (these are (gene, phenotype, ancestry) tuples for which no conditioning is carried out.
    # Before doing this, we should ensure that there are no (gene, trait) pairs for which we haven't performed 
    # conditioning, but should have - unlikely, but I think could happen if snakemake runs out of memory.
    dt <- dt %>% mutate(
        N = ifelse(is.na(N), N_old, N),
        N_case = ifelse(is.na(N_case), N_case_old, N_case),
        N_ctrl = ifelse(is.na(N_ctrl), N_ctrl_old, N_ctrl),
        N_eff = ifelse(is.na(N_eff), N_eff_old, N_eff))

	dt <- dt %>% filter(
		ifelse(is.na(N_ctrl), 
			N >= as.integer(args$case_control_threshold),
			(N_ctrl >= as.integer(args$case_control_threshold) & 
			 N_case >= as.integer(args$case_control_threshold))
		))
	dt <- dt %>% filter(Group != "Cauchy")

	# Update the Neff where we have the weights
	if (!args$no_Neff) {
		dt_weights <- fread(args$Neff_weights_file) %>% 
			rename(phenotypeID = pheno,
				pop = ancestry,
				biobank = dataset) %>% select(phenotypeID, pop, biobank, nglmm)
		setkeyv(dt_weights, c("phenotypeID", "biobank", "pop"))
		dt <- merge(dt, dt_weights, all.x=TRUE)
		dt <- dt %>% 
			mutate(N_eff = ifelse(is.na(nglmm), N_eff, nglmm))
    }
    setDT(dt)
    
    # Ensure valid Pvalues are generated
    dt <- dt %>% filter(
        (!is.na(Pvalue)) & (!is.na(Pvalue_SKAT)) & (!is.na(Pvalue_Burden))
    )
    # If the Pvalue_cond are NA, then then there were no variants to condition on
    dt <- dt %>% mutate(
        Pvalue_cond = ifelse(is.na(Pvalue_cond), Pvalue, Pvalue_cond),
        Pvalue_SKAT_cond = ifelse(is.na(Pvalue_SKAT_cond), Pvalue_SKAT, Pvalue_SKAT_cond),
        BETA_Burden_cond = ifelse(is.na(Pvalue_Burden_cond), BETA_Burden, BETA_Burden_cond),
        SE_Burden_cond = ifelse(is.na(Pvalue_Burden_cond), SE_Burden, SE_Burden_cond),
        Pvalue_Burden_cond = ifelse(is.na(Pvalue_Burden_cond), Pvalue_Burden, Pvalue_Burden_cond),
    )

    # Nudge P-values that are close to exactly 1 (for Stouffer).
    dt <- dt %>% mutate(
        Pvalue = ifelse(Pvalue > 0.99, 0.99, Pvalue),
        Pvalue_SKAT = ifelse(Pvalue_SKAT > 0.99, 0.99, Pvalue_SKAT),
        Pvalue_Burden = ifelse(Pvalue_Burden > 0.99, 0.99, Pvalue_Burden),
        Pvalue_cond = ifelse(Pvalue_cond > 0.99, 0.99, Pvalue_cond),
        Pvalue_Burden_cond = ifelse(Pvalue_Burden_cond > 0.99, 0.99, Pvalue_cond),
        Pvalue_SKAT_cond = ifelse(Pvalue_SKAT_cond > 0.99, 0.99, Pvalue_cond)
    )

    # dt <- dt %>% mutate(Pvalue_cond = ifelse((Pvalue_cond < 1e-10 & Pvalue > 0.001), Pvalue, Pvalue_cond))
    # dt <- dt %>% mutate(Pvalue_Burden_cond = ifelse((Pvalue_Burden_cond < 1e-10 & Pvalue_Burden > 0.001), Pvalue_Burden, Pvalue_Burden_cond))
    # dt <- dt %>% mutate(Pvalue_SKAT_cond = ifelse((Pvalue_SKAT_cond < 1e-10 & Pvalue_SKAT > 0.001), Pvalue_SKAT, Pvalue_SKAT_cond))

    dt_meta_cond <- list()
    for (cond in c(FALSE, TRUE)) {
        dt_meta <- list()
        for (test in tests)
        {
            cat(paste0(test, "...\n"))
            dt_meta[[test]] <- list()

            Pvalue_col <- ifelse(test == "SKAT-O", "Pvalue", paste0("Pvalue_", test))
            Pvalue_col <- ifelse(cond, paste0(Pvalue_col, "_cond"), Pvalue_col)
            Pvalue_out <- ifelse(cond, "Pvalue_cond", "Pvalue")

            BETA_col <- ifelse(cond, "BETA_Burden_cond", "BETA_Burden")
            SE_col <- ifelse(cond, "SE_Burden_cond", "SE_Burden")

            dt_to_test <- dt %>% group_by(Region, Group, max_MAF, phenotypeID)

            # Weighted Fisher's meta-analysis of p-values
            dt_meta[[test]][["weighted Fisher"]] <- run_weighted_fisher(
                dt_to_test, "N_eff", Pvalue_col, Pvalue_out,
                two_tail = ifelse(test == "Burden", TRUE, FALSE),
                input_beta = ifelse(test == "Burden", BETA_col, NULL)) %>% 
            mutate(Stat = NA, type="Weighted Fisher")

            # Stouffer's Z - Make sure P-values match, Stat= weighted_Z_Burden_Stouffer
            dt_meta[[test]][["Stouffer"]] <- run_stouffer(
                dt_to_test, "N_eff", "Stat", Pvalue_col, Pvalue_out,
                two_tail = ifelse(test == "Burden", TRUE, FALSE),
                input_beta = ifelse(test == "Burden", BETA_col, NULL)) %>% 
            mutate(type="Stouffer")
            dt_meta[[test]][["Stouffer"]] <- data.table(
                dt_meta[[test]][["Stouffer"]],
                key = c("Region", "Group", "max_MAF", "phenotypeID"))

            if (test == "Burden") {
                # And also run the inverse-variance weighted meta-analysis
                dt_meta[[test]][["inverse_variance_weighted"]] <- run_inv_var(
                    dt %>% group_by(Region, Group, max_MAF, phenotypeID), BETA_col, SE_col,
                    BETA_col, SE_col, Pvalue_out) %>%
                mutate(type="Inverse variance weighted")
                dt_meta[[test]][["inverse_variance_weighted"]] <- data.table(
                    dt_meta[[test]][["inverse_variance_weighted"]],
                    key = c("Region", "Group", "max_MAF", "phenotypeID"))

                # And also evaluate the heterogeneity P-values for both inverse variance weighted
                # and Stouffers
                dt_meta[[test]][["Stouffer"]] <- merge(dt_meta[[test]][["Stouffer"]],
                    data.table(run_heterogeneity_test(
                        weights(dt, FALSE, se_name=SE_col, n_eff_name="N_eff") %>% 
                            group_by(Region, Group, max_MAF, phenotypeID),
                        input_beta=BETA_col, output_meta_beta=BETA_col),
                    key=c("Region", "Group", "max_MAF", "phenotypeID"))
                )
                dt_meta[[test]][["inverse_variance_weighted"]] <- merge(dt_meta[[test]][["inverse_variance_weighted"]],
                    data.table(run_heterogeneity_test(
                        weights(dt, TRUE, se_name=SE_col, n_eff_name="N_eff") %>% 
                            group_by(Region, Group, max_MAF, phenotypeID),
                        input_beta=BETA_col, output_meta_beta=ifelse(cond, "BETA_meta_cond", "BETA_meta")),
                    key=c("Region", "Group", "max_MAF", "phenotypeID"))
                )
            }
            dt_meta[[test]] <- rbindlist(dt_meta[[test]], use.names=TRUE, fill=TRUE) %>% mutate(class=test)
        }
        dt_meta_cond[[as.character(cond)]] <- rbindlist(dt_meta, fill=TRUE)
    }

    dt_meta_cond[["TRUE"]] <- dt_meta_cond[["TRUE"]] %>% rename(
        Stat_cond = Stat, Pvalue_het_cond = Pvalue_het, sum_weights_cond = sum_weights,
        chisq_het_cond = chisq_het, df_cond = df)
    setkeyv(dt_meta_cond[["TRUE"]], c("Region", "Group", "max_MAF", "phenotypeID", "type", "class"))
    setkeyv(dt_meta_cond[["FALSE"]], c("Region", "Group", "max_MAF", "phenotypeID", "type", "class"))

    dt_meta <- merge(dt_meta_cond[["TRUE"]], dt_meta_cond[["FALSE"]])

    # Now, let's compare these results for the unconditioned versions, against the originals
    cat("writing the meta-analysed and conditioned results...\n")
    fwrite(dt_meta, file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/conditioned_meta_analysis_results.tsv.gz", quote=FALSE, sep="\t")

    dt_meta <- dt_meta %>% filter(!is.na(Pvalue) & !is.na(Pvalue_cond)) %>% 
        filter((class == "Burden" & type == "Inverse variance weighted") |
            (class %in% c("SKAT", "SKAT-O") & type == "Stouffer")) %>%
        mutate(Pvalue = ifelse(Pvalue > 0.99, 0.99, Pvalue),
            Pvalue_cond = ifelse(Pvalue_cond > 0.99, 0.99, Pvalue_cond),
            weights = 1)
    meta_cauchy <- run_cauchy(dt_meta %>% group_by(Region, phenotypeID), "weights", "Cauchy_stat", "Pvalue", "Cauchy_Pvalue")
    meta_cauchy_cond <- run_cauchy(dt_meta %>% group_by(Region, phenotypeID), "weights","Cauchy_stat_cond", "Pvalue_cond", "Cauchy_Pvalue_cond")
    setDT(meta_cauchy)
    setDT(meta_cauchy_cond)

    setkeyv(meta_cauchy, c("Region", "phenotypeID"))
    setkeyv(meta_cauchy_cond, c("Region", "phenotypeID"))
    meta_cauchy <- merge(meta_cauchy, meta_cauchy_cond)

    cat("writing the Cauchy of the meta-analysed results across region, test, and group...\n")
    fwrite(meta_cauchy, file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/conditioned_meta_analysis_cauchy_results.tsv.gz", quote=FALSE, sep="\t")
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--gene_results_data_dir", default="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/step2_conditioning",
	required=FALSE, help="Location of the folder containing the input into the meta-analysis")
parser$add_argument("--variant_results_data_dir", default="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/step2_conditioning",
	required=FALSE, help="Location of the folder containing the input into the meta-analysis")
parser$add_argument("--case_control_threshold", default=100, required=FALSE,
	help="Minimum number of cases")
parser$add_argument("--out_dir", default="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs",
	required=FALSE, help="Output folder path")
parser$add_argument("--inflation_file", default="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/inflation_summaries.tsv.gz",
	required=FALSE, help="Inflation file")
parser$add_argument("--phenotypeID", required=FALSE, default=NULL,
	help="The phenotype ID to run meta-analysis on. Note: must match the naming in input folder.")
parser$add_argument("--no_Neff", required=FALSE, default=FALSE, action='store_true')
parser$add_argument("--Neff_weights_file", required=FALSE,
    default="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/Neff/Neff_weights_may25.tsv.gz",
    help="File to pass effective sample sizes")
args <- parser$parse_args()

main(args)

# Determine the associations that were significant
sig <- unique(dt_meta %>% filter(
    ((type == "Inverse variance weighted" & Pvalue < 6.7e-7) | 
     (type == "Stouffer" & Pvalue < 2.5e-7)),
    Group %in% c("pLoF", "damaging_missense_or_protein_altering","pLoF;damaging_missense_or_protein_altering" )) %>% select(Region, phenotypeID))

# Non-sig after conditioning
not_sig <- unique(dt_meta %>% filter(
    ((type == "Inverse variance weighted" & Pvalue < 6.7e-7 & Pvalue_cond >= 6.7e-7) | 
     (type == "Stouffer" & Pvalue < 2.5e-7 & Pvalue_cond >= 2.5e-7)),
    Group %in% c("pLoF", "damaging_missense_or_protein_altering","pLoF;damaging_missense_or_protein_altering" )) %>% select(Region, phenotypeID))

# Still sig after conditioning
still_sig <- unique(dt_meta %>% filter(
    ((type == "Inverse variance weighted" & Pvalue < 6.7e-7 & Pvalue_cond < 6.7e-7) | 
     (type == "Stouffer" & Pvalue < 2.5e-7 & Pvalue_cond < 2.5e-7)),
    Group %in% c("pLoF", "damaging_missense_or_protein_altering","pLoF;damaging_missense_or_protein_altering" )) %>% select(Region, phenotypeID))

setdiff(sig, still_sig)

library(googlesheets4)
dt_cond <- fread("conditioned_meta_analysis_results.tsv.gz")
dt_cts <- read_sheet("https://docs.google.com/spreadsheets/d/1LMtyhOIWgyPfs_UGv2CA6yh8jnXV7gCjN_5zRSi7ORw/edit?gid=70010180#gid=70010180", sheet="Table S14")
dt_binary <- read_sheet("https://docs.google.com/spreadsheets/d/1LMtyhOIWgyPfs_UGv2CA6yh8jnXV7gCjN_5zRSi7ORw/edit?gid=70010180#gid=70010180", sheet="Table S15")
dt <- rbind(dt_cts, dt_binary, use_names=TRUE) %>% filter(`meta analyzed` == "ALL") %>% rename(Region = `Gene ID`, phenotypeID = `Phenotype ID`, max_MAF = `max MAF`, Group = Mask, type = `Meta analysis method`)
setkeyv(dt_cond, c("Region", "Group", "max_MAF", "phenotypeID", "type", "class"))
setDT(dt)
setkeyv(dt, c("Region", "Group", "max_MAF", "phenotypeID", "type", "class"))

dt_merge <- merge(dt, dt_cond)
plot(-log10(dt_merge$`Pvalue.x`), -log10(dt_merge$`Pvalue.y`))




