library(dplyr)
library(data.table)
source("utils/r_options.r")

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--initial_sample_qc_file", required=TRUE, help="Path to INITIAL_SAMPLE_QC_FILE output from 03_0_initial_sample_qc.py")
parser$add_argument("--sample_information", required=TRUE,
    help=paste0("Path to sample information file (aka phenotype file) - this should contain two columns PCT_CHIMERAS and ",
        "PCT_CONTAMINATION, the chimeric read % and freemix contamination %, taken from the GATK/picard metadata. Also include any ",
        "factor you would like to split on and edit the commented code in this file to plot.")
parser$add_argument("--sample_list_initial_qc", required=TRUE, help="Path to initial sample list output .tsv file")
parser$add_argument("--sample_summary_count", required=TRUE, help="Path to initial sample filter summary counts .tsv file")
args <- parser$parse_args()

# Run the plotting again to ensure that the thresholds are as in the plots.
source("03_initial_1_sample_qc_plot.r")

SAMPLE_LIST_INITIAL_QC <- args$sample_list_initial_qc
SAMPLE_SUMMARY_COUNT <- args$sample_summary_count

df_out <- filter(df, call_rate > T_sample_callRate) %>%
	filter(PCT_CONTAMINATION < T_pct_contamination) %>%
	filter(PCT_CHIMERAS < T_pct_chimeras) %>%
	filter(dp_stats.mean > T_dpMean) %>%
	filter(gq_stats.mean > T_gqMean)

df_out <- df_out %>% select(s)
print(dim(df_out))

fwrite(df_out, file=SAMPLE_LIST_INITIAL_QC, quote=FALSE, row.names=FALSE, col.names=FALSE)

# Create the table too
df_summary_count <- data.table(
	"Filter" = c("Samples after initial filter",
			   paste0("Sample call rate < ", T_sample_callRate),
			   paste0("% FREEMIX contamination > ", T_pct_contamination),
			   paste0("% chimeric reads > ", T_pct_chimeras),
			   paste0("Mean DP < ", T_dpMean),
			   paste0("Mean GQ < ", T_gqMean),
			   "Samples after sample QC filters"),
	"Samples" = c(nrow(df),
			    nrow(filter(df, call_rate <= T_sample_callRate)),
				nrow(filter(df, PCT_CONTAMINATION >= T_pct_contamination)),
				nrow(filter(df, PCT_CHIMERAS >= T_pct_chimeras)),
				nrow(filter(df, dp_stats.mean <= T_dpMean)),
				nrow(filter(df, gq_stats.mean <= T_gqMean)),
				nrow(df_out))
	)

fwrite(df_summary_count, file=SAMPLE_SUMMARY_COUNT, quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
