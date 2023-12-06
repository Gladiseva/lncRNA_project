## Installation and Loading Packages
# Function to check and install a package
install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

# List of packages
required_packages <- c(
  "BiocManager", "rhdf5", "sleuth", "biomaRt",
  "devtools", "ggplot2", "ggrepel", "readxl",
  "svglite", "readr", "dplyr", "pheatmap"
)

# Install and load packages
invisible(sapply(required_packages, install_if_missing))

# Load libraries
library(readr)
library(dplyr)
library(rhdf5)
library(devtools)
library(sleuth)
library(BiocManager)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(readxl)
library(svglite)
library(pheatmap)

set.seed(123)

### Get results from SLEUTH
# specify where the kallisto results are stored.
sample_id <- dir(file.path(".", "results"))

# paths to the kallisto results indexed by the sample IDs is collated with
kallisto_dirs <- file.path(".", "results", sample_id)

# load file that describes the experimental design
sample2condition <- read.table(file.path(".", "NSCLC_exp_design.txt"),
                               header = TRUE,
                               stringsAsFactors = FALSE)

# directories must be appended in a new column
sample2condition <- dplyr::mutate(sample2condition, path = kallisto_dirs)

# Important: before run prepare_annotation_table.R
annotation_data <- read.csv("annotation_table_from_ALLandENS.csv")

# transformation_function because by default the transformation of counts is natural log
so_transcript_level <- sleuth_prep(sample2condition,
                             target_mapping = annotation_data,
                             read_bootstrap_tpm = TRUE,
                             extra_bootstrap_summary = TRUE,
                             transformation_function = function(x) log2(x + 0.5))
so_transcript_level <- sleuth_fit(so_transcript_level, ~condition, "full")
so_transcript_level <- sleuth_fit(so_transcript_level, ~1, "reduced")
# OUTPUT: NA values were found during variance shrinkage estimation LOESS
# Transcript level: These are the target ids with NA values: ENST00000361624.2, ENST00000387347.2

#likelihood ratio test (LRT)
so_transcript_level <- sleuth_lrt(so_transcript_level, "reduced", "full")
models(so_transcript_level)

# Test significant differences between conditions using the Wald test
# Wald test for differential expression of isoforms. var oe -> observed2expected
oe_transcript_level <- sleuth_wt(so_transcript_level, which_beta = "conditionparaclonal")
sleuth_results_oe_tl <- sleuth_results(oe_transcript_level,
                                    test = "conditionparaclonal",
                                    show_all = TRUE)

# top 20 significant genes with a q-value <= 0.05.
sleuth_significant_tl <- dplyr::filter(sleuth_results_oe_tl, qval <= 0.05)
sleuth_sihnificant_top_20_tl <- head(sleuth_significant_tl, 20)

# save outputs to csv
write.csv(sleuth_results_oe_tl, file = "sleuth_results_full_tl.csv", row.names = FALSE)
write.csv(sleuth_significant_tl, file = "sleuth_results_significant_tl.csv", row.names = FALSE)

# Print or further process the unique gene/transcript biotypes
unique_gene_biotypes <- unique(sleuth_results_oe_tl$gene_biotype)
unique_transcript_biotypes <- unique(sleuth_results_oe_tl$transcript_biotype)

# lncRNA
lncRNA_result_tl <- subset(sleuth_results_oe_tl, gene_biotype == "lncRNA" | transcript_biotype == "lncRNA")
lncRNA_result_top_20_tl <- head(lncRNA_result_tl, 20)
write.csv(lncRNA_result_tl, file = "sleuth_results_lncRNA_tl.csv", row.names = FALSE)

# protein_coding
protein_coding_result_tl <- subset(sleuth_results_oe_tl, gene_biotype == "protein_coding" | transcript_biotype == "protein_coding")
protein_coding_result_top_20_tl <- head(protein_coding_result_tl, 20)
write.csv(protein_coding_result_tl, file = "sleuth_results_protein_coding_tl.csv", row.names = FALSE)

# NAs
NA_result_tl <- sleuth_results_oe_tl[grepl("MSTRG", sleuth_results_oe_tl$target_id), ]
NA_result_top_20_tl <- head(NA_result_tl, 20)
write.csv(NA_result_tl, file = "sleuth_results_NAs_tl.csv", row.names = FALSE)

#normalized expression values and estimated counts as the expression units
counts_per_replicate_tl <- sleuth_to_matrix(oe_transcript_level, "obs_norm", "est_counts")
write.csv(counts_per_replicate_tl, file = "counts_per_replicate_tl.csv", row.names = TRUE)

so_gene_level <- sleuth_prep(sample2condition,
                             target_mapping = annotation_data,
                             read_bootstrap_tpm = TRUE,
                             extra_bootstrap_summary = TRUE,
                             gene_mode = TRUE, # diferent from tl
                             aggregation_column = 'ens_gene', # diferent from tl
                             transformation_function = function(x) log2(x + 0.5))
so_gene_level <- sleuth_fit(so_gene_level, ~condition, "full")
so_gene_level <- sleuth_fit(so_gene_level, ~1, "reduced")
# OUTPUT: NA values were found during variance shrinkage estimation LOESS
# Gene level: NA values: ENSG00000004848, ENSG00000140873, ENSG00000206052, ENSG00000272894
# NA values: ENSG00000140873, ENSG00000206052, ENSG00000272894, ENSG00000198804

#likelihood ratio test (LRT)
so_gene_level <- sleuth_lrt(so_gene_level, "reduced", "full")
models(so_gene_level)

# Test significant differences between conditions using the Wald test
oe_gene_level <- sleuth_wt(so_gene_level, which_beta = "conditionparaclonal")
sleuth_results_oe_gl <- sleuth_results(oe_gene_level,
                                    test = "conditionparaclonal",
                                    show_all = TRUE)

# top 20 significant genes with a q-value <= 0.05.
sleuth_significant_gl <- dplyr::filter(sleuth_results_oe_gl, qval <= 0.05)
sleuth_sihnificant_top_20_gl <- head(sleuth_significant_gl, 20)

# save outputs to csv
write.csv(sleuth_results_oe_gl, file = "sleuth_results_full_gl.csv", row.names = FALSE)
write.csv(sleuth_significant_gl, file = "sleuth_results_significant_gl.csv", row.names = FALSE)

#normalized expression values and estimated counts as the expression units
counts_per_replicate <- sleuth_to_matrix(oe_gene_level, "obs_norm", "scaled_reads_per_base")
write.csv(counts_per_replicate, file = "counts_per_replicate_gl.csv", row.names = TRUE)

### Exploratory analysis
# Create a volcano plot
create_volcano_plot <- function(data, title = "Volcano Plot") {
  ggplot(data, aes(x = b, y = -log10(pval))) +
    geom_point(aes(color = qval < 0.05), alpha = 0.5, size = 2) +
    geom_text_repel(aes(label = ifelse(!is.na(get("ext_gene")), get("ext_gene"), get("ens_gene"))), 
                    data = subset(data, qval < 0.05), 
                    box.padding = 0.5, point.padding = 0.2, segment.color = "grey50") +
    scale_color_manual(values = c("blue", "red")) +
    labs(title = title, x = "Log2-Fold Change (b)", y = "-log10(p-value)") +
    theme_minimal()
}

# plot for all values
# create_volcano_plot(sleuth_results_oe_tl, title = "All transcripts")

# lncRNA_result (ext_gene as label)
create_volcano_plot(lncRNA_result_tl, title = "lncRNAs, transcripts")

# protein_coding_result (ext_gene as label)
create_volcano_plot(protein_coding_result_tl,
                    title = "Protein coding genes, transcripts")

# NA_result (target_id as label)
create_volcano_plot(NA_result_tl,
                    title = "Novel genes, transcripts")

# exploratory analysis
# interactive visualization
sleuth_live(oe_transcript_level)

# PCA to visualize sample relationships based on gene expression
plot_pca(so_transcript_level, color_by = 'condition')
# Group density plot to visualize the distribution of gene expression levels across conditions
plot_group_density(so_transcript_level,
                   use_filtered = TRUE,
                   units = "est_counts",
                   trans = "log",
                   grouping = setdiff(colnames(so_transcript_level$sample_to_covariates),
                                      "sample"), offset = 1)

# check if found genes are present in the paper
df_one <- read.csv("sleuth_results_significant_tl.csv")
df_two <- read.csv("parental_vs_para_from_paper.csv")

unique_ens_gene <- unique(df_one$ens_gene)
unique_ensembl_gene_id <- unique(df_two$ensembl_gene_id)

cat("Unique values in df_one$ens_gene:", length(unique_ens_gene), "\n")
cat("Unique values in df_two$ensembl_gene_id:", length(unique_ensembl_gene_id), "\n")

# Display values that are in df_one but not in df_two
ens_gene_not_in_df_two <- setdiff(unique_ens_gene, unique_ensembl_gene_id)
cat("Values in df_one$ens_gene not found in df_two$ensembl_gene_id:", length(ens_gene_not_in_df_two), "\n")
ens_gene_not_in_df_two
# Remove elements containing "MSTRG"
filtered_genes <- ens_gene_not_in_df_two[!grepl("MSTRG", ens_gene_not_in_df_two)]

# Print the result
print(filtered_genes)
length(filtered_genes)


# Read the data from the CSV file (adjust the file path accordingly)
## for all significant top 20
data <- read.csv("counts_per_replicate_tl.csv", row.names = 1)
data_filtered <- data[rownames(data) %in% sleuth_sihnificant_top_20_tl$target_id, ]

data_matrix <- as.matrix(data_filtered)

# Create a heatmap using pheatmap
pheatmap(data_matrix, 
         cluster_cols = FALSE, # Turn off column clustering
         scale = "row",        # Scale rows (genes) by default
         main = "Significant transcripts top 20",
         annotation_col = NULL)  # Turn off column annotations

## for lncRNAs top 20
data <- read.csv("counts_per_replicate_tl.csv", row.names = 1)
data_filtered <- data[rownames(data) %in% lncRNA_result_top_20_tl$target_id, ]

data_matrix <- as.matrix(data_filtered)

# Create a heatmap using pheatmap
pheatmap(data_matrix, 
         cluster_cols = FALSE, # Turn off column clustering
         scale = "row",        # Scale rows (genes) by default
         main = "lncRNAs top 20",
         annotation_col = NULL)  # Turn off column annotations

## for protein coding top 20
data <- read.csv("counts_per_replicate_tl.csv", row.names = 1)
data_filtered <- data[rownames(data) %in% protein_coding_result_top_20_tl$target_id, ]

data_matrix <- as.matrix(data_filtered)

# Create a heatmap using pheatmap
pheatmap(data_matrix, 
         cluster_cols = FALSE, # Turn off column clustering
         scale = "row",        # Scale rows (genes) by default
         main = "Protein coding, top 20",
         annotation_col = NULL)  # Turn off column annotations

## for NA_result_top_20
data <- read.csv("counts_per_replicate_tl.csv", row.names = 1)
data_filtered <- data[rownames(data) %in% NA_result_top_20_tl$target_id, ]

data_matrix <- as.matrix(data_filtered)

# Create a heatmap using pheatmap
pheatmap(data_matrix, 
         cluster_cols = FALSE, # Turn off column clustering
         scale = "row",        # Scale rows (genes) by default
         main = "Novel transcripts, top 20",
         annotation_col = NULL)  # Turn off column annotations

# plot bootstrap , to account for kallisto estimation ENST00000257555.11
target_id <- "ENST00000257555.11"  # Replace with the actual transcript or gene identifier
units <- "est_counts"
color_by <- "condition"  # Replace with the actual grouping variable in your sample metadata
x_axis_angle <- 50
divide_groups <- TRUE

# Call the plot_bootstrap function
plot_bootstrap(oe_transcript_level, target_id, units = units, color_by = color_by,
               x_axis_angle = x_axis_angle, divide_groups = divide_groups)
