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
sleuth_object <- sleuth_prep(sample2condition,
                             target_mapping = annotation_data,
                             read_bootstrap_tpm = TRUE,
                             extra_bootstrap_summary = TRUE,
                             gene_mode = TRUE, # remove if only trascript level needed
                             aggregation_column = 'ens_gene', # remove if only trascript level needed
                             transformation_function = function(x) log2(x + 0.5))
sleuth_object <- sleuth_fit(sleuth_object, ~condition, "full")
sleuth_object <- sleuth_fit(sleuth_object, ~1, "reduced")
# OUTPUT: NA values were found during variance shrinkage estimation LOESS
# Transcript level: These are the target ids with NA values: ENST00000361624.2, ENST00000387347.2
# Gene level: NA values: ENSG00000004848, ENSG00000140873, ENSG00000206052, ENSG00000272894
# NA values: ENSG00000140873, ENSG00000206052, ENSG00000272894, ENSG00000198804
#likelihood ratio test (LRT)
sleuth_object <- sleuth_lrt(sleuth_object, "reduced", "full")
models(sleuth_object)

# Test significant differences between conditions using the Wald test
# Wald test for differential expression of isoforms. var oe -> observed2expected
oe <- sleuth_wt(sleuth_object, which_beta = "conditionparaclonal")
sleuth_results_oe <- sleuth_results(oe,
                                    test = "conditionparaclonal",
                                    show_all = TRUE)

# top 20 significant genes with a q-value <= 0.05.
sleuth_significant <- dplyr::filter(sleuth_results_oe, qval <= 0.05)
sleuth_sihnificant_top_20 <- head(sleuth_significant, 20)

# save outputs to csv
write.csv(sleuth_results_oe, file = "sleuth_results_full.csv", row.names = FALSE)
write.csv(sleuth_significant, file = "sleuth_results_significant.csv", row.names = FALSE)

# Print or further process the unique gene/transcript biotypes
unique_gene_biotypes <- unique(sleuth_results_oe$gene_biotype)
unique_transcript_biotypes <- unique(sleuth_results_oe$transcript_biotype)

# lncRNA
lncRNA_result <- subset(sleuth_results_oe, gene_biotype == "lncRNA" | transcript_biotype == "lncRNA")
lncRNA_result_top_20 <- head(lncRNA_result, 20)
write.csv(lncRNA_result, file = "sleuth_results_lncRNA.csv", row.names = FALSE)

# protein_coding
protein_coding_result <- subset(sleuth_results_oe, gene_biotype == "protein_coding" | transcript_biotype == "protein_coding")
protein_coding_result_top_20 <- head(protein_coding_result, 20)
write.csv(protein_coding_result, file = "sleuth_results_protein_coding.csv", row.names = FALSE)

# NAs
NA_result <- sleuth_results_oe[grepl("MSTRG", sleuth_results_oe$target_id), ]
NA_result_top_20 <- head(NA_result, 20)
write.csv(NA_result, file = "sleuth_results_NAs.csv", row.names = FALSE)

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
create_volcano_plot(sleuth_results_oe, title = "All genes")

# lncRNA_result (ext_gene as label)
create_volcano_plot(lncRNA_result, title = "Volcano Plot - lncRNAs")

# protein_coding_result (ext_gene as label)
create_volcano_plot(protein_coding_result,
                    title = "Volcano Plot - Protein coding genes")

# NA_result (target_id as label)
create_volcano_plot(NA_result,
                    title = "Volcano Plot - Novel genes")

# exploratory analysis
# interactive visualization
sleuth_live(oe)

# PCA to visualize sample relationships based on gene expression
plot_pca(sleuth_object, color_by = 'condition')
# Group density plot to visualize the distribution of gene expression levels across conditions
plot_group_density(sleuth_object,
                   use_filtered = TRUE,
                   units = "est_counts",
                   trans = "log",
                   grouping = setdiff(colnames(sleuth_object$sample_to_covariates),
                                      "sample"), offset = 1)

# check if found genes are present in the paper
df_one <- read.csv("sleuth_results_significant.csv")
df_two <- read.csv("parental_vs_para_from_paper.csv")

unique_ens_gene <- unique(df_one$ens_gene)
unique_ensembl_gene_id <- unique(df_two$ensembl_gene_id)

cat("Unique values in df_one$ens_gene:", length(unique_ens_gene), "\n")
cat("Unique values in df_two$ensembl_gene_id:", length(unique_ensembl_gene_id), "\n")

# Display values that are in df_one but not in df_two
ens_gene_not_in_df_two <- setdiff(unique_ens_gene, unique_ensembl_gene_id)
cat("Values in df_one$ens_gene not found in df_two$ensembl_gene_id:", length(ens_gene_not_in_df_two), "\n")
ens_gene_not_in_df_two

#normalized expression values and estimated counts as the expression units
#counts_per_replicate <- sleuth_to_matrix(oe, "obs_norm", "est_counts")
counts_per_replicate <- sleuth_to_matrix(oe, "obs_norm", "scaled_reads_per_base")
write.csv(counts_per_replicate, file = "counts_per_replicate.csv", row.names = TRUE)
# Read the data from the CSV file (adjust the file path accordingly)
## for all significant top 20
data <- read.csv("counts_per_replicate_gene_level.csv", row.names = 1)
data_filtered <- data[rownames(data) %in% sleuth_sihnificant_top_20$target_id, ]

data_matrix <- as.matrix(data_filtered)

# Create a heatmap using pheatmap
pheatmap(data_matrix, 
         cluster_cols = FALSE, # Turn off column clustering
         scale = "row",        # Scale rows (genes) by default
         main = "Significant genes top 20",
         annotation_col = NULL)  # Turn off column annotations

## for lncRNAs top 20
data <- read.csv("counts_per_replicate.csv", row.names = 1)
data_filtered <- data[rownames(data) %in% lncRNA_result_top_20$target_id, ]

data_matrix <- as.matrix(data_filtered)

# Create a heatmap using pheatmap
pheatmap(data_matrix, 
         cluster_cols = FALSE, # Turn off column clustering
         scale = "row",        # Scale rows (genes) by default
         main = "lncRNAs top 20",
         annotation_col = NULL)  # Turn off column annotations

## for protein coding top 20
data <- read.csv("counts_per_replicate.csv", row.names = 1)
data_filtered <- data[rownames(data) %in% protein_coding_result_top_20$target_id, ]

data_matrix <- as.matrix(data_filtered)

# Create a heatmap using pheatmap
pheatmap(data_matrix, 
         cluster_cols = FALSE, # Turn off column clustering
         scale = "row",        # Scale rows (genes) by default
         main = "Protein coding genes top 20",
         annotation_col = NULL)  # Turn off column annotations

## for NA_result_top_20
data <- read.csv("counts_per_replicate.csv", row.names = 1)
data_filtered <- data[rownames(data) %in% NA_result_top_20$target_id, ]

data_matrix <- as.matrix(data_filtered)

# Create a heatmap using pheatmap
pheatmap(data_matrix, 
         cluster_cols = FALSE, # Turn off column clustering
         scale = "row",        # Scale rows (genes) by default
         main = "Novel genes top 20",
         annotation_col = NULL)  # Turn off column annotations

# plot bootstrap , to account for kallisto estimation ENST00000257555.11
target_id <- "ENST00000257555.11"  # Replace with the actual transcript or gene identifier
units <- "est_counts"
color_by <- "condition"  # Replace with the actual grouping variable in your sample metadata
x_axis_angle <- 50
divide_groups <- TRUE

# Call the plot_bootstrap function
plot_bootstrap(oe, target_id, units = units, color_by = color_by,
               x_axis_angle = x_axis_angle, divide_groups = divide_groups)
