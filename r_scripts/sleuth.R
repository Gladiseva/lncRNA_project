## Instalation
# Install and load packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("rhdf5", "sleuth", "biomaRt"))
if (!require("devtools")) install.packages("devtools")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggrepel")) install.packages("ggrepel")
if (!require("readxl")) install.packages("readxl")
install.packages("svglite")

library("rhdf5")
library("devtools")
library("sleuth")
library("BiocManager")
library("biomaRt")
library("ggplot2")
library("ggrepel")
library("readxl")
library(svglite)
install.packages("pheatmap")
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

# sleuth object will store the information about the experiment,
# and details of the model to be used for differential testing and the results.
# It is prepared and used with four commands that
#(1) load the kallisto processed data into the object
#(2) estimate parameters for the sleuth response error measurement (full) model
#(3) estimate parameters for the sleuth reduced model
#(4) perform differential analysis (testing) using the likelihood ratio test.

# add gene names from ENSEMBL using biomaRt (GRCh38)
# collect gene names with
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = "https://ensembl.org")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name", "gene_biotype",
                                     "transcript_biotype"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
# By default the transformation of counts is natural log
sleuth_object <- sleuth_prep(sample2condition,
                             target_mapping = t2g,
                             read_bootstrap_tpm = TRUE,
                             extra_bootstrap_summary = TRUE,
                             transformation_function = function(x) log2(x + 0.5))
sleuth_object <- sleuth_fit(sleuth_object, ~condition, "full")
sleuth_object <- sleuth_fit(sleuth_object, ~1, "reduced")
# OUTPUT: NA values were found during variance shrinkage estimation LOESS
# These are the target ids with NA values: ENST00000361624.2, ENST00000387347.2
sleuth_object <- sleuth_lrt(sleuth_object, "reduced", "full")

models(sleuth_object)

# Test significant differences between conditions using the Wald test
# Wald test for differential expression of isoforms. var oe -> observed2expected
oe <- sleuth_wt(sleuth_object, which_beta = "conditionparental")
sleuth_results_oe <- sleuth_results(oe,
                                    test = "conditionparental",
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
print(unique_gene_biotypes)
print(unique_transcript_biotypes)

# lncRNA
lncRNA_result <- subset(sleuth_results_oe, gene_biotype == "lncRNA" | transcript_biotype == "lncRNA")
lncRNA_result_top_20 <- head(lncRNA_result, 20)
write.csv(lncRNA_result, file = "sleuth_results_lncRNA.csv", row.names = FALSE)

lncRNA_result_significant <- subset(sleuth_significant, gene_biotype == "lncRNA" | transcript_biotype == "lncRNA")
lncRNA_significant_result_top_20 <- head(lncRNA_result_significant, 20)
# protein_coding
protein_coding_result <- subset(sleuth_results_oe, gene_biotype == "protein_coding" | transcript_biotype == "protein_coding")
protein_coding_result_top_20 <- head(protein_coding_result, 20)
write.csv(protein_coding_result, file = "sleuth_results_protein_coding.csv", row.names = FALSE)

# NAs
NA_result <- subset(sleuth_results_oe, is.na(gene_biotype) | is.na(transcript_biotype))
NA_result_top_20 <- head(NA_result, 20)
write.csv(NA_result, file = "sleuth_results_NAs.csv", row.names = FALSE)

# Create a volcano plot
create_volcano_plot <- function(data, label_column, significance_threshold = 0.05, title = "Volcano Plot") {
  ggplot(data, aes(x = b, y = -log10(pval))) +
    geom_point(aes(color = qval < significance_threshold), alpha = 0.5, size = 2) +
    geom_text_repel(aes(label = get(label_column)), 
                    data = subset(data, qval < significance_threshold), 
                    box.padding = 0.5, point.padding = 0.2, segment.color = "grey50") +
    scale_color_manual(values = c("blue", "red")) +
    labs(title = title, x = "Log2-Fold Change (b)", y = "-log10(p-value)") +
    theme_minimal()
}

# all (ext_gene as label)
create_volcano_plot(sleuth_results_oe, label_column = "ext_gene")

# lncRNA_result (ext_gene as label)
create_volcano_plot(lncRNA_result, label_column = "ext_gene")

# protein_coding_result (ext_gene as label)
create_volcano_plot(protein_coding_result, label_column = "ext_gene")

# NA_result (target_id as label)
create_volcano_plot(NA_result, label_column = "target_id")

# exploratory analysis
# interactive visualization
sleuth_live(oe)

# Plot bootstrap estimates for transcript
plot_bootstrap(sleuth_object, "ENST00000223642.3", units = "est_counts", color_by = "condition")
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


# occurencies per gene
sleuth_genes <- sleuth_gene_table(oe, 'conditionparental', test_type ='wt',
                                  which_group = 'ext_gene')
kallisto_output <- read_kallisto("results/3_2/abundance.tsv", read_bootstrap = FALSE)
#normalized expression values and estimated counts as the expression units
counts_per_replicate <- sleuth_to_matrix(oe, "obs_norm", "est_counts")
write.csv(counts_per_replicate, file = "counts_per_replicate.csv", row.names = TRUE)
# Read the data from the CSV file (adjust the file path accordingly)
## for all significant top 20
data <- read.csv("counts_per_replicate.csv", row.names = 1)
data_filtered <- data[rownames(data) %in% sleuth_sihnificant_top_20$target_id, ]

data_matrix <- as.matrix(data_filtered)

# Create a heatmap using pheatmap
pheatmap(data_matrix, 
         cluster_cols = FALSE, # Turn off column clustering
         scale = "row",        # Scale rows (genes) by default
         main = "Heatmap Example",
         annotation_col = NULL)  # Turn off column annotations

## for lncRNAs top 20
data <- read.csv("counts_per_replicate.csv", row.names = 1)
data_filtered <- data[rownames(data) %in% lncRNA_result_top_20$target_id, ]

data_matrix <- as.matrix(data_filtered)

# Create a heatmap using pheatmap
pheatmap(data_matrix, 
         cluster_cols = FALSE, # Turn off column clustering
         scale = "row",        # Scale rows (genes) by default
         main = "Heatmap Example",
         annotation_col = NULL)  # Turn off column annotations

## for protein coding top 20
data <- read.csv("counts_per_replicate.csv", row.names = 1)
data_filtered <- data[rownames(data) %in% protein_coding_result_top_20$target_id, ]

data_matrix <- as.matrix(data_filtered)

# Create a heatmap using pheatmap
pheatmap(data_matrix, 
         cluster_cols = FALSE, # Turn off column clustering
         scale = "row",        # Scale rows (genes) by default
         main = "Heatmap Example",
         annotation_col = NULL)  # Turn off column annotations

## for NA_result_top_20
data <- read.csv("counts_per_replicate.csv", row.names = 1)
data_filtered <- data[rownames(data) %in% NA_result_top_20$target_id, ]

data_matrix <- as.matrix(data_filtered)

# Create a heatmap using pheatmap
pheatmap(data_matrix, 
         cluster_cols = FALSE, # Turn off column clustering
         scale = "row",        # Scale rows (genes) by default
         main = "Heatmap Example",
         annotation_col = NULL)  # Turn off column annotations
