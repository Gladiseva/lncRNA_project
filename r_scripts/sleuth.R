## Instalation
BiocManager::install("rhdf5")
install.packages("devtools")
devtools::install_github("pachterlab/sleuth")
install.packages("BiocManager")
BiocManager::install("biomaRt")
library("sleuth")
library(ggplot2)
install.packages("ggrepel")
library(ggrepel)
install.packages("readxl")
library(readxl)
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
                         host = "ensembl.org")
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
# NA values were found during variance shrinkage estimation LOESS
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
head(sleuth_significant, 20)

sleuth_results_oe
write.csv(sleuth_results_oe, file = "sleuth_results_full.csv", row.names = FALSE)
write.csv(sleuth_significant, file = "sleuth_results_significant.csv", row.names = FALSE)

# Print or further process the unique gene biotypes
unique_gene_biotypes <- unique(sleuth_results_oe$gene_biotype)
print(unique_gene_biotypes)

# Print or further process the unique transcript biotypes
unique_transcript_biotypes <- unique(sleuth_results_oe$transcript_biotype)
print(unique_transcript_biotypes)

# lncRNA
lncRNA_result <- subset(sleuth_significant, gene_biotype == "lncRNA" | transcript_biotype == "lncRNA")
head(lncRNA_result)
write.csv(lncRNA_result, file = "sleuth_results_lncRNA_significant.csv", row.names = FALSE)

# protein_coding
protein_coding_result <- subset(sleuth_significant, gene_biotype == "protein_coding" | transcript_biotype == "protein_coding")
head(protein_coding_result)
write.csv(protein_coding_result, file = "sleuth_results_protein_coding_significant.csv", row.names = FALSE)

# NAs
NA_coding_result <- subset(sleuth_significant, is.na(gene_biotype) | is.na(transcript_biotype))
head(NA_coding_result)
dim(NA_coding_result)
write.csv(NA_coding_result, file = "sleuth_results_NAs_significant.csv", row.names = FALSE)

sleuth_live(oe)

# Create a volcano plot
significance_threshold <- 0.05
ggplot(sleuth_significant, aes(x = b, y = -log10(pval))) +
  geom_point(aes(color = qval < significance_threshold), alpha = 0.5, size = 2) +  # Highlight significant points
  geom_text_repel(aes(label = ext_gene), 
                  data = subset(sleuth_significant, qval < significance_threshold), 
                  box.padding = 0.5, point.padding = 0.2, segment.color = "grey50") +  # Add gene names for significant points
  scale_color_manual(values = c("blue", "red")) +       # Customize colors
  labs(title = "Volcano Plot", x = "Log2-Fold Change (b)", y = "-log10(p-value)") +
  theme_minimal()
ggplot(lncRNA_result, aes(x = b, y = -log10(pval))) +
  geom_point(aes(color = qval < significance_threshold), alpha = 0.5, size = 2) +  # Highlight significant points
  geom_text_repel(aes(label = ext_gene), 
                  data = subset(lncRNA_result, qval < significance_threshold), 
                  box.padding = 0.5, point.padding = 0.2, segment.color = "grey50") +  # Add gene names for significant points
  scale_color_manual(values = c("blue", "red")) +       # Customize colors
  labs(title = "Volcano Plot", x = "Log2-Fold Change (b)", y = "-log10(p-value)") +
  theme_minimal()
ggplot(protein_coding_result, aes(x = b, y = -log10(pval))) +
  geom_point(aes(color = qval < significance_threshold), alpha = 0.5, size = 2) +  # Highlight significant points
  geom_text_repel(aes(label = ext_gene), 
                  data = subset(protein_coding_result, qval < significance_threshold), 
                  box.padding = 0.5, point.padding = 0.2, segment.color = "grey50") +  # Add gene names for significant points
  scale_color_manual(values = c("blue", "red")) +       # Customize colors
  labs(title = "Volcano Plot", x = "Log2-Fold Change (b)", y = "-log10(p-value)") +
  theme_minimal()
ggplot(NA_coding_result, aes(x = b, y = -log10(pval))) +
  geom_point(aes(color = qval < significance_threshold), alpha = 0.5, size = 2) +  # Highlight significant points
  geom_text_repel(aes(label = target_id), 
                  data = subset(NA_coding_result, qval < significance_threshold), 
                  box.padding = 0.5, point.padding = 0.2, segment.color = "grey50") +  # Add gene names for significant points
  scale_color_manual(values = c("blue", "red")) +       # Customize colors
  labs(title = "Volcano Plot", x = "Log2-Fold Change (b)", y = "-log10(p-value)") +
  theme_minimal()


# exploratory analysis
plot_bootstrap(sleuth_object, "ENST00000223642.3", units = "est_counts", color_by = "condition")
sleuth_live(sleuth_object)
plot_pca(sleuth_object, color_by = 'condition')
plot_group_density(sleuth_object,
                   use_filtered = TRUE,
                   units = "est_counts",
                   trans = "log",
                   grouping = setdiff(colnames(sleuth_object$sample_to_covariates),
                                      "sample"), offset = 1)

#comparison with data from paper
df = read_excel("1-s2.0-S147655861830232X-mmc3.xlsx", sheet=4)
write.csv(df, gsub("xlsx", "4.csv", "1-s2.0-S147655861830232X-mmc3.xlsx"), row.names=FALSE)

#df_one <- read.csv("sleuth_results_full.csv")
df_one <- read.csv("sleuth_results_significant.csv")
df_two <- read.csv("parental_vs_para_from_paper.csv")

# Merge the two data frames based on the specified columns
merged_df <- merge(df_one, df_two, by.x = "ens_gene", by.y = "ensembl_gene_id", all.x = TRUE)

# Count the number of entries from the first CSV that have a match in the second CSV
matching_entries <- sum(!is.na(merged_df$hgnc_symbol))

# Print the result
cat("Number of entries from the first CSV found in the second CSV:", matching_entries, "\n")


unique_ens_gene <- unique(df_one$ens_gene)
unique_ensembl_gene_id <- unique(df_two$ensembl_gene_id)

cat("Unique values in df_one$ens_gene:", length(unique_ens_gene), "\n")
cat("Unique values in df_two$ensembl_gene_id:", length(unique_ensembl_gene_id), "\n")

# Display values that are in df_one but not in df_two
ens_gene_not_in_df_two <- setdiff(unique_ens_gene, unique_ensembl_gene_id)
cat("Values in df_one$ens_gene not found in df_two$ensembl_gene_id:", length(ens_gene_not_in_df_two), "\n")
ens_gene_not_in_df_two

# Display values that are in df_two but not in df_one
ens_gene_not_in_df_one <- setdiff(unique_ensembl_gene_id, unique_ens_gene)
cat("Values in df_two$ens_gene not found in df_one$ens_gene:", length(ens_gene_not_in_df_one), "\n")
ens_gene_not_in_df_one

