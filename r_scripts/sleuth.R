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
# Assuming sleuth_results_oe is your data frame
write.csv(sleuth_results_oe, file = "sleuth_results_full.csv", row.names = FALSE)

# Print or further process the unique gene biotypes
unique_gene_biotypes <- unique(sleuth_results_oe$gene_biotype)
print(unique_gene_biotypes)

# Print or further process the unique gene biotypes
unique_transcript_biotypes <- unique(sleuth_results_oe$transcript_biotype)
print(unique_transcript_biotypes)

# lncRNA
lncRNA_result <- subset(sleuth_results_oe, gene_biotype == "lncRNA" | transcript_biotype == "lncRNA")
head(lncRNA_result)
write.csv(lncRNA_result, file = "sleuth_results_lncRNA.csv", row.names = FALSE)

# protein_coding
protein_coding_result <- subset(sleuth_results_oe, gene_biotype == "protein_coding" | transcript_biotype == "protein_coding")
head(protein_coding_result)
write.csv(protein_coding_result, file = "sleuth_results_protein_coding.csv", row.names = FALSE)

sleuth_live(oe)

# Create a volcano plot
significance_threshold <- 0.05
ggplot(sleuth_results_oe, aes(x = b, y = -log10(pval))) +
  geom_point(aes(color = qval < significance_threshold), alpha = 0.5, size = 2) +  # Highlight significant points
  geom_text_repel(aes(label = ext_gene), 
                  data = subset(sleuth_results_oe, qval < significance_threshold), 
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
