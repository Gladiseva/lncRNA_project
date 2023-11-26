## Instalation
#source("http://bioconductor.org/biocLite.R")
#BiocManager::install("rhdf5")
#install.packages("devtools")
#devtools::install_github("pachterlab/sleuth")
library("sleuth")

#specify where the kallisto results are stored.
sample_id <- dir(file.path(".", "results"))

#paths to the kallisto results indexed by the sample IDs is collated with
kallisto_dirs <- file.path(".", "results", sample_id)

#load file that describes the experimental design
sample2condition <- read.table(file.path(".", "NSCLC_exp_design.txt"),
                               header = TRUE,
                               stringsAsFactors=FALSE)

#directories must be appended in a new column
sample2condition <- dplyr::mutate(sample2condition, path = kallisto_dirs)

#sleuth object will store not only the information about the experiment,
#but also details of the model to be used for differential testing, and the results.
#It is prepared and used with four commands that 
#(1) load the kallisto processed data into the object 
#(2) estimate parameters for the sleuth response error measurement (full) model 
#(3) estimate parameters for the sleuth reduced model
#(4) perform differential analysis (testing) using the likelihood ratio test.
sleuth_object <- sleuth_prep(sample2condition, extra_bootstrap_summary = TRUE)
sleuth_object <- sleuth_fit(sleuth_object, ~condition, 'full')
sleuth_object <- sleuth_fit(sleuth_object, ~1, 'reduced')
#2 NA values were found during variance shrinkage estimation due to mean observation values outside of the range used for the LOESS fit.
#The LOESS fit will be repeated using exact computation of the fitted surface to extrapolate the missing values.
#These are the target ids with NA values: ENST00000361624.2, ENST00000387347.2
sleuth_object <- sleuth_lrt(sleuth_object, 'reduced', 'full')
models(sleuth_object)

sleuth_table <- sleuth_results(sleuth_object, 'reduced:full', 'lrt', show_all = FALSE)
#top 20 significant genes with a (Benjamini-Hochberg multiple testing corrected) q-value <= 0.05.
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 20)

plot_bootstrap(sleuth_object, "ENST00000223642.3", units = "est_counts", color_by = "condition")

# add gene names from ENSEMBL using biomaRt
install.packages("BiocManager")
BiocManager::install("biomaRt")

#collect gene names with
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
sleuth_object <- sleuth_prep(sample2condition,
                             target_mapping = t2g,
                             read_bootstrap_tpm = TRUE,
                             extra_bootstrap_summary = TRUE)
sleuth_object <- sleuth_fit(sleuth_object, ~condition, 'full')
sleuth_object <- sleuth_fit(sleuth_object, ~1, 'reduced')
sleuth_object <- sleuth_lrt(sleuth_object, 'reduced', 'full')

#exploratory analysis
sleuth_live(sleuth_object)
plot_pca(sleuth_object, color_by = 'condition')
plot_group_density(sleuth_object,
                   use_filtered = TRUE,
                   units = "est_counts",
                   trans = "log",
                   grouping = setdiff(colnames(sleuth_object$sample_to_covariates),
                                                     "sample"), offset = 1)
