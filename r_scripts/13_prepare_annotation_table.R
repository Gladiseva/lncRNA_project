#Installation
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("biomaRt")

if (!require("dplyr")) install.packages("dplyr")
library(dplyr)

# get gene data from ENSEMBL using biomaRt (GRCh38)
# collect gene names with
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = "ensembl.org")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name", "gene_biotype",
                                     "transcript_biotype"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# Assuming you have your custom code to read GTF data
gtf_path <- "ALL_merged.gtf"

# Read the GTF file, skipping metadata lines
gtf_data <- readLines(gtf_path)
gtf_data <- gtf_data[!grepl("^#", gtf_data)]

# Extract necessary information from GTF data
gene_info <- data.frame(
  target_id = sub('.*transcript_id "(.*?)";.*', '\\1', gtf_data),
  ens_gene = sub('.*gene_id "(.*?)";.*', '\\1', gtf_data)
)

gene_info <- unique(gene_info)

gene_info <- gene_info %>%
  mutate(modified_target_id = sub("\\.\\d+$", "", target_id))

annotation_data <- gene_info %>%
  left_join(t2g, by = c("modified_target_id" = "target_id"))

annotation_data <- annotation_data %>%
  mutate(ens_gene = coalesce(ens_gene.y, ens_gene.x)) %>%
  select(target_id, ens_gene, modified_target_id, ext_gene, gene_biotype, transcript_biotype)

# Save to a CSV file or any desired format
write.csv(annotation_data, "annotation_table_from_ALLandENS.csv", row.names = FALSE)
