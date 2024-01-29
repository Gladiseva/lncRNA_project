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

df_one <- read.csv("sleuth_results_significant_gl.csv")
df_two <- read.csv("parental_vs_para_from_paper.csv")
# get gene ids from my data
gene_id_and_log2fch_my <- data.frame(
  ensembl_gene_id = df_one$target_id,
  log2FoldChange = ifelse(df_one$b > 0, "increase", "decrease")
)

# get gene ids from paper
gene_id_and_log2fch_paper <- data.frame(
  ensembl_gene_id = df_two$ensembl_gene_id,
  log2FoldChange = ifelse(df_two$log2FoldChange > 0, "decrease", "increase")
)

# merge left join my data with paper data, on ensembl_gene_id
merged_data <- merge(unique(gene_id_and_log2fch_my), gene_id_and_log2fch_paper, by = "ensembl_gene_id", all.x = TRUE)

# count differences in expression
differences <- merged_data[merged_data$log2FoldChange.x != merged_data$log2FoldChange.y, ]
differences <- differences[complete.cases(differences), ]
dim(differences)[1]

