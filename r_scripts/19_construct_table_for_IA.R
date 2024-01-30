# Read sleuth_results_significant_tl.csv
all_tl <- read.csv("sleuth_results_full_tl.csv")

# Keep only the specified columns
columns_selected_stl <- all_tl[, c("target_id", "ens_gene", "ext_gene", "gene_biotype", "transcript_biotype", "qval", "b")]

# Read counts_per_replicate_tl.csv
counts_per_replicate <- read.csv("counts_per_replicate_tl.csv")

# Merge the data frames based on "target_id" (assuming it's the common column)
merged_data <- merge(columns_selected_stl, counts_per_replicate, by.x = "target_id", by.y = "transcript_id", all.x = TRUE)
# lncRNAs
lncRNA_rows <- merged_data[merged_data$gene_biotype == "lncRNA" | merged_data$transcript_biotype == "lncRNA", ]
overlap3prime <- read.table("annotated_overlap3prime.bed", header = FALSE, col.names = c("Column1", "Column2", "Column3", "Column4", "Column5", "Column6"))
overlap5prime <- read.table("annotated_overlap3prime.bed", header = FALSE, col.names = c("Column1", "Column2", "Column3", "Column4", "Column5", "Column6"))

# Extract values from the 4th column (overlap3prime.bed)
values_from_column4 <- overlap3prime$Column4
lncRNA_rows$overlap3prime <- lncRNA_rows$target_id %in% values_from_column4

# Extract values from the 4th column (overlap5prime.bed)
values_from_column4 <- overlap5prime$Column4
lncRNA_rows$overlap5prime <- lncRNA_rows$target_id %in% values_from_column4
lncRNA_rows <- subset(lncRNA_rows, !(is.na(qval) & is.na(b)))
write.csv(lncRNA_rows, file = "lncRNA_integrative_analysis.csv", row.names = FALSE)
# Novel
MSTRG_rows <- merged_data[grep("MSTRG", merged_data$target_id), ]
MSTRG_rows <- MSTRG_rows[, c("target_id", "ens_gene", "qval", "b", "X3_2","X3_4","X3_7","P1","P2","P3")]
# Read the data from the file
novel_intergenic <- read.table("novel_intergenic.bed", header = FALSE, col.names = c("Column1", "Column2", "Column3", "Column4", "Column5", "Column6"))
overlap3prime <- read.table("overlap3prime.bed", header = FALSE, col.names = c("Column1", "Column2", "Column3", "Column4", "Column5", "Column6"))
overlap5prime <- read.table("overlap3prime.bed", header = FALSE, col.names = c("Column1", "Column2", "Column3", "Column4", "Column5", "Column6"))

# Extract values from the 4th column (novel_intergenic.bed)
values_from_column4 <- novel_intergenic$Column4
MSTRG_rows$intragenic <- MSTRG_rows$target_id %in% values_from_column4

# Extract values from the 4th column (overlap3prime.bed)
values_from_column4 <- overlap3prime$Column4
MSTRG_rows$overlap3prime <- MSTRG_rows$target_id %in% values_from_column4

# Extract values from the 4th column (overlap5prime.bed)
values_from_column4 <- overlap5prime$Column4
MSTRG_rows$overlap5prime <- MSTRG_rows$target_id %in% values_from_column4

#Coding potential
novel_data <- read.table("novel_coding_potential.dat", header = FALSE, stringsAsFactors = FALSE)

novel_data$short_target_id <- sub("::.*", "", novel_data$V1)
# Assuming MSTRG_rows is your existing data frame
for (i in 1:nrow(MSTRG_rows)) {
  target_id <- MSTRG_rows$target_id[i]
  novel_row <- novel_data[novel_data$short_target_id == target_id, ]
  
  if (nrow(novel_row) > 0) {
    MSTRG_rows$coding_potential[i] <- novel_row$V5
  }
  MSTRG_rows$cutoff_CPAT[i] <- MSTRG_rows$coding_potential[i] > 0.364
}
MSTRG_rows <- subset(MSTRG_rows, !(is.na(qval) & is.na(b)))
write.csv(MSTRG_rows, file = "novel_integrative_analysis.csv", row.names = FALSE) 

library(ggplot2)
par(mfrow = c(2, 2))

pie_data <- data.frame(
  Condition = c("Intragenic with both overlaps", "Intragenic with no overlap"),
  Count = c(
    sum(MSTRG_rows$overlap5prime == TRUE & MSTRG_rows$overlap3prime == TRUE & MSTRG_rows$intragenic == TRUE),
    sum(MSTRG_rows$overlap5prime == FALSE & MSTRG_rows$overlap3prime == FALSE & MSTRG_rows$intragenic == TRUE)
  )
)

# Create the pie chart
ggplot(pie_data, aes(x = "", y = Count, fill = Condition)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(title = "Pie Chart for Overlap Conditions",
       fill = "Condition") +
  theme_void()

pie_data <- data.frame(
  Condition = c("Intragenic", "Not intragenic"),
  Count = c(
    sum(MSTRG_rows$intragenic == TRUE),
    sum(MSTRG_rows$intragenic == FALSE)
  )
)

# Create the pie chart
ggplot(pie_data, aes(x = "", y = Count, fill = Condition)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(title = "Pie Chart for Overlap Conditions",
       fill = "Condition") +
  theme_void()

pie_data <- data.frame(
  Condition = c("Both overlaps", "No overlap"),
  Count = c(
    sum(MSTRG_rows$overlap5prime == TRUE & MSTRG_rows$overlap3prime == TRUE),
    sum(MSTRG_rows$overlap5prime == FALSE & MSTRG_rows$overlap3prime == FALSE)
  )
)

# Create the pie chart
ggplot(pie_data, aes(x = "", y = Count, fill = Condition)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(title = "Pie Chart for Overlap Conditions",
       fill = "Condition") +
  theme_void()

pie_data <- data.frame(
  Condition = c("Both overlaps", "No overlaps"),
  Count = c(
    sum(lncRNA_rows$overlap5prime == TRUE & lncRNA_rows$overlap3prime == TRUE),
    sum(lncRNA_rows$overlap5prime == FALSE & lncRNA_rows$overlap3prime == FALSE)
  )
)

# Create the pie chart
ggplot(pie_data, aes(x = "", y = Count, fill = Condition)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(title = "Pie Chart for Overlap Conditions",
       fill = "Condition") +
  theme_void()
