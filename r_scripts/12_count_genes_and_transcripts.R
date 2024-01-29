install.packages("jsonlite")
library(jsonlite)

# Function to process each folder
process_folder <- function(folder_path) {
  # Load the abundance.tsv file
  abundance <- read.table(file.path(folder_path, "abundance.tsv"), header = TRUE, sep = "\t")
  
  filtered_abundance_ENS <- subset(abundance, grepl("ENS", target_id) & tpm > 0.0)
  
  # Calculate the length of unique target_id
  total_unique_transcripts <- length(unique(filtered_abundance_ENS$target_id))
  total_unique_genes <- length(unique(sub("\\.\\d+$", "", filtered_abundance_ENS$target_id)))
  
  
  filtered_abundance_MST <- subset(abundance, grepl("MST", target_id) & tpm > 0.0)
  
  # Calculate the length of unique target_id
  total_novel_unique_transcripts <- length(unique(filtered_abundance_MST$target_id))
  total_novel_unique_genes <- length(unique(sub("\\.\\d+$", "", filtered_abundance_MST$target_id)))
  
  # Create a data frame with the results
  result_df <- data.frame(
    Folder = basename(folder_path),
    Total_Genes = total_unique_genes,
    Total_Transcripts = total_unique_transcripts,
    Novel_Transcripts = total_novel_unique_transcripts,
    Novel_genes = total_novel_unique_genes
  )
  
  return(result_df)
}

# Get a list of folders in the "results" directory
results_folders <- list.dirs("results", full.names = TRUE, recursive = FALSE)

# Process each folder and bind the results
all_results <- do.call(rbind, lapply(results_folders, process_folder))

# Save the results to a CSV file
write.csv(all_results, "summary_kallisto_results_2.csv", row.names = FALSE)
