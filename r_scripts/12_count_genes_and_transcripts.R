install.packages("jsonlite")
library(jsonlite)

# Function to process each folder
process_folder <- function(folder_path) {
  # Load the abundance.tsv file
  abundance <- read.table(file.path(folder_path, "abundance.tsv"), header = TRUE, sep = "\t")
  
  filtered_abundance_ENS <- subset(abundance, grepl("ENS", target_id) & tpm > 0.0)
  
  # Calculate the length of unique target_id
  total_transcripts_genes <- dim(unique(filtered_abundance_ENS))[1]
  
  filtered_abundance_MST <- subset(abundance, grepl("MST", target_id) & tpm > 0.0)
  
  # Calculate the length of unique target_id
  total_novel <- dim(unique(filtered_abundance_MST))[1]
  
  # Create a data frame with the results
  result_df <- data.frame(
    Folder = basename(folder_path),
    Total_Transcripts_Genes = total_transcripts_genes,
    Novel_Transcripts = total_novel
  )
  
  return(result_df)
}

# Get a list of folders in the "results" directory
results_folders <- list.dirs("results", full.names = TRUE, recursive = FALSE)

# Process each folder and bind the results
all_results <- do.call(rbind, lapply(results_folders, process_folder))

# Save the results to a CSV file
write.csv(all_results, "summary_kallisto_results.csv", row.names = FALSE)
