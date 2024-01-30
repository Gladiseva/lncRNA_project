novel_data <- read.csv("novel_integrative_analysis.csv")
novel_data$qval <- as.numeric(novel_data$qval)
novel_data$b <- as.numeric(novel_data$b)

# Sort the data
sorted_data <- novel_data[order(-novel_data$intragenic, novel_data$qval, -abs(novel_data$b)), ]

write.csv(sorted_data, file = "SORTED_novel_integrative_analysis.csv", row.names = FALSE) 

lncRNA_data <- read.csv("lncRNA_integrative_analysis.csv")
lncRNA_data$qval <- as.numeric(lncRNA_data$qval)
lncRNA_data$b <- as.numeric(lncRNA_data$b)

# Sort the data
sorted_data <- lncRNA_data[order(lncRNA_data$qval, -abs(lncRNA_data$b)), ]

write.csv(sorted_data, file = "SORTED_lncRNA_integrative_analysis.csv", row.names = FALSE) 