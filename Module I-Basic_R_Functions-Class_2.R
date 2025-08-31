# ===================================================================
#               AI and Biotechnology / Bioinformatics
# ===================================================================
# -------------------------------------------------------------------
#             AI and Omics Research Internship (2025)
# -------------------------------------------------------------------
#                Module I: Getting Started with R
# -------------------------------------------------------------------
#                           Assignment 2
# ===================================================================
# Define input folder
input_dir <- "C:/Users/Lenovo/OneDrive/Documents/AI_Omics_Internship_2025/Module_I/raw_data"
output_dir <- "C:/Users/Lenovo/OneDrive/Documents/AI_Omics_Internship_2025/Module_I/results"

# # Import dataset and Read file
data1 <- read.csv(file.path(input_dir, "DEGs_Data_1.csv"), header = TRUE)
data2 <- read.csv(file.path(input_dir, "DEGs_Data_2.csv"), header = TRUE)


# Files to process (must match your real CSV file names)
files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")

# Container to store results
result_list <- list()

# Loop over datasets
for (file_name in files_to_process) {
  cat("\nProcessing:", file_name, "\n")
  
  # 1. Read CSV
  input_file_path <- file.path(input_dir, file_name)
  data <- read.csv(input_file_path, header = TRUE)
  
  # 2. Replace missing padj with 1
  na_padj_before <- sum(is.na(data$padj))
  if (na_padj_before > 0) {
    data$padj[is.na(data$padj)] <- 1
    cat("Replaced", na_padj_before, "missing padj values with 1.\n")
  }
  
  # 3. Classify genes
  data$status <- mapply(classify_gene, data$logFC, data$padj)
  
  # 4. Add a column with the source filename
  data$source_file <- file_name
  
  # 5. Save processed file
  output_file_path <- file.path(output_dir, paste0("Processed_", file_name))
  write.csv(data, output_file_path, row.names = FALSE)
  cat("Saved:", output_file_path, "\n")
  
  # 6. Summaries
  cat("Status counts for", file_name, ":\n")
  print(table(data$status))
  
  # 7. Store in list
  result_list[[file_name]] <- data
}

results_1 <- result_list[[1]] 
results_2 <- result_list[[2]]
