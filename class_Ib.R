# ================================================================
# AI and Omics Research Internship 2025 - Class Ib
# Author: Youssra Boumait
# Task: Load and clean patient_info.csv dataset
# ================================================================

# ✅ 1. S

# ✅ 2. Create subfolders (run only once)
dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")

# ✅ 3. Load the dataset
patient_data <- read.csv("raw_data/patient_info.csv", stringsAsFactors = FALSE)

# ✅ 4. Inspect the structure
head(patient_data)
str(patient_data)
summary(patient_data)

# ✅ 5. Convert variable types
patient_data$gender <- as.factor(patient_data$gender)
patient_data$diagnosis <- as.factor(patient_data$diagnosis)

# ✅ 6. Create binary factor for smoker (Yes = 1, No = 0)
patient_data$smoker_bin <- ifelse(patient_data$smoker == "Yes", 1, 0)
patient_data$smoker_bin <- as.factor(patient_data$smoker_bin)

# ✅ 7. Save cleaned dataset
write.csv(patient_data, "clean_data/patient_info_clean.csv", row.names = FALSE)

# ✅ 8. Save workspace (optional)
save.image(file = "results/full_workspace.RData")
