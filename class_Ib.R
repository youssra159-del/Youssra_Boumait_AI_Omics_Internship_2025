# ================================================================
# AI and Omics Research Internship 2025 - Class Ib
# Author: Youssra Boumait
# ================================================================

# 1. Set Working Directory
getwd()


#2. Create subfolders 
dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")
dir.create("Histogram")

# 3. Download "patient_info.csv" dataset from GitHub repository
patient_data <- read.csv("raw_data/patient_info.csv")

# 4. Inspect the structure
head(patient_data)
str(patient_data)
summary(patient_data)

# 5. Identify variables with incorrect or inconsistent data types
patient_data$gender <- as.factor(patient_data$gender)
patient_data$diagnosis <- as.factor(patient_data$diagnosis)

# 6. Créer une variable binaire pour "smoker"
patient_data$smoker_bin <- ifelse(patient_data$smoker == "Yes", 1, 0)
patient_data$smoker_bin <- as.factor(patient_data$smoker_bin)

# 7. Save cleaned dataset
write.csv(patient_data, "clean_data/patient_info_clean.csv", row.names = FALSE)

# 8. Saving entire Workspace
save.image(file = "Youssra_Boumait_Class_Ib_Assignment.RData")
