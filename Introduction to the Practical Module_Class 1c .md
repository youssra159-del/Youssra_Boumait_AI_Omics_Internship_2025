# ====================== 1.1. Check Cholesterol level (using if) ================
colesterol_level <- 240

if (colesterol_level > 230) {
  print("High Cholesterol")
}


# ====================== 2. Blood Pressure Status (using if...else)==========
# Write an if…else statement to check if blood pressure is normal.
# If it’s less than 120, print: “Blood Pressure is normal”
# If false then print: “Blood Pressure is high”

Systolic_bp <- 130


if (Systolic_bp < 120) {
  print("Blood Pressure is normal")
} else {
  print("Blood Pressure is high")
}

  
# ============ 3. Automating Data Type Conversion with for loop ===============
  
  # Use patient_info.csv data and metadata.csv
# Import the datasets
patient_info <- read.csv("patient_info.csv")
metadata <- read.csv("metadata.csv")
# Import the datasets
patient_info <- read.csv("patient_info.csv")
metadata <- read.csv("metadata.csv")

# Create copies
data1 <- patient_info
data2 <- Metadata

# Identify columns to convert (example: adjust according to dataset structure)
factor_cols1 <- c("gender", "smoker", "diagnosis")
factor_cols2 <- c("name", "height", "gender")

# Convert to factors using for loop
for (col in factor_cols1) {
  data1[[col]] <- as.factor(data1[[col]])
}

# Convert them to factors
for (col in factor_cols2) {
  data2[[col]] <- as.factor(data2[[col]])
}

# Verify changes
str(data1)
str(data2)


# ================4. Converting Factors to Numeric Codes===============
# Application on patient_info data
binary_cols1 <- c("smoker", "gender", "diagnosis")


for (col in binary_cols1) {
  if (col == "smoker") {
    data1[[col]] <- ifelse(data1[[col]] == "Yes", 1, 0)
  } else if (col == "gender") {
    data1[[col]] <- ifelse(data1[[col]] == "Male", 1, 0)
  } else if (col == "diagnosis") {
    data1[[col]] <- ifelse(data1[[col]] == "Cancer", 1, 0)
  }
}

# Verify
str(data1)
head(data1)


# Verify
str(patient_info)  # original
str(data1)         # modified

# Application on metadata

# Identify character columns in Metadata
factor_cols2 <- c("name", "height", "gender")

# Convert them to factors
for (col in factor_cols2) {
  data2[[col]] <- as.factor(data2[[col]])
}

# Verify changes
str(data2)


# Verify changes
str(data2)
str(Metadata)
