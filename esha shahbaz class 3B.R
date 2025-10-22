getwd()# Change the path below to where your .CEL.gz files are stored


# Check files
list.files(pattern = ".CEL.gz$")




# Install and load R.utils if not installed
if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
library(R.utils)

# Unzip (gunzip) all .CEL.gz files in current working directory
cel_gz_files <- list.files(pattern = ".CEL.gz$", full.names = TRUE)
sapply(cel_gz_files, gunzip, remove = FALSE)

# Check that .CEL files are now visible
list.files(pattern = ".CEL$")





library(affy)

# Read the CEL files into an AffyBatch object
raw_data <- ReadAffy()

# Check summary of the dataset
raw_data



# Delete all compressed CEL files in the current folder
gz_files <- list.files(pattern = "\\.CEL\\.gz$", ignore.case = TRUE, full.names = TRUE)

# Confirm how many .gz files exist
length(gz_files)

# Delete them
file.remove(gz_files)

# Check again to confirm deletion
list.files(pattern = "\\.CEL\\.gz$", ignore.case = TRUE)





# Load the QC package
library(arrayQualityMetrics)

# Create a folder to store the QC report
if (!dir.exists("QC_Raw_Data")) dir.create("QC_Raw_Data")

# Run QC on raw AffyBatch data
arrayQualityMetrics(
  expressionset = raw_data,
  outdir = "QC_Raw_Data",
  force = TRUE,
  do.logtransform = TRUE
)






# Perform RMA normalization
normalized_data <- rma(raw_data)

# Extract normalized expression matrix
exprs_data <- exprs(normalized_data)

# Check dimensions: probes x samples
dim(exprs_data)

# Optional: preview the first few rows
head(exprs_data)

# Save normalized expression matrix to CSV
write.csv(exprs_data, "GSE124226_RMA_Normalized_Expression.csv")






# Load matrixStats for rowMedians function
if (!requireNamespace("matrixStats", quietly = TRUE)) install.packages("matrixStats")
library(matrixStats)

# Compute median intensity for each probe across all samples
row_median <- rowMedians(as.matrix(exprs_data))

# Visualize the distribution
hist(row_median, breaks = 100, main = "Median Intensity Distribution", xlab = "Median log2 Expression", freq = FALSE)
abline(v = 3.5, col = "red", lwd = 2)  # example threshold, adjust if needed

# Filter probes above threshold
threshold <- 3.5
filtered_data <- exprs_data[row_median > threshold, ]

# Check new dimensions
dim(filtered_data)



# Define groups based on sample names
sample_names <- colnames(filtered_data)
groups <- ifelse(grepl("PCOS", sample_names, ignore.case = TRUE), "PCOS", "Control")

# Check group assignment
table(groups)
groups








setwd("C:/Users/khagga/OneDrive/Desktop/AI_Omics_Internship_2025/ModuleIII/RawData")
getwd()
exprs_data <- read.csv("GSE124226_RMA_Normalized_Expression.csv", row.names = 1)
# -------------------------------
# Plots for Normalized Expression
# -------------------------------

# Load required library
library(ggplot2)

# Confirm exprs_data is loaded
if(!exists("exprs_data")) stop("Normalized expression data 'exprs_data' not found!")

# Open a new graphics window (for Windows)
dev.new(width=12, height=8)

# ------------------------------
# 1. Boxplot
# ------------------------------
boxplot(exprs_data,
        main = "Boxplot of RMA Normalized Expression",
        las = 2,        # rotate sample names
        col = "lightblue",
        ylab = "Log2 Expression")

# Save boxplot automatically
png("Boxplot_Normalized.png", width=1200, height=800)
boxplot(exprs_data,
        main = "Boxplot of RMA Normalized Expression",
        las = 2,
        col = "lightblue",
        ylab = "Log2 Expression")
dev.off()

# ------------------------------
# 2. PCA plot
# ------------------------------
# Perform PCA (samples as rows)
pca <- prcomp(t(exprs_data), scale. = TRUE)

# Base R plot
plot(pca$x[,1], pca$x[,2],
     col = ifelse(grepl("PCOS", colnames(exprs_data)), "red", "blue"),
     pch = 19,
     xlab = paste0("PC1 (", round(summary(pca)$importance[2,1]*100,1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca)$importance[2,2]*100,1), "%)"),
     main = "PCA of RMA Normalized Expression")
text(pca$x[,1], pca$x[,2], labels = colnames(exprs_data), pos = 3, cex = 0.8)

# Save PCA plot automatically
png("PCA_Normalized.png", width=1200, height=800)
plot(pca$x[,1], pca$x[,2],
     col = ifelse(grepl("PCOS", colnames(exprs_data)), "red", "blue"),
     pch = 19,
     xlab = paste0("PC1 (", round(summary(pca)$importance[2,1]*100,1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca)$importance[2,2]*100,1), "%)"),
     main = "PCA of RMA Normalized Expression")
text(pca$x[,1], pca$x[,2], labels = colnames(exprs_data), pos = 3, cex = 0.8)
dev.off()

