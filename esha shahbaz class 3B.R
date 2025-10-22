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


