# Load necessary libraries
library(readr)   # for reading RDS files
library(dplyr)   # for data manipulation
library(ggplot2) # for visualization
library(DESeq2)
library(SummarizedExperiment)
library(edgeR)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")

browseVignettes("SummarizedExperiment")
# Load the RDS file
IBDdata <- readRDS("IBDse.rds")
PSdata <- readRDS("PSse.rds")


# Display the structure of the data
str(IBDdata)
dim(IBDdata)
dim(PSdata)

IBDdata

colData(IBDdata)
colData(PSdata)

# Summary statistics
summary(IBDdata)
summary(PSdata)

hist(assay(IBDdata)[,1], main = "Distribution of Counts in File 1", xlab = "Counts")
hist(assay(PSdata)[,1], main = "Distribution of Counts in File 2", xlab = "Counts")

head(IBDdata)
if (any(assay(IBDdata) < 0)) {
  stop("Some values in the assay matrix are negative.")
}
dds1 <- DESeqDataSet(assay(IBDdata), design = ~ case_control)



# Load necessary libraries
library(WGCNA)
library(SummarizedExperiment)

# Assuming you have loaded your SummarizedExperiment objects se1 and se2

# Step 1: Data Preparation
# Extract count data from SummarizedExperiment objects
counts1 <- assay(IBDdata)
counts2 <- assay(PSdata)
metadata1 <- metadata(IBDdata)
metadata <- colData(IBDdata)
# Step 2: Data Preprocessing
# You may need to perform normalization and transformation steps here if necessary
# For RNA-seq data, commonly used methods include variance stabilizing transformation (VST) or regularized logarithm transformation (rlog)
counts1 <- vst(counts1, blind=FALSE)
# Step 3: Constructing a Co-expression Network
# Combine count data from both experiments
all_counts <- cbind(counts1, counts2)

# Calculate correlations
correlation_matrix <- cor(counts1)
correlation_matrix <- t(correlation_matrix)
# Construct adjacency matrix using soft-thresholding
powers <- pickSoftThreshold(correlation_matrix)

# Build the network
adjacency_matrix <- adjacency(correlation_matrix, power = powers$power)

# Step 4: Module Detection
# Perform hierarchical clustering
geneTree <- hclust(as.dist(1 - correlation_matrix), method = "average")

# Module detection using dynamic tree cut
modules <- cutreeDynamic(dendro = geneTree, distM = as.matrix(1 - correlation_matrix), deepSplit = 2, minClusterSize = 30)

# Step 5: Module-Trait Relationships
# If you have sample trait information, you can explore module-trait relationships using correlation analysis
# For example, you can correlate module eigengenes with sample traits

# Step 6: Visualization and Interpretation
# Visualize the network, modules, and module-trait relationships
# You can use functions like plotDendroAndColors(), moduleEigengenes(), etc.

# Example plots
plotDendroAndColors(geneTree, modules)
moduleEigengenes <- moduleEigengenes(counts1, colors = modules)$eigengenes

plotEigengeneNetworks(moduleEigengenes, adjacency_matrix)




modules





par(mfrow = c(100, 100), oma = c(2, 2, 2, 2))  # Adjust margins as needed
plotEigengeneNetworks(moduleEigengenes, adjacency_matrix)
par(mfrow = c(1, 1))  # Reset margins to default









# Visualizations
# Example: Histogram of a numeric variable
ggplot(data, aes(x = numeric_variable)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Numeric Variable")

# Example: Box plot of a numeric variable
ggplot(data, aes(x = factor(categorical_variable), y = numeric_variable)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  labs(title = "Box Plot of Numeric Variable by Category")

# Example: Scatter plot of two numeric variables
ggplot(data, aes(x = numeric_variable1, y = numeric_variable2)) +
  geom_point(color = "blue") +
  labs(title = "Scatter Plot of Two Numeric Variables")

# Handle missing values
# Example: Remove rows with missing values
data_clean <- na.omit(data)

# Outlier detection
# Example: Box plot to identify outliers
ggplot(data, aes(x = "", y = numeric_variable)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  labs(title = "Box Plot of Numeric Variable")

# Correlation matrix
correlation_matrix <- cor(data_clean[, c("numeric_variable1", "numeric_variable2")])
print(correlation_matrix)
