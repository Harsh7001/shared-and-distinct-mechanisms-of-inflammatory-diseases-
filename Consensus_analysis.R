#### Consensus analysis

# Load necessary libraries
library(WGCNA)
library(VennDiagram)
library(grid)
options(stringsAsFactors = FALSE)

# Ensure multiData is defined
# multiData should be a list with two elements (IBD_me_matrix and PS_me_matrix), each containing a data frame with gene expression data
multiData <- list(
  IBD = list(data = IBD_me_matrix),
  PS = list(data = PS_me_matrix)
)

# Define soft-thresholding powers
powers <- c(1:10, seq(from = 12, to=30, by=2))

# Determine soft-thresholding power
sft1 <- pickSoftThreshold(t(IBD_me_matrix), powerVector = powers)
sft2 <- pickSoftThreshold(t(PS_me_matrix), powerVector = powers)





# Plotting for IBD dataset
par(mfrow = c(1, 2))
plot(sft1$fitIndices[, 1], -sign(sft1$fitIndices[, 3]) * sft1$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence - IBD")
text(sft1$fitIndices[, 1], -sign(sft1$fitIndices[, 3]) * sft1$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.90, col = "red")  # Threshold for scale-free topology

plot(sft1$fitIndices[, 1], sft1$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity - IBD")
text(sft1$fitIndices[, 1], sft1$fitIndices[, 5],
     labels = powers, col = "red")

# Plotting for PS dataset
par(mfrow = c(1, 2))
plot(sft2$fitIndices[, 1], -sign(sft2$fitIndices[, 3]) * sft2$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence - PS")
text(sft2$fitIndices[, 1], -sign(sft2$fitIndices[, 3]) * sft2$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.90, col = "red")  # Threshold for scale-free topology

plot(sft2$fitIndices[, 1], sft2$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity - PS")
text(sft2$fitIndices[, 1], sft2$fitIndices[, 5],
     labels = powers, col = "red")








# Function to align genes across datasets
alignGenes <- function(multiExpr) {
  commonGenes <- Reduce(intersect, lapply(multiExpr, function(set) colnames(set$data)))
  alignedSets <- lapply(multiExpr, function(set) {
    set$data <- set$data[, commonGenes, drop = FALSE]
    set
  })
  return(alignedSets)
}

# Align genes in multiData
alignedMultiData <- alignGenes(multiData)



# Check dimensions of aligned data
print(sapply(alignedMultiData, function(set) dim(set$data)))



# Create a Venn diagram
genes_IBD <- colnames(IBD_me_matrix)
genes_PS <- colnames(PS_me_matrix)
venn.plot <- venn.diagram(
  x = list(IBD = genes_IBD, PS = genes_PS),
  category.names = c("IBD", "PS"),
  filename = NULL,
  output = TRUE
)
grid.draw(venn.plot)

common_genes <- intersect(genes_IBD, genes_PS)
length(common_genes)

# Print the common genes
print(common_genes)


# Run the blockwiseConsensusModules function
consensusModules <- blockwiseConsensusModules(
  multiExpr = alignedMultiData,
  checkMissingData = TRUE,
  maxBlockSize = 30000,
  corType = "pearson",
  power = 6,  # Ensure this power is appropriate based on your sft1 and sft2 results
  networkType = "unsigned",
  TOMType = "unsigned",
  saveIndividualTOMs = TRUE,
  saveConsensusTOMs = FALSE,
  verbose = 3
)

# Examine results
moduleColors <- consensusModules$colors
moduleEigengenes <- consensusModules$multiMEs

# Ensure the length of the module colors vector matches the number of genes in the dendrogram
print(length(moduleColors))
print(length(consensusModules$dendrograms[[1]]$order))

# Plot dendrogram of one block with module colors
plotDendroAndColors(consensusModules$dendrograms[[1]], moduleColors, 
                    "Module colors", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05)

# Inspect the structure of the consensusModules object
str(consensusModules)

# Print a summary of the consensus modules
print(consensusModules)

# Calculate and inspect module eigengenes
consensusMEs <- consensusModules$multiMEs
head(consensusMEs)
## Inspect the structure of consensusMEs
str(consensusMEs)

# Ensure all values in consensusMEs are numeric
consensusMEs <- as.data.frame(lapply(consensusMEs, as.numeric))




# Extract the data frames from consensusMEs
data1 <- consensusMEs$dataset1$data
data2 <- consensusMEs$dataset2$data

# Combine the data frames
combined_data <- rbind(data1, data2)

# Ensure all values are numeric
combined_data <- as.data.frame(lapply(combined_data, as.numeric))

# Calculate the dissimilarity matrix
MEs_diss <- 1 - cor(combined_data)

# Perform hierarchical clustering
METree <- hclust(as.dist(MEs_diss), method = "average")

# Plot the dendrogram
plot(METree, main = "Clustering of Module Eigengenes", sub = "", xlab = "", cex.main = 2)




# Additional checks
# Make sure the moduleColors vector is correctly assigned and matches the genes
table(moduleColors)

# Assuming you have run blockwiseConsensusModules and have module colors
moduleColors_IBD <- consensusModules$colors[[1]]
moduleColors_PS <- consensusModules$colors[[2]]

# Create a table to show overlap of module memberships
table(moduleColors_IBD, moduleColors_PS)


###################





###
# Load the WGCNA library
library(WGCNA)

# Assuming counts1, counts2, module_assignment_IBD, module_assignment_PS,
# IBD_me_matrix, and PS_me_matrix are already defined in your environment.

# Transpose the count matrices if needed
dat1 <- t(counts1)
dat2 <- t(counts2)

# Transpose the module eigengene matrices
MEs1 <- t(IBD_me_matrix)
MEs2 <- t(PS_me_matrix)

# Define dataset names
name1 <- "IBD"
name2 <- "PS"

# Run the overlapTableUsingKME function
overlap_table <- overlapTableUsingKME(
  dat1 = dat1,
  dat2 = dat2,
  colorh1 = module_assignment_IBD,
  colorh2 = module_assignment_PS,
  name1 = name1,
  name2 = name2,
  cutoffMethod = "assigned",
  datIsExpression = TRUE
)

# Print the resulting overlap table
print(overlap_table)
print(dim(dat1))
print(dim(dat2))
print(length(module_assignment_IBD))
print(length(module_assignment_PS))
print(dim(MEs1))
print(dim(MEs2))






#############

overlap_matrix <- matrix(NA, nrow = moduleNames(IBD_me), ncol = moduleNames(PS_me))
overlap_matrix <- matrix(0, nrow = length(moduleNames(IBD_me)), ncol = length(moduleNames(PS_me)))
rownames(overlap_matrix) <- moduleNames(IBD_me)
colnames(overlap_matrix) <- moduleNames(PS_me)

for (m_ibd in moduleNames(IBD_me)) {
  for (m_ps in moduleNames(PS_me)) {
    # list of genes in the IBD module
    module_genes_ibd <- module_assignment_IBD[names(module_assignment_IBD) == m_ibd]
    # list of genes in the PS module
    module_genes_ps <- module_assignment_PS[names(module_assignment_PS) == m_ps]
    
    # Overlap = number of overlapping genes / size of the smaller module
    overlap <- sum(module_genes_ibd %in% module_genes_ps) / min(length(module_genes_ibd), length(module_genes_ps))
    
    # Store the overlap in the matrix
    overlap_matrix[which(moduleNames(IBD_me) == m_ibd), which(moduleNames(PS_me) == m_ps)] <- overlap
  }
}

# Print overlap matrix
print(overlap_matrix)





library(ggplot2)
library(reshape2) # for melt function

overlap_matrix_melted <- melt(overlap_matrix)
names(overlap_matrix_melted) <- c("IBD_Module", "PS_Module", "Overlap")
ggplot(overlap_matrix_melted, aes(x = IBD_Module, y = PS_Module, fill = Overlap)) +
  geom_tile() +
  scale_fill_gradient(low = "skyblue", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill = "Overlap Score", title = "Overlap of Genes between IBD and PS Modules")



# Find the maximum value in the overlap matrix
max_overlap_value <- max(overlap_matrix, na.rm = TRUE)

# Identify the position of the maximum value
max_position <- which(overlap_matrix == max_overlap_value, arr.ind = TRUE)

# Get the module names corresponding to the maximum overlap
highest_overlapping_modules <- data.frame(
  IBD_Module = rownames(overlap_matrix)[max_position[, 1]],
  PS_Module = colnames(overlap_matrix)[max_position[, 2]],
  Overlap_Value = rep(max_overlap_value, nrow(max_position))
)

# Print the highest overlapping modules and their overlap value
print(highest_overlapping_modules)


# Number of top overlaps you want to retrieve
top_n = 15  # Change this number based on how many top overlaps you want

# Flatten the matrix to a dataframe and include module names
overlap_data <- as.data.frame(as.table(overlap_matrix))
names(overlap_data) <- c("IBD_Module", "PS_Module", "Overlap_Value")

# Filter out NA values if any
overlap_data <- subset(overlap_data, !is.na(Overlap_Value))

# Sort the data to find the top overlaps
top_overlaps <- overlap_data[order(-overlap_data$Overlap_Value), ]

# Retrieve the top n overlaps
top_overlaps <- head(top_overlaps, top_n)

# Print the top overlapping modules
print(top_overlaps)















