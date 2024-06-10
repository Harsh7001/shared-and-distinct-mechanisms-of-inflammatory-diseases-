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


devtools::install_github("jackgisby/ReducedExperiment")


###################
centrality <- getCentrality(object = multiData)
###################
# Example definition of the overlap function
overlap <- function(m1, m2) {
  # Assuming m1 and m2 are vectors/matrices and we want to find overlapping elements
  common_elements <- intersect(m1, m2)
  print(common_elements)
  return(common_elements)
}

# Example usage with your loop
for (m1 in PS_me_matrix) {
  for (m2 in IBD_me_matrix) {
    overlap(m1, m2)
  }
}

# Print sample data from PS_me_matrix and IBD_me_matrix
print(head(PS_me))
print(head(IBD_me))
# Define a simple overlap function to find common elements
overlap <- function(m1, m2) {
  # Assuming m1 and m2 are vectors/matrices and we want to find overlapping elements
  common_elements <- intersect(m1, m2)
  print(common_elements)
  return(common_elements)
}
for (i in 1:nrow(PS_me)) {
  for (j in 1:nrow(IBD_me)) {
    overlap(PS_me[i, ], IBD_me[j, ])
  }
}
for (m1 in module_assignment_IBD) {
  for (m2 in module_assignment_PS) {
    overlap(m1, m2)
  }
}


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






###########

IBD_me
PS_me

for (m_ibd in moduleNames(IBD_me)) {
  for (m_ps in moduleNames(PS_me)) {
    # list of genes in the ibd module
    module_genes_ibd <- which(module_assignment_IBD == m_ibd)
    # list of genes in the PS module
    module_genes_ps <- which(module_assignment_PS == m_ps)
    # Overlap = number of overlapping genes / size of the smaller module
    overlap <- sum(module_genes_ibd %in% module_genes_ps) / min(length(module_genes_ibd), length(module_genes_ps))
    print(overlap)
  }
}


overlap_matrix # ibd modules on x, ps modules on y


#############

overlap_matrix <- matrix(NA, nrow = moduleNames(IBD_me), ncol = moduleNames(PS_me)))

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

























