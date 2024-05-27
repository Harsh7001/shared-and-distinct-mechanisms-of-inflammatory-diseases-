# Load necessary libraries
library(WGCNA)
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

# Run the blockwiseConsensusModules function
consensusModules <- blockwiseConsensusModules(
  multiExpr = alignedMultiData,
  checkMissingData = TRUE,
  maxBlockSize = 30000,
  corType = "pearson",
  power = 1,  # Ensure this power is appropriate based on your sft1 and sft2 results
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


# Visualize the network of eigengenes
MEs_diss <- 1 - cor(consensusMEs)
METree <- hclust(as.dist(MEs_diss), method = "average")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Additional checks
# Make sure the moduleColors vector is correctly assigned and matches the genes
table(moduleColors)
