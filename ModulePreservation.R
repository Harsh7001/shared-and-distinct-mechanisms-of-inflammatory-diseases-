
# Load the WGCNA package
library(WGCNA)
BiocManager::install("WGCNA")
BiocManager::install("org.Hs.eg.db")
install.packages("gtable")

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(SingleCellExperiment)
if (!requireNamespace("enrichplot", quietly = TRUE)) {
  install.packages("enrichplot")
}
library(enrichplot)

# network analysis & visualization package:
library(igraph)
# PS data
PS_me
# IBD data  
IBD_me

##########

# Convert your ModularExperiment objects to expression matrices

IBD_me_matrix <- as.matrix(IBD_me)
PS_me_matrix <- as.matrix(PS_me)

IBD_me_matrix <- assay(IBD_me, "transformed")
PS_me_matrix <- assay(PS_me, "transformed")

#########################
eigengene_PS
eigengene_IBD

# Check the class and dimensions of the input data
class(eigengene_IBD)
dim(eigengene_IBD)

class(eigengene_PS)
dim(eigengene_PS)

#########################
# Assuming your data is organized as described above
# Create lists containing your data
multiData <- list(dataset1 = list(data = t(IBD_me_matrix)), dataset2 = list(data = t(PS_me_matrix)))

# Create lists containing module assignments for each dataset
multiColor <- list(dataset1 = setNames(names(module_assignment_IBD), module_assignment_IBD), dataset2 = setNames(names(module_assignment_PS), module_assignment_PS))
# Call the modulePreservation function with your data
modulePreservationResult <- modulePreservation(
  multiData = multiData,
  multiColor = multiColor,
  dataIsExpr = TRUE,
  nPermutations = 100,
  quickCor = 0,
  verbose= 50
)

# Check the results
print(modulePreservationResult)



# Extract the preservation statistics
preservationStats <- modulePreservationResult$preservation$Z$ref.dataset1$inColumnsAlsoPresentIn.dataset2

# Convert the preservation statistics to a data frame for easier handling
preservationStats_df <- as.data.frame(preservationStats)

# Print the preservation statistics for inspection
print(preservationStats_df)

# Plot the preservation statistics
# A common plot is the Zsummary vs. module size
moduleSizes <- preservationStats_df$moduleSize
Zsummary <- preservationStats_df$Zsummary.pres

# Create a basic scatter plot
plot(moduleSizes, Zsummary, 
     xlab = "Module Size", 
     ylab = "Zsummary Preservation",
     main = "Module Preservation Statistics",
     pch = 19)

# Add a horizontal line at Zsummary = 2 and 10 for interpretation
abline(h = 2, col = "blue", lty = 2)
abline(h = 10, col = "red", lty = 2)

# Add text labels to the plot
text(moduleSizes, Zsummary, labels = rownames(preservationStats_df), pos = 4, cex = 0.7)

# You may also want to plot the medianRank statistics
medianRank <- preservationStats_df$medianRank.pres

# Create a basic scatter plot for medianRank
plot(moduleSizes, medianRank, 
     xlab = "Module Size", 
     ylab = "Median Rank Preservation",
     main = "Median Rank Preservation Statistics",
     pch = 19)

# Add text labels to the plot
text(moduleSizes, medianRank, labels = rownames(preservationStats_df), pos = 4, cex = 0.7)
dev.list()
# Interpretation:
# Higher Zsummary values (above 10) indicate strong evidence of preservation.
# Median Rank is another measure where lower values indicate better preservation.

# Save the module preservation results to a file for later use
write.csv(preservationStats_df, file = "module_preservation_stats.csv", row.names = TRUE)





# Filter for strongly preserved modules
strongly_preserved <- preservationStats_df[preservationStats_df$Zsummary.pres > 5, ]
print(strongly_preserved)

# Filter for weakly or non-preserved modules
weakly_preserved <- preservationStats_df[preservationStats_df$Zsummary.pres < 5, ]
print(weakly_preserved)


#####################################
# Load necessary libraries for enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)

# Example: Perform GO enrichment analysis for a specific module
# Extract genes belonging to "module_1"
module_genes <- names(multiColor$dataset1[multiColor$dataset1 == "module_1"])

# Verify the genes extracted
print(module_genes)

ego <- enrichGO(gene = module_genes, 
                OrgDb = org.Hs.eg.db, 
                keyType = "ENSEMBL", 
                ont = "BP", 
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

# View the GO enrichment results
head(ego)
# Compute term similarity matrix
ego <- pairwise_termsim(ego)
summary(ego)

barplot(ego, showCategory = 20, title = "GO Enrichment Analysis - Bar Plot")
dotplot(ego, showCategory = 20, title = "GO Enrichment Analysis - Dot Plot")
emapplot(ego, showCategory = 50, title = "GO Enrichment Analysis - Enrichment Map")
cnetplot(ego, categorySize = "pvalue", foldChange = NULL, showCategory = 5, title = "GO Enrichment Analysis - Network Plot")

#####################################

# Extract the unique modules
modules <- unique(multiColor$dataset1)

# Loop through each module to create and plot the network
for (module in modules) {
  
  # Get the genes in the current module
  module_genes <- names(multiColor$dataset1[multiColor$dataset1 == "module_32"])
  
  # Ensure the gene names are present in the row names of IBD_me_matrix
  valid_genes <- module_genes[module_genes %in% rownames(IBD_me_matrix)]
  
  # Check if valid_genes has at least 2 genes
  if (length(valid_genes) < 2) {
    next # Skip if not enough genes to form a network
  }
  
  # Extract the expression data for these genes
  module_expr_matrix <- IBD_me_matrix[valid_genes, ]
  
  # Compute the correlation matrix
  cor_matrix <- cor(t(module_expr_matrix))
  
  # Apply a threshold to create an adjacency matrix
  threshold <- 0.7
  adj_matrix <- ifelse(abs(cor_matrix) > threshold, 1, 0)
  
  # Convert the adjacency matrix to a graph
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
  
  # Plot the network
  plot(g, vertex.label = valid_genes, vertex.size = 5, edge.arrow.size = 0.5, main = paste("Network for", module))
}




save.image(file='yoursession.RData')
