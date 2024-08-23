
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
module_genes <- names(multiColor$dataset1[multiColor$dataset1 == "module_7"])

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




save.image(file='yoursession.RData')

# Assuming 'IBDdata' is your object and you want to export 'colData@listData'
metadata <- IBDdata@colData@listData
# Export metadata to CSV format
write.csv(metadata, file = "metadata.csv", row.names = FALSE)


assay <-assay(IBDdata)
write.csv(assay, file = "assay.csv", row.names = TRUE)

PS_assay <- assay(PSdata)
write.csv(PS_assay, file = "PS_assay.csv", row.names = TRUE)


# Assuming 'IBDdata' is your object and you want to export 'colData@listData'
PS_metadata <- PSdata@colData@listData
PS_metadata <- as.data.frame(PS_metadata)
# Export metadata to CSV format
write.csv(PS_metadata, file = "PS_metadata.csv", row.names = FALSE)


assay_ps <- assay(PSdata)
# PCA for UC vs Control

# Perform PCA
pca_ps <- prcomp(t(assay_ps), scale. = TRUE, center = TRUE)

# Prepare data for ggplot
pca_data_ps <- as.data.frame(pca_ps$x)
pca_data_ps$sample_id <- rownames(pca_data_ps)
pca_data_ps <- dplyr::left_join(pca_data_ps, PS_metadata, by = c("sample_id" = "sample_id"))

# Plot PCA
ggplot(pca_data_ps, aes(x = PC1, y = PC2, color = case_control)) +
  geom_point(size = 5) +
  labs(title = "PCA of PS Samples (FDR < 0.001)",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal()

# Perform PCA
pca_ibd <- prcomp(t(assay), scale. = TRUE, center = TRUE)

# Prepare data for ggplot
pca_data_ibd <- as.data.frame(pca_ibd$x)
pca_data_ibd$sample_id <- rownames(pca_data_ibd)
pca_data_ibd <- dplyr::left_join(pca_data_ibd, metadata, by = c("sample_id" = "sample_id"))

# Plot PCA
ggplot(pca_data_ibd, aes(x = PC1, y = PC2, color = case_control)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Control_NonI" = "grey", "CD_NonI" = "blue", "UC_NonI" = "blue", "CD_I" = "red", "UC_I" = "red")) +
  labs(title = "PCA of IBD Samples (FDR < 0.001)",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal()

# Prepare data for ggplot
pca_data_ibd <- as.data.frame(pca_ibd$x)
pca_data_ibd$sample_id <- rownames(pca_data_ibd)
pca_data_ibd <- dplyr::left_join(pca_data_ibd, metadata, by = c("sample_id" = "sample_id"))

# Create a new label column
pca_data_ibd$label <- dplyr::case_when(
  pca_data_ibd$case_control %in% c("CD_NonI", "UC_NonI") ~ "Non_Inflamed",
  pca_data_ibd$case_control %in% c("CD_I", "UC_I") ~ "Inflamed",
  pca_data_ibd$case_control == "Control_NonI" ~ "Control"
)

# Plot PCA with relabeled categories
ggplot(pca_data_ibd, aes(x = PC1, y = PC2, color = label)) +
  geom_point(size = 5) +
  scale_color_manual(values = c("Control" = "grey", "Non_Inflamed" = "blue", "Inflamed" = "red")) +
  labs(title = "PCA of IBD Samples (FDR < 0.001)",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal()

