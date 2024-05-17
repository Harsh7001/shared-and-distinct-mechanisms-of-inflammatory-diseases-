
# Load the WGCNA package
library(WGCNA)
BiocManager::install("WGCNA")

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(SingleCellExperiment)

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
multiColor <- list(dataset1 = module_assignment_IBD, dataset2 = module_assignment_PS)

# Call the modulePreservation function with your data
modulePreservationResult <- modulePreservation(
  multiData = multiData,
  multiColor = multiColor,
  dataIsExpr = TRUE,
  nPermutations = 1
  
)



print(length(multiColor))
print(dim(multiColor[[1]]))

# Check the results
print(modulePreservation)

print(dim(counts1))
a <- list(module_assignment_IBD)
print(dim(counts2))
print(length(module_assignment_PS))


# Extract preserved modules based on predefined thresholds
preserved_modules <- rownames(modulePreservation)[modulePreservation$preservation & modulePreservation$Zsummary > 2]

# Extract genes in preserved modules
genes <- unlist(genesInModules(ps_modules)[[preserved_modules]])

# Enrichment analysis
enrichResult <- enrichGO(gene = genes, ...)

# Plot preservation
plotPreservation(modulePreservation)
