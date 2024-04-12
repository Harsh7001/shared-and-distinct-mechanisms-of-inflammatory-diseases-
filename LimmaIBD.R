# Install and load necessary packages
install.packages("EnhancedVolcano")
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

library(limma)

# Assuming you already have a SummarizedExperiment object named "exp_data"
# Check the structure of the SummarizedExperiment object
print(IBDdata)

# Prepare experimental design (replace with your actual design)
design <- model.matrix(~ disease, data = colData(IBDdata))  # 'Group' is your experimental variable

# Get the count matrix from the SummarizedExperiment object
counts_matrix <- assay(IBDdata)

# Perform differential expression analysis using limma
fit <- lmFit(counts_matrix, design)
fit <- eBayes(fit)
results <- decideTests(fit)
results

# View top differentially expressed genes
topGenes <- topTable(fit, coef=1, number=10)  # Assuming you're interested in the first coefficient
print(topGenes)


colnames(colData(IBDdata))

# Volcano plot
volcanoPlot(fit, coef=1, highlight=10)  # Replace coef with the appropriate coefficient if needed
# Create volcano plot
# Create volcano plot
EnhancedVolcano(fit, contrast=1, lab=rownames(results), pCutoff=0.05, FCcutoff=1.5)
