# Install and load necessary packages
install.packages("fastICA")
library(fastICA)

# Assuming your SummarizedExperiment object is named "se"

# Extract count matrix from SummarizedExperiment object
IBDcounts <- assay(IBDdata)

# Transpose the count matrix if necessary (ICA usually operates on samples as rows)
#IBDcounts <- t(IBDcounts)

# Perform ICA
IBD_ica_result <- fastICA(IBDcounts, n.comp = 8)  # Adjust the number of components as needed

# Access the independent components (ICs)
IBD_independent_components <- IBD_ica_result$S

# Optionally, you can visualize the ICs or perform further analysis
pheatmap(IBD_independent_components, scale = "row", cluster_rows = TRUE)

dim(IBDindependent_components)