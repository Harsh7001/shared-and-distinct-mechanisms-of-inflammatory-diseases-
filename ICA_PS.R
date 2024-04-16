# Install and load necessary packages
install.packages("fastICA")
library(fastICA)

# Assuming your SummarizedExperiment object is named "se"

# Extract count matrix from SummarizedExperiment object
counts <- assay(PSdata)

# Transpose the count matrix if necessary (ICA usually operates on samples as rows)
counts <- t(counts)

# Perform ICA
ica_result <- fastICA(counts, n.comp = 2)  # Adjust the number of components as needed

# Access the independent components (ICs)
independent_components <- ica_result$S

# Optionally, you can visualize the ICs or perform further analysis
