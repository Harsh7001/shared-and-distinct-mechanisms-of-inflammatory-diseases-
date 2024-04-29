# Install and load necessary packages
install.packages("fastICA")
library(fastICA)

# Assuming your SummarizedExperiment object is named "se"

# Extract count matrix from SummarizedExperiment object
countsIBD <- assay(IBDdata)

# Transpose the count matrix if necessary (ICA usually operates on samples as rows)
#countsIBD <- t(countsIBD)

# Perform ICA
ica_result <- fastICA(countsIBD, n.comp = 8)  # Adjusting the number of components as needed

# Access the independent components (ICs)
independent_components <- ica_result$S

# Optionally, you can visualize the ICs or perform further analysis
# Example heatmap of the first independent component
library(pheatmap)
pheatmap(independent_components, scale = "row", cluster_rows = TRUE)

dim(independent_components)
par(mfcol = c(2,2))
# Plot original gene expression profiles
plot(countsIBD[,1], type="l", main="Original Gene 1") 
plot(countsIBD[,2], type="l", main="Original Gene 2")

# Plot ICA components  
plot(ica_result$S[,1], type="l", main="ICA Component 1")
plot(ica_result$S[,2], type="l", main="ICA Component 2")
#---------------------------------------------------
# 1: un-mixing two mixed independent uniforms
#---------------------------------------------------

a <- fastICA(countsIBD, 2, alg.typ = "parallel", fun = "logcosh", alpha = 1, 
             method = "C", row.norm = FALSE, maxit = 200, 
             tol = 0.0001, verbose = TRUE)
par(mfrow = c(1, 3))
plot(a$X, main = "Pre-processed data")
plot(a$X %*% a$K, main = "PCA components")
plot(a$S, main = "ICA components")


#--------------------------------------------
# 2: un-mixing two independent signals
#--------------------------------------------

# Generate mixed signals
S <- cbind(sin((1:1000)/20), rep((((1:200)-100)/100), 5))
A <- matrix(c(0.291, 0.6557, -0.5439, 0.5572), 2, 2)
X <- S %*% A

# Perform ICA using the count data
a <- fastICA(countsIBD, 2, alg.typ = "parallel", fun = "logcosh", alpha = 1, 
             method = "R", row.norm = FALSE, maxit = 200, 
             tol = 0.0001, verbose = TRUE)

par(mfcol = c(2, 3))

# Plot original signals (gene expression levels)
plot(1:nrow(countsIBD), countsIBD[, 1], type = "l", main = "Original Signal 1", 
     xlab = "", ylab = "")
plot(1:nrow(countsIBD), countsIBD[, 2], type = "l", main = "Original Signal 2", 
     xlab = "", ylab = "")

# Plot mixed signals
plot(1:1000, X[, 1], type = "l", main = "Mixed Signal 1", 
     xlab = "", ylab = "")
plot(1:1000, X[, 2], type = "l", main = "Mixed Signal 2", 
     xlab = "", ylab = "")

# Plot ICA source estimates
plot(1:nrow(countsIBD), a$S[, 1], type = "l", main = "ICA Source Estimate 1", 
     xlab = "", ylab = "")
plot(1:nrow(countsIBD), a$S[, 2], type = "l", main = "ICA Source Estimate 2", 
     xlab = "", ylab = "")









# Load packages
library(fastICA) 
library(SummarizedExperiment)
library(pheatmap)

# Load packages and data
library(fastICA)
library(SummarizedExperiment) 


# Select two genes/features to mix  
gene1 <- countsIBD[,1]
gene2 <- countsIBD[,2] 

# Mix the selected genes
mixing_matrix <- matrix(c(0.7, 0.3, -0.1, 0.9), 2, 2)
mixed_genes <- cbind(gene1, gene2) %*% mixing_matrix

# Perform ICA on mixed genes
ica_result <- fastICA(mixed_genes, 2) 

par(mfrow=c(2,2))

# Plot original genes  
plot(gene1, main="Original Gene 1")
plot(gene2, main="Original Gene 2")

# Plot mixed genes
plot(mixed_genes, main="Mixed Genes") 

# Plot estimated sources from ICA
plot(ica_result$S, main="ICA Estimates")

plot(ica_sim$S, main="ICA Estimates")


