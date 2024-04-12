metadata(rowData(se))
head(rownames(se))
dplyr::glimpse(se)
dim(se)
dimnames(se)
se[1:2, 2:1]
colData(se)
se$extract_protocol
sessionInfo()

# Install devtools package if not already installed
install.packages("devtools")

# Load devtools package
library(devtools)

# Install the ReducedExperiment package from GitHub
   
BiocManager::install("biomaRt")
# Install missing dependencies
install.packages("biomaRt")

# Retry installing ReducedExperiment package from GitHub
devtools::install_github("jackgisby/ReducedExperiment")
data$tissue
library(ReducedExperiment)


a <- ReducedExperiment(se)
