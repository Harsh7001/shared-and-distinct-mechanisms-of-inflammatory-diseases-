# Load devtools package
library(devtools)
library(ggplot2)

# Install the ReducedExperiment package from GitHub

BiocManager::install("clusterProfiler")
# Install missing dependencies
install.packages("biomaRt")

# Retry installing ReducedExperiment package from GitHub
devtools::install_github("jackgisby/ReducedExperiment")
data$tissue
library(ReducedExperiment)


IBD_me <- identify_modules(IBDdata, verbose=0, powers=1:30)
plotDendro(IBD_me)
IBD_me

reduced(IBD_me)[1:5, 1:5]
assignments(IBD_me)[1:5]
loadings(IBD_me)[1:5]



colData(IBDdata)
f <- "~ case_control + sex  + severity "
associations_me <- associate_components(IBD_me, method = "lm", formula = f)
levels(IBD_me$individual_id)
