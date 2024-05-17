counts1 <- assay(IBDdata)
counts2 <- assay(PSdata)


metadata(rowData(se))
head(rownames(IBDdata))
dplyr::glimpse(IBDdata)
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
library(ggplot2)

# Install the ReducedExperiment package from GitHub
   
BiocManager::install("clusterProfiler")
# Install missing dependencies
install.packages("biomaRt")

# Retry installing ReducedExperiment package from GitHub
devtools::install_github("jackgisby/ReducedExperiment")
data$tissue
library(ReducedExperiment)


stability_res <- estimate_stability(
  PSdata, 
  seed=1, 
  n_runs=30, 
  min_components=10, max_components=60, by=2, 
  mean_stability_threshold=0.85
)

plot_stability(stability_res, plot_path = "stabilityPS.png")

PS_fe <- estimate_factors(
  PSdata, 
  nc=35, 
  seed=1, 
  use_stability=TRUE, 
  n_runs=30, 
  stability_threshold=0.25
)

PS_fe

# get reduced data for first 5 samples and factors
reduced(PS_fe)[1:5, 1:5]




# get loadings for first 5 genes and factors
loadings(PS_fe)[1:5, 1:5]

PS_me <- identify_modules(PSdata, verbose=0, powers=1:30)
plotDendro(PS_me)
PS_me

reduced(PS_me)[1:5, 1:5]
assignments(PS_me)[1:5]
loadings(PS_me)[1:5]

eigengene_PS <- reduced(PS_me)
module_assignment_PS <- assignments(PS_me)


colData(PSdata)
f1 <- "~ case_control"
associations_mePS <- associate_components(PS_me, method = "lm", formula = f1)
associations_fePS <- associate_components(PS_fe, method = "lm", formula = f1)
levels(PS_me$sra_accession)


associations_fePS$summaries[associations_fePS$summaries$term == "case_controlPS" ,]

ggplot(
  associations_fePS$summaries[associations_fePS$summaries$term == "case_controlPS" ,],
  aes(-log10(adj_pvalue), reorder(component, -log10(adj_pvalue)), fill = estimate)
) +
  geom_point(pch = 21, size = 3) +
  xlab("-log10(p-value)") +
  ylab("Factor name (ordered by p-value)") +
  geom_vline(xintercept = -log10(0.05)) +
  scale_fill_gradient2(low = "#3A3A98", high = "#832424")
##factors that are "upregulated"in PS patients are shown in red, whereas "downregulated" factors are in blue. The vertical line indicates significance (i.e., adjusted p-valuesless than 0.05).

ggplot(
  associations_mePS$summaries[associations_mePS$summaries$term == "case_controlPS" ,],
  aes(-log10(adj_pvalue), reorder(component, -log10(adj_pvalue)), fill = estimate)
) +
  geom_point(pch = 21, size = 3) +
  xlab("-log10(p-value)") +
  ylab("Module name (ordered by p-value)") +
  geom_vline(xintercept = -log10(0.05)) +
  scale_fill_gradient2(low = "#3A3A98", high = "#832424")

aligned_features <- getAlignedFeatures(PS_fe, format = "data.frame")
head(aligned_features)


factor_enrich <- runEnrich(PS_fe, as_dataframe = TRUE)
module_enrich <- runEnrich(PS_me, as_dataframe = TRUE)

ggplot(
  factor_enrich[factor_enrich$component == "factor_17" & factor_enrich$p.adjust < 0.05 ,][1:35, ],
  aes(-log10(p.adjust), reorder(substr(ID, 1, 45), -log10(p.adjust)), fill = enrichmentScore)
) +
  geom_point(pch = 21, size = 3) +
  xlab("-log10(p-value)") +
  ylab("Pathway name (ordered by p-value)") +
  geom_vline(xintercept = -log10(0.05)) +
  scale_fill_gradient2(low = "#3A3A98", high = "#832424")

ggplot(
  module_enrich[module_enrich$component == "module_5" & module_enrich$p.adjust < 0.05 ,][1:35, ],
  aes(-log10(p.adjust), reorder(substr(ID, 1, 45), -log10(p.adjust)))
) +
  geom_point(pch = 21, size = 3) +
  xlab("-log10(p-value)") +
  ylab("Pathway name (ordered by p-value)")

