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


stability_resIBD <- estimate_stability(
  IBDdata, 
  seed=1, 
  n_runs=30, 
  min_components=10, max_components=60, by=2, 
  mean_stability_threshold=0.85
)

plot_stability(stability_resIBD, plot_path = "stabilityIBD.png")

IBD_fe <- estimate_factors(
  IBDdata, 
  nc=35, 
  seed=1, 
  use_stability=TRUE, 
  n_runs=30, 
  stability_threshold=0.25
)

IBD_fe

# get reduced data for first 5 samples and factors
reduced(PS_fe)[1:5, 1:5]


# get loadings for first 5 genes and factors
loadings(PS_fe)[1:5, 1:5]


IBD_me <- identify_modules(IBDdata, verbose=0, powers=1:30)
plotDendro(IBD_me)
IBD_me

reduced(IBD_me)[1:5, 1:5]
assignments(IBD_me)[1:5]
loadings(IBD_me)[1:5]



colData(IBDdata)
f <- "~ case_control + sex  + severity "
associations_meIBD <- associate_components(IBD_me, method = "lm", formula = f)
associations_feIBD <- associate_components(IBD_fe, method = "lm", formula = f1)
levels(IBD_me$individual_id)


associations_feIBD$summaries[associations_feIBD$summaries$term == "case_controlUC_I" ,]

ggplot(
  associations_feIBD$summaries[associations_feIBD$summaries$term == "case_controlUC_I" ,],
  aes(-log10(adj_pvalue), reorder(component, -log10(adj_pvalue)), fill = estimate)
) +
  geom_point(pch = 21, size = 3) +
  xlab("-log10(p-value)") +
  ylab("Factor name (ordered by p-value)") +
  geom_vline(xintercept = -log10(0.05)) +
  scale_fill_gradient2(low = "#3A3A98", high = "#832424")
##factors that are "upregulated"in IBD patients are shown in red, whereas "downregulated" factors are in blue. The vertical line indicates significance (i.e., adjusted p-valuesless than 0.05).

ggplot(
  associations_meIBD$summaries[associations_meIBD$summaries$term == "case_controlUC_I" ,],
  aes(-log10(adj_pvalue), reorder(component, -log10(adj_pvalue)), fill = estimate)
) +
  geom_point(pch = 21, size = 3) +
  xlab("-log10(p-value)") +
  ylab("Module name (ordered by p-value)") +
  geom_vline(xintercept = -log10(0.05)) +
  scale_fill_gradient2(low = "#3A3A98", high = "#832424")

aligned_featuresIBD <- getAlignedFeatures(IBD_fe, format = "data.frame")
head(aligned_featuresIBD)


factor_enrichIBD <- runEnrich(IBD_fe, as_dataframe = TRUE)
module_enrichIBD <- runEnrich(IBD_me, as_dataframe = TRUE)

ggplot(
  factor_enrich[factor_enrichIBD$component == "factor_28" & factor_enrichIBD$p.adjust < 0.05 ,][1:70, ],
  aes(-log10(p.adjust), reorder(substr(ID, 1, 45), -log10(p.adjust)), fill = enrichmentScore)
) +
  geom_point(pch = 21, size = 3) +
  xlab("-log10(p-value)") +
  ylab("Pathway name (ordered by p-value)") +
  geom_vline(xintercept = -log10(0.05)) +
  scale_fill_gradient2(low = "#3A3A98", high = "#832424")

ggplot(
  module_enrich[module_enrichIBD$component == "module_38" & module_enrichIBD$p.adjust < 0.05 ,][1:35, ],
  aes(-log10(p.adjust), reorder(substr(ID, 1, 45), -log10(p.adjust)))
) +
  geom_point(pch = 21, size = 3) +
  xlab("-log10(p-value)") +
  ylab("Pathway name (ordered by p-value)")
