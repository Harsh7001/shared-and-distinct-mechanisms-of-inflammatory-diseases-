Project title: Investigating shared and distinct mechanisms of inflammatory diseases 

About the project: While transcriptomics has seen considerable use in investigating inflammatory diseases, most studies of gene expression consider only a single disease in isolation. Relatively few studies have attempted to investigate shared signatures of inflammatory diseases, in part due to the challenges involved in integrating expression across different experiments. We gained access to data from the IMID-Bio-UK network (https://www.gla.ac.uk/research/az/imid/aboutus/) and curated public datasets. This study will investigate the transcriptomes of inflammatory disease both in bulk RNA-seq and single-cell data. Signatures of disease and drug response will be identified by applying dimensionality reduction methods, such as weighted gene correlation network analysis (WGCNA) [1] and independent component analysis (ICA) [2]. The resulting signatures will then be compared across multiple datasets to differentiate between shared and distinct mechanisms of inflammatory disorders that might be targeted by new or repositioned therapies. Machine learning approaches could also be applied to compare and integrate these data; for instance, we could alternatively apply autoencoders to identify latent variables and compare diseases.

 

Plan 

Core project: comparing signatures of IMIDs 

    Data ingestion and exploration of IMID data (likely looking at Atopic Dermatitis in addition to existing Psoriasis and IBD datasets) 

    Comparing gene-level signatures of disease 

    Using dimensionality reduction to simplify and interpret IMID data 

    Applying dimensionality reduction across datasets to compare IMID signatures 

Potential follow-up work 

    Apply and compare alternative approaches for interpretable dimensionality reduction (e.g., regularised factor analysis, deep learning) 

    Predict disease outcomes (using projected signatures) 

    Compare IMID signatures to signatures of drug response (using CMap/L1000 data) 

    Add functionality to R package for applying dimensionality reduction 

 

Initial steps 

    Reading into relevant methods/concepts (RNA-seq, WGCNA/ICA, etc.) 

    Started with a couple of datasets (Inflammatory bowel disease and Psoriasis) 

    Had a look at the studies linked to the data 

    Performed some exploratory analysis (e.g., PCA) 

    Batch effect correction:  

    1. Combat- It required prerequisite knowledge of batches of data collection, which is only available for partial dataset, therefore, it didn't pan out. 

    2.SVA- It gives us an estimate of batches based on the counts data, but due to normalized counts it seems to bug out and not give accurate results. 

    We concluded to do the analysis without doing batch corrections as most of the data used should have accounted for batches by respective labs. 

    Performed some standard analyses (e.g., differential expression) 

    1. Deseq2 - did not work since our dataset contains normalized counts and Deseq2 requires raw counts. 

    2. limma -  

    Compare results to those published in the original studies 

    Then move on to apply dimensionality reduction 

    Started by applying WGCNA – could start with ReducedExperiment package, but worthwhile to understand and/or attempt to apply underlying WGCNA package 

    Attempt to characterise and interpret the resultant modules in each dataset 

    Compare modules across datasets 

    ICA - 

    EdgeR  

 

 

