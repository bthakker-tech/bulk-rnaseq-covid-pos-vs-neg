Bulk RNA-seq Analysis of SARS-CoV-2 Positive vs Negative Samples
(Python preprocessing + R/DESeq2 statistical analysis)

Project Overview

This project performs an end-to-end bulk RNA-seq analysis to characterize host transcriptional responses associated with SARS-CoV-2 infection. Gene expression profiles from SARS-CoV-2 positive (POS) and negative (NEG) human samples were analyzed to identify differentially expressed genes and enriched biological pathways.

The workflow follows a hybrid Python → R pipeline:

Python was used for data inspection, quality control, filtering, normalization for visualization, exploratory PCA, and metadata parsing.
R (Bioconductor) was used for statistical modeling with DESeq2 and pathway-level interpretation using Gene Set Enrichment Analysis (GSEA).

Dataset Source

The dataset analyzed in this project was obtained from the NCBI Gene Expression Omnibus (GEO).
GEO accession: GSE152075

Organism: Human

Data type: Bulk RNA-seq

Study focus: Host transcriptional response to SARS-CoV-2 infection

🔗 GEO link:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152075

Raw gene-level count data and the GEO series matrix file were downloaded directly from GEO.

Project Structure
bulk-rnaseq-covid-pos-vs-neg/
├── data/
│   └── raw/
│       ├── counts_filtered.csv
│       └── metadata.csv
├── scripts/
│   ├── python/
│   │   └── 01_preprocessing_qc.py
│   └── r/
│       └── 02_deseq2_gsea_pipeline.R
├── results/
│   ├── python_qc/
│   │   ├── pca_by_positivity.png
│   │   └── pca_by_batch.png
│   ├── deseq2/
│   │   ├── pos_vs_neg_unadjusted_DEG.csv
│   │   └── pos_vs_neg_adj_gender_age_DEG.csv
│   ├── heatmaps/
│   └── gsea/
└── README.md

Analysis Objectives

Perform quality control and exploratory analysis of bulk RNA-seq data
Identify differentially expressed genes (DEGs) between SARS-CoV-2 POS and NEG samples
Assess robustness of results after adjusting for gender and age
Identify biological pathways upregulated and downregulated in infected samples using GSEA
Produce a fully reproducible, documented analysis pipeline

Methods Summary
1. Data Preprocessing and Quality Control (Python)
Initial data handling and exploratory analysis were performed in Python.

Steps performed:

Loaded raw gene-level count matrix from GEO
Verified data integrity and checked for missing values
Removed genes with zero counts across all samples
Filtered lowly expressed genes (genes with >1 count in fewer than 5 samples removed)
Performed log2-CPM normalization for visualization only
Conducted principal component analysis (PCA) to examine global expression patterns
Exploratory PCA analyses included:
PCA colored by SARS-CoV-2 positivity status
PCA colored by sequencing batch
These analyses were used to assess sample separation and potential batch effects prior to differential expression analysis.
Metadata (positivity status, Ct value, age, gender, batch) was parsed directly from the GEO series matrix file, cleaned, and aligned with the expression data.

The processed outputs were exported as:

counts_filtered.csv
metadata.csv

2. Differential Expression Analysis (R / DESeq2)

Differential expression analysis was performed in R using DESeq2.
Models evaluated:

Unadjusted model:
~ positivity

Adjusted model:
~ gender + age + positivity

Key details:

Raw counts (not normalized values) were used as required by DESeq2
POS vs NEG contrasts were explicitly specified
Genes with missing adjusted p-values were removed
Differentially expressed genes (DEGs) were defined as:
Adjusted p-value < 0.05
|log2 fold change| ≥ 1
A comparison between unadjusted and adjusted models was performed to assess the stability of DEGs after covariate adjustment.

3. Robustness Assessment

The unadjusted model identified 4,444 DEGs
The adjusted model (gender + age) identified 4,284 DEGs
Approximately 82% of DEGs overlapped between the two models
This high overlap indicates that the transcriptional response associated with SARS-CoV-2 infection is largely independent of gender and age.

4. Gene Set Enrichment Analysis (GSEA)

Pathway-level interpretation was performed using clusterProfiler.

Approach:

Genes were ranked by log2 fold change
Gene symbols were mapped to ENTREZ IDs
GSEA was conducted using GO Biological Process (BP) terms
Directional GSEA was performed separately for:
Genes upregulated in POS samples
Genes downregulated in POS samples

Key Results
Differential Expression

Top upregulated genes in SARS-CoV-2 positive samples included interferon-stimulated and immune response genes such as:

CXCL10

IFIT1

IFIT2

These results are consistent with known antiviral host responses.

Pathway Enrichment

POS-upregulated pathways:

Antiviral defense response
Innate immune response
Cytokine-mediated signaling

POS-downregulated pathways:

Translation and ribosome biogenesis
Mitochondrial respiration
Oxidative phosphorylation and ATP synthesis
This pattern reflects immune activation accompanied by suppression of core metabolic and translational processes during viral infection.

Visual Outputs

The pipeline generates and saves:

PCA plots from Python preprocessing (results/python_qc/)
Volcano plots from DESeq2
Heatmaps of top DEGs
GSEA dot plots for POS-upregulated and POS-downregulated pathways
All figures are saved automatically to the results/ directory.

Limitations

Ct values were excluded from adjusted models due to substantial missing data
Bulk RNA-seq does not resolve cell-type-specific expression changes
Observational study design limits causal inference

Reproducibility

Python used for preprocessing and exploratory analysis
R (Bioconductor) used for statistical modeling and enrichment analysis
Output directories are created automatically by scripts
The pipeline is designed to be run from the project root directory

References

Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology. 2014;15:550.
https://doi.org/10.1186/s13059-014-0550-8

Wu T, Hu E, Xu S, et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. Nature Protocols. 2021;16:2769–2785.
https://doi.org/10.1038/s41596-021-00516-2

Lieberman NAP, Peddu V, Xie H, et al. In vivo antiviral host transcriptional response to SARS-CoV-2 by viral load, sex, and age. Nature Communications. 2020;11:6319.
https://doi.org/10.1038/s41467-020-20175-9
