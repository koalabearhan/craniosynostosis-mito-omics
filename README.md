# craniosynostosis-mito-omics
R scripts for craniosynostosis transcriptomics, mitochondrial gene analysis, machine-learning feature selection, molecular subtyping, GSVA, ROC/LASSO modeling, and cranial suture proteomics.

## Study overview

The workflow integrates public human craniosynostosis microarray datasets with mitochondrial gene annotation and proteomics-related analyses from a cranial suture model. The main analyses include:

- integration and batch correction of public GEO transcriptomic datasets;
- differential expression analysis between craniosynostosis and control samples;
- identification and annotation of mitochondria-related differentially expressed genes;
- machine-learning-based selection of core mitochondrial genes;
- ROC and LASSO-based model assessment;
- core-gene correlation and co-expression analyses;
- consensus clustering-based molecular subtyping;
- GSVA-based pathway activity comparison between subtypes;
- subtype-associated differential expression and enrichment analysis;
- proteomics-related candidate and pathway analyses.

## Data sources

The transcriptomic analyses use public datasets from the Gene Expression Omnibus:

- GSE27976
- GSE50796

Mitochondrial gene annotation is based on MitoCarta 3.0.

The associated proteomics dataset has been submitted to the iProX repository (ProteomeXchange Consortium) under project ID IPX0015035000, ProteomeXchange ID PXD072912.

## Software environment

The analyses were performed in R. 

## Script descriptions

| Script | Purpose |
|---|---|
| `01_data_integration.R` | Integrates transcriptomic datasets, maps probes to gene symbols, merges shared genes, and performs batch correction or normalization. |
| `02_deg_analysis.R` | Performs differential expression analysis between craniosynostosis and control samples. |
| `03_deg_enrichment.R` | Performs GO and KEGG enrichment analysis for differentially expressed genes. |
| `04_mitodeg_analysis.R` | Identifies mitochondria-related differentially expressed genes by intersecting DEGs with MitoCarta 3.0. |
| `05_machine_learning.R` | Applies machine-learning feature selection methods, including LASSO, random forest, and support vector machine-based selection. |
| `06_core_gene_correlation.R` | Performs correlation and co-expression analyses centered on the selected core genes. |
| `07_roc_analysis.R` | Generates ROC curves and calculates AUC values for individual genes and gene-panel models. |
| `08_lasso_model.R` | Builds or refits LASSO-based models using selected mitochondrial genes. |
| `09_molecular_subtyping.R` | Performs consensus clustering and defines molecular subtypes based on core-gene expression. |
| `10_subtype_gene_de_expression.R` | Performs differential gene expression analysis between molecular subtypes. |
| `11_subtype_clinical_heatmap.R` | Generates heatmaps showing subtype-related gene expression and clinical or sample annotations. |
| `12_gsva_analysis.R` | Performs GSVA pathway activity scoring and compares pathway scores between subtypes. |
| `13_pca_analysis.R` | Performs principal component analysis for data quality assessment, batch correction evaluation, or subtype visualization. |
| `14_subtype_difference_analysis.R` | Summarizes subtype-associated differences across genes, pathways, or sample-level features. |
| `15_subtype_enrichment_analysis.R` | Performs functional enrichment analysis for subtype-associated genes. |
| `16_proteomics_analysis.R` | Analyzes processed proteomics output, including candidate proteins, mitochondrial annotations, and pathway-level summaries. |
