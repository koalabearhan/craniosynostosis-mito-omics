# craniosynostosis-mito-omics
Scripts for craniosynostosis transcriptomics, mitochondrial gene analysis, machine-learning feature selection, molecular subtyping, GSVA, ROC/LASSO modeling, and cranial suture proteomics.

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
| `14_subtype_difference_analysis.R` | Summarizes subtype-associated differences across genes, pathways, or sample-level features. |
| `15_subtype_enrichment_analysis.R` | Performs functional enrichment analysis for subtype-associated genes. |
| `16_proteomics_analysis.R` | Analyzes processed proteomics output, including candidate proteins, mitochondrial annotations, and pathway-level summaries. |
