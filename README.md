# PDAC_Atlas
Analysis workflow for all data present in our PDAC Atlas Manuscript titled, "Transfer learning associates CAFs with EMT and inflammation in tumor cells in human tumors and organoid co-culture in pancreatic ductal adenocarcinoma"

## Atlas Workflow
Scripts for the creation of the atlas, and its downstream analyses can be found under the /Atlas Directory.

### Data Required to follow along
* PDAC atlas engineering_051622_BK.Rmd:
  1. RAW Counts from all manuscripts
* PDAC atlas analysis_051622_BK.Rmd:
  1. cds_combined_8_SHARED_GENES_QC_filtered_harmonized_preprocessed_aligned_manuscript_UMAP_learn_graph.rds
* PDAC_MHCII_DifferentialExprs.Rmd:
  1. cds_duct_SHARED_GENES_QC_filtered_harmonized_preprocessed_aligned_manuscript_UMAP_learn_graph.rds
* PDAC_CAF_MHCII_DifferentialExprs.Rmd:
  1. cds_combined_8_SHARED_GENES_QC_filtered_harmonized_preprocessed_aligned_manuscript_UMAP_learn_graph.rds
* 06222022_STEELE_PYSCENIC_HG38only.sh & Domino_PDAC_EpCAF_STEELE_HG38only.Rmd: 
  1. cds_combined_8_SHARED_GENES_QC_filtered_harmonized_preprocessed_aligned_manuscript_UMAP_learn_graph.rds
* 06222022_PENG_PYSCENIC_HG38only.sh & Domino_PDAC_EpCAF_PENG_HG38only.Rmd: 
  1. cds_combined_8_SHARED_GENES_QC_filtered_harmonized_preprocessed_aligned_manuscript_UMAP_learn_graph.rds
  
## Mutliseq Workflow
Scripts for the analysis of the 12HR PDO-CAF Coculture scRNA-seq data can be found under the /Multiseq Directory.

### Data Required to follow along
* Preprocess_fromCore:
  1. RAW FASTQ
* 12H_MultiSeq_01_ProjR_Inflammation.Rmd:
  1. dgeData.rda
  2. epiMat-e74511f6-8a1d-4929-ab42-1cb67c5a4fb6-result-8pattern.rds
* 12H_MultiSeq_02_Moffit.Rmd:
  1. 12H_MultiSeq_ProjR_InflammationData.rda

## Bulk RNAseq Workflow
Scripts for the analysis of the 3 untreated PDO-CAF bulk RNA-seq data can be found under the /Bulkseq Directory.

### Data Required to follow along
* 00_Extract_Sample_Names.Rmd:
  1. Raw Data Folder names
* 01_salmonQuant.sh & 01_SG_Novogene_12202022_SalmonSetup.R:
  1. RAW FASTQ files from Novogene
* 02_SalmonQuant_to_txi.Rmd:
  1. Output of Salmon.sh (Salmon_Quants/)
  2. tgMap.tsv
* 03_PDACAtlas_CoGAPS_Projection.Rmd:
  1. txi_data_all.rds
  2. SG_Novogene_Annot.csv
  3. epiMat-e74511f6-8a1d-4929-ab42-1cb67c5a4fb6-result-8pattern.rds
