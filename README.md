# PDAC_Atlas
Analysis workflow for all data present in our PDAC Atlas Manuscript titled, "Transfer learning associates CAFs with EMT and inflammation in tumor cells in human tumors and organoid co-culture in pancreatic ductal adenocarcinoma"

## Atlas Workflow
Scripts for the creation of the atlas, and its downstream analyses can be found under the /Atlas Directory.

## Mutliseq Workflow
Scripts for the analysis of the 12HR PDO-CAF Coculture scRNA-seq data can be found under the /Multiseq Directory.

### Data Required to follow along
* Preprocess_fromCore:
  1. RAW FASTQ
* 12H_MultiSeq_01:
  1. dgeData.rda
  2. epiMat-e74511f6-8a1d-4929-ab42-1cb67c5a4fb6-result-8pattern.rds
* 12H_MultiSeq_02:
  1. 12H_MultiSeq_ProjR_InflammationData.rda
* 12H_MultiSeq_03:
  1. 12H_MultiSeq_ProjR_InflammationData.rda
* Domino_01: 
  1. cds_combined_8_SHARED_GENES_QC_filtered_harmonized_preprocessed_aligned_manuscript_UMAP_learn_graph(6).rds
  2. auc_PDAC_PENG_EpCAF.csv
  3. all cpdb files (human signalling database directory)
  4. regulons_PDAC_PENG_EpCAF.csv
  
## Bulk RNAseq Workflow
Scripts for the analysis of the 3 untreated PDO-CAF bulk RNA-seq data can be found under the /Bulkseq Directory.
