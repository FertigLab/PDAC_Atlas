# Author: Jacob Mitchell

# Preprocessing of cds object for creation of matrix, geneNames, and sampleNames
# to be passed into CoGAPS

# Import Libraries
library(monocle3)
library(CoGAPS)
library(Matrix)
library(dplyr)

# increase memory limit to avoid crashing during manipulation of expression matrix
if (memory.limit() < 15000){
  memory.limit(size = 15000)
}


# Read in complete PDAC dataset
cds <- readRDS("data/cds_combined_8_SHARED_GENES_QC_filtered_harmonized_preprocessed_aligned_manuscript_UMAP_learn_graph(6).rds")

# Subset Data
# Manuscripts: Peng, Steele
# Assigned Cell Types: Epithelial_cancer, Epithelial_normal, Epithelial_unspecified

cds_cogaps <- cds[,(pData(cds)$manuscript == "Peng" | 
                      pData(cds)$manuscript == "Steele") & 
                    (!is.na(pData(cds)$TN_assigned_cell_type_immune_broad)) & 
                    (pData(cds)$TN_assigned_cell_type_immune_broad == "Epithelial_cancer" |
                       pData(cds)$TN_assigned_cell_type_immune_broad == "Epithelial_normal" |
                       pData(cds)$TN_assigned_cell_type_immune_broad == "Epithelial_unspecified")]

# Save the subset cds object
saveRDS(cds_cogaps, file = "CoGAPS/cds_combined_epithelial_Peng&Steele_8_SHARED_GENES_QC_filtered_harmonized_preprocessed_aligned_manuscript_UMAP_learn_graph(6).rds")

# Save the expression matrix
mat <- counts(cds_cogaps)

# Use feature data from the cds object to convert to gene short names
cds_gene_meta <- data.frame(fData(cds_cogaps))
mat@Dimnames[[1]] <- cds_gene_meta$gene_short_name

# Preprocessing on expression matrix for CoGAPS

# Remove cells with no signal (all expression values in column are 0)
mat <- mat[, apply(mat, 2, max) > 0]
# Remove genes with no signal/constant signal (stdev of row is 0)
mat <- mat[apply(mat, 1, sd) > 0, ]
# Remove genes with no signal (max of row is 0)
mat <- mat[apply(mat, 1, max) > 0, ]
# Log transform based on the counts data being whole numbers 
mat <- log2(mat+1)

# Check min of transformed data
min(mat)
# Recheck the dim
dim(mat)
range(mat)
mean(mat)

# Convert to Sparse Matrix
sparse.mat <- Matrix(mat, sparse = TRUE)

# Save the matrix
writeMM(obj = sparse.mat, file = "CoGAPS/epiMat.mtx")

# Save row and column names for use in CoGAPS parameter creation
geneNames<- rownames(mat)
sampleNames<- colnames(mat)

saveRDS(geneNames, file = "CoGAPS/geneNames.rds")
saveRDS(sampleNames, file = "CoGAPS/sampleNames.rds")


