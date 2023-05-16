###########################
### clear the workspace
rm(list=ls())


###########################
### Reads the coammand line arguments
args <- commandArgs(TRUE)


###########################
### set values to null
dge_path <- ""


###########################
### go through the command line arguments and look for matches to the
### defined flags. When those are found, set the appropriate value.
i <- 1
while (i < length(args)) {
  opt <- args[i]
  val <- args[i+1]
  if (opt == "-i") { dge_path <- val }
  i <- i+2
}


###########################
### make sure that the value is set, and if not, complain and quit
if (dge_path == "")
{
  message("Usage: Rscript UMAP.R -i dge_path")
  q()
}


###########################
### load library
library(Seurat)
library(dplyr)
library(cowplot)


###########################
### load data
load("deMULTIplex_objs.rda")
rm(list=c("aggrData", "bar.table", "bar.table.full", "bar.tsne", "final.calls"))


###########################
### read count data
umis <- Read10X(data.dir=dge_path, strip.suffix=TRUE)

### create a seurat object
samplename <- sub("_GE", "", sub(".*/", "", sub(".filtered_feature_bc_matrix", "", dge_path)))
dgeData <- CreateSeuratObject(counts = as.matrix(umis), project = samplename)

### assign MULTIseq classificatio
## Check order
if(all(umis.barcodes.merged$FeatureIndex==umis@Dimnames[[2]])) 
{
	print("Order looks good")
}
dgeData$MULTIseqClassification <- umis.barcodes.merged$Samplename


###########################
### QC and selecting cells for further analysis
dgeData[["percent.mt"]] <- PercentageFeatureSet(dgeData, pattern = "^MT-")

# Visualize QC metrics as a violin plot
pdf("violin_plot.pdf", height=6, width=18)
VlnPlot(dgeData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf("scatter_plot.pdf", height=6, width=12)
plot1 <- FeatureScatter(dgeData, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(dgeData, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

# Normalizing the data
dgeData <- NormalizeData(dgeData, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features
dgeData <- FindVariableFeatures(dgeData, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dgeData), 10)

# plot variable features with and without labels
pdf("variable_feature_plot.pdf", height=6, width=12)
plot1 <- VariableFeaturePlot(dgeData)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()

# Scaling the data
all.genes <- rownames(dgeData)
dgeData <- ScaleData(dgeData, features = all.genes)

# Perform linear dimensional reduction
dgeData <- RunPCA(dgeData, features = VariableFeatures(object = dgeData))

pdf("dimensional_reduction_plot.pdf", height=6, width=12)
VizDimLoadings(dgeData, dims = 1:2, reduction = "pca")
DimPlot(dgeData, reduction = "pca")
dev.off()

pdf("heatmap_plot.pdf", height=10, width=10)
DimHeatmap(dgeData, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(dgeData, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()


###########################
### Determine the ‘dimensionality’ of the dataset
pdf("elbow_plot.pdf", height=10, width=10)
ElbowPlot(dgeData, ndims=50)
dev.off()


###########################
### Cluster the cells
dgeData <- FindNeighbors(dgeData, dims = 1:30)
dgeData <- FindClusters(dgeData, resolution = 0.5)


###########################
### Run non-linear dimensional reduction (UMAP/tSNE)
dgeData <- RunUMAP(dgeData, dims = 1:30)


###########################
### make umap plot
dgeData$dummy.id <- "Background"

mySel_CAF1 <- list()
mySel_CAF1[[1]] <- grep("CAF1", dgeData$MULTIseqClassification)
names(mySel_CAF1) <- paste0(paste0("CAF1:", length(mySel_CAF1[[1]])), " cells")

mySel_CC1 <- list()
mySel_CC1[[1]] <- grep("CC1", dgeData$MULTIseqClassification)
names(mySel_CC1) <- paste0(paste0("CC1:", length(mySel_CC1[[1]])), " cells")

mySel_Org1 <- list()
mySel_Org1[[1]] <- grep("Org1", dgeData$MULTIseqClassification)
names(mySel_Org1) <- paste0(paste0("Org1:", length(mySel_Org1[[1]])), " cells")

mySel_Panc1 <- list()
mySel_Panc1[[1]] <- grep("Panc1", dgeData$MULTIseqClassification)
names(mySel_Panc1) <- paste0(paste0("Panc1:", length(mySel_Panc1[[1]])), " cells")

mySel_Negative <- list()
mySel_Negative[[1]] <- grep("Negative", dgeData$MULTIseqClassification)
names(mySel_Negative) <- paste0(paste0("Negative:", length(mySel_Negative[[1]])), " cells")

mySel <- list()
mySel[[1]] <- grep("CAF1", dgeData$MULTIseqClassification)
mySel[[2]] <- grep("CC1", dgeData$MULTIseqClassification)
mySel[[3]] <- grep("Org1", dgeData$MULTIseqClassification)
mySel[[4]] <- grep("Panc1", dgeData$MULTIseqClassification)
names(mySel) <- c(paste0(paste0("CAF1:", length(mySel[[1]])), " cells"), paste0(paste0("CC1:", length(mySel[[2]])), " cells") , paste0(paste0("Org1:", length(mySel[[3]])), " cells"),
			paste0(paste0("Panc1:", length(mySel[[4]])), " cells"))

pdf("umap_plot.pdf", height=7.5, width=15)
p1 <- DimPlot(dgeData, reduction = "umap")
p2 <- DimPlot(dgeData, reduction = "umap", split.by = "MULTIseqClassification")
p3 <- DimPlot(dgeData, reduction = "umap", group.by = "dummy.id", cols="grey", cells.highlight=mySel, cols.highlight=c("blue", "red", "gold", "darkcyan"))
p4 <- DimPlot(dgeData, reduction = "umap", group.by = "dummy.id", cols="grey", cells.highlight=mySel_CAF1, sizes.highlight=2)
p5 <- DimPlot(dgeData, reduction = "umap", group.by = "dummy.id", cols="grey", cells.highlight=mySel_CC1, sizes.highlight=2)
p6 <- DimPlot(dgeData, reduction = "umap", group.by = "dummy.id", cols="grey", cells.highlight=mySel_Org1, sizes.highlight=2)
p7 <- DimPlot(dgeData, reduction = "umap", group.by = "dummy.id", cols="grey", cells.highlight=mySel_Panc1, sizes.highlight=2)
p8 <- DimPlot(dgeData, reduction = "umap", group.by = "dummy.id", cols="grey", cells.highlight=mySel_Negative, sizes.highlight=2)
plot_grid(p1, p2)
plot_grid(p1, p3)
plot_grid(p1, p4)
plot_grid(p1, p5)
plot_grid(p1, p6)
plot_grid(p1, p7)
plot_grid(p1, p8)
dev.off()


###########################
### make violin & UMAP plots for FAP & CDH1 genes
pdf("FAP_CDH1.pdf", height=7.5, width=15)
VlnPlot(dgeData, features = c("FAP", "CDH1"))
FeaturePlot(dgeData, features = c("FAP", "CDH1"))
dev.off()

###########################
### Finding differentially expressed features (cluster biomarkers)
# find markers for every cluster compared to all remaining cells, report only the positive ones
dgeData.markers <- FindAllMarkers(dgeData, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- dgeData.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

pdf("features_heatmap_plot.pdf", height=25, width=25)
DoHeatmap(dgeData, features = top10$gene) + NoLegend()
dev.off()

# find all markers of a single cluster compared to all other cells.
numOfClusters <- 12
for (i in 1:numOfClusters)
{
	cluster.markers <- FindMarkers(dgeData, ident.1 = i-1, min.pct = 0.25)
	cluster.markers$gene <- rownames(cluster.markers)
	cluster.markers <- cluster.markers[,c(6,1:5)]
	write.table(cluster.markers, paste0(paste0("cluster", i-1), "_markers.txt"), row.names=F, sep="\t", quote=F)
}


###########################
### save data
save(dgeData, file="dgeData.rda")


###########################
### sessionInfo and clean and quit
sessionInfo()
date()
rm(list=ls())


###########################
####quit
q("yes")
