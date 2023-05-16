# Author: Jacob Mitchell

# Collect pattern markers from the CoGAPS result
# Carry out overrepresentation analysis with MsigDB Hallmarks for each pattern's markers
# Add pattern weights to the metadata of the cds object
# Plot pattern weights on the UMAP embedding of the epithelial cells
# Plot pattern weights in Epithelial cell types in control vs tumor subjects

# import necessary libraries
library(ggplot2, quietly = TRUE)
library(ggpubr, quietly = TRUE)
library(monocle3, quietly = TRUE)
library(CoGAPS, quietly = TRUE)
library("fgsea", quietly = TRUE)
library("msigdbr", quietly = TRUE)
library(biomaRt, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(forcats, quietly = TRUE)

## Ensure that the directories are ready for the output files #################
if (!dir.exists("results/figures")){
  dir.create("results/figures")
}
if (!dir.exists("results/figures/Pattern_MsigDB_ORA")){
  dir.create("results/figures/Pattern_MsigDB_ORA")
}
if (!dir.exists("results/MsigDB_ORA")){
  dir.create("results/MsigDB_ORA")
}
if (!dir.exists("results/figures/UMAP")){
  dir.create("results/figures/UMAP")
}
if (!dir.exists("results/figures/Pattern_Weights")){
  dir.create("results/figures/Pattern_Weights")
}

# read in the CoGAPS results object and cds object subset to epithelial cells
cogaps <- readRDS("CoGAPS_output/epiMat-e74511f6-8a1d-4929-ab42-1cb67c5a4fb6-result-8pattern.rds")
cds <- readRDS("CoGAPS/cds_combined_epithelial_Peng&Steele_8_SHARED_GENES_QC_filtered_harmonized_preprocessed_aligned_manuscript_UMAP_learn_graph(6).rds")

# List of MsigDB hallmarks and the genes in each set
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- hallmark_df %>% split(x = .$gene_symbol, f = .$gs_name)

# obtain universe of all human hgnc gene symbols from ensembl
mart <- useEnsembl('genes', dataset = "hsapiens_gene_ensembl") #GRCh38
Hs_genes <- getBM("hgnc_symbol", mart = mart)


## Find Pattern Markers #######################################################
# Extract the pattern matrix
patMat <- cogaps@sampleFactors
# Create temporary CDS object
tempCds <- cds
# Merge with the column data in the cds
tempCds@colData <- cbind(colData(tempCds), patMat[colnames(tempCds), ])
# Extract lists of pattern marker genes
# use "threshold" parameters to limit to the highest ranked markers
patternMarkerResults <- patternMarkers(cogaps, threshold = "cut", lp = NA, axis = 1)
names(patternMarkerResults$PatternMarkers) <- colnames(patternMarkerResults$PatternMarkerRanks)

# Save the pattern markers object as an RDS for downstream use
saveRDS(patternMarkerResults,
        file = "results/Pattern_Marker_Results.rds")

# save monocle object with incorporated pattern weights
saveRDS(tempCds,
        file = "results/cds_epithelial_CoGAPS_weights.rds")

# List of each gene as a pattern marker
PMlist <- patternMarkerResults$PatternMarkers

patternmarks <- data.frame(lapply(patternMarkerResults$PatternMarkers, "length<-",
                                  max(lengths(patternMarkerResults$PatternMarkers))))

# Save the dataframe of pattern markers
write.csv(patternmarks,
          file = "results/Pattern_Marker_Genes.csv",
          na = "")

# Create a tumor proximity identifier labeling cells by whether they came from
# a tumor or non-tumor subject (TN) and the epithelial cell type
# (TN_assigned_cell_type_immune_broad)

# pattern enrichment in Adjacent normal epithelium in of patients with tumors
# create a new meta factor for TN + cell type immune broad
colData(tempCds)$TumorProx <- paste(as.character(colData(tempCds)$TN),"_", 
                                    as.character(colData(tempCds)$TN_assigned_cell_type_immune_broad), sep = "")

# metadata dataframe created after the introduction of new categories
cds_cell_meta <- data.frame(pData(tempCds))

# Defining tissue factors for plotting of epithelial type barplots
# Tissue vector for setting order of how the epithelial celltypes are plotted
Tissue <- factor(cds_cell_meta[,"TumorProx"])
levels(Tissue)[1:6] <- c("Control - Cancer",
                         "Control - Ductal Epithelium",
                         "Control - Unspecified",
                         "Tumor - Cancer",
                         "Tumor - Ductal Epithelium",
                         "Tumor - Unspecified")

# Subset of metadata including only the ductal epithelium in control and tumor
tnde <- cds_cell_meta[cds_cell_meta$TumorProx %in% c("N_Epithelial_normal","T_Epithelial_normal"),]
Tissue2 <- factor(tnde[,"TN"])
levels(Tissue2)[1:2] <- c("Control", "Tumor")

## Overall Results ############################################################

# Plot UMAP of cells colored by cell type
celltype <- plot_cells(tempCds, color_cells_by = "TN_assigned_cell_type_immune_broad",
                       label_cell_groups = FALSE, show_trajectory_graph = FALSE)
ggsave(paste0("results/figures/UMAP/Cell_Type_UMAP.pdf"),
       plot = celltype,
       width = unit(8, "in"),
       height = unit(5, "in"))

# Plot UMAP of cells colored by Tumor Proximity
prox <- plot_cells(tempCds, color_cells_by = "TumorProx",
                   label_cell_groups = FALSE, show_trajectory_graph = FALSE) +
  scale_color_manual(values=c("#FFC252",
                              "#FFA500",
                              "#BF7C00",
                              "#D10000",
                              "#8B0000",
                              "#4F0000"),
                     labels = c("Control - Cancer",
                                "Control - Ductal Epithelium",
                                "Control - Unspecified",
                                "Tumor - Cancer",
                                "Tumor - Ductal Epithelium",
                                "Tumor - Unspecified"),
                     guide_legend(title="Pancreas Tissue")) +
  guides(color=guide_legend(title="Pancreas Tissue", override.aes = list(size=4)))
ggsave(paste0("results/figures/UMAP/Tumor_Proximity_UMAP.pdf"),
       plot = prox,
       width = unit(8, "in"),
       height = unit(5, "in"))

# Per Pattern Figure Plotting Loop ############################################

for (pattern in names(PMlist)) {
  
  # Over-representation analysis of pattern markers ###########################
  suppressWarnings(
    result <- fora(pathways = hallmark_list,
                   genes = PMlist[[pattern]],
                   universe = Hs_genes$hgnc_symbol,
                   maxSize=2038)
  )
  
  # Append ratio of overlap/# of genes in hallmark (k/K)
  result$"k/K" <- result$overlap/result$size
  
  # Append log-transformed HB-adjusted q value (-10 * log(adjusted p value))
  result$"neg.log.q" <- (-10) * log10(result$padj)
  
  # Append shortened version of hallmark's name without "HALLMARK_"
  result$MsigDB_Hallmark <- substr(result$pathway, 10, nchar(result$pathway))
  
  # Reorder dataframe by ascending adjusted p value
  result <- mutate(result, MsigDB_Hallmark=fct_reorder(MsigDB_Hallmark, - padj))
  
  # rearrange columns to move overlap genes to the end and convert to vector for
  # compatibility with saving as csv
  result <- relocate(result, "overlapGenes", .after = "MsigDB_Hallmark")
  result <- relocate(result, "neg.log.q", .after = "padj")
  result <- result %>% mutate(overlapGenes = sapply(overlapGenes, toString))
  
  # save pattern marker overlaps as a csv file
  write.csv(result,
            file = paste0("results/MsigDB_ORA/", pattern, "_MsigDB_OverRepresentationAnalysis.csv"))
  
  # for plotting, limit the results to significant over-representation
  result <- result[result$padj < 0.05,]
  
  #plot and save the waterfall plot of ORA p-values
  plot <- ggplot(result, aes_string(y = "neg.log.q", x = "MsigDB_Hallmark", fill = "MsigDB_Hallmark")) +
    ## Specifies barplot
    geom_col() +
    ## Rename y axis
    ylab("-10*log10(FDR q-value)") + 
    ## Flips the coordinates
    coord_flip() +
    ## Makes the background white
    theme_minimal() +
    ## Add title
    ggtitle(paste0(pattern, ": Overrepresented MsigDB Hallmarks")) +
    ## This creates the dotted line at .05 value 
    geom_hline(yintercept=c(13.0103), linetype="dotted") + # Add veritcle line to show significances
    ## Adds the q values
    geom_text(aes(label=format(signif(padj, 4))), hjust = -.04) +
    ## Removes legend
    theme(legend.position = "none") +
    ## specifies limits 
    ylim(0, ceiling(max(result$"neg.log.q")) + (max(result$"neg.log.q")/4))
  ggsave(paste0("results/figures/Pattern_MsigDB_ORA/", pattern, "_ORA.pdf"),
         plot = plot,
         width = unit(10, "in"),
         height = unit(6, "in"),
         device = "pdf")
  
  # UMAP Embedding of CoGAPS Pattern Weight ###################################
  umap <- plot_cells(tempCds, color_cells_by = pattern,
                     label_cell_groups = FALSE, show_trajectory_graph = FALSE)
  ggsave(paste0("results/figures/UMAP/", pattern, "_UMAP.pdf"),
         plot = umap,
         width = unit(8, "in"),
         height = unit(5, "in"))
  
  # Jittered violin plots of pattern weights in tumor proximity groups ########

  plot <- ggplot(cds_cell_meta, 
                 aes(factor(cds_cell_meta[,"TumorProx"]), cds_cell_meta[,pattern]))
  plot + 
    geom_boxplot(aes(fill = Tissue), outlier.shape = NA) +
    geom_jitter(shape = 21, aes(fill = Tissue),
                width = 0.1,
                show.legend = FALSE) +
    labs(fill = "Pancreas Tissue") +
    scale_fill_manual(values = c("Control - Cancer" = "#FFC252",
                                 "Control - Ductal Epithelium" = "#FFA500",
                                 "Control - Unspecified" = "#BF7C00",
                                 "Tumor - Cancer" = "#D10000",
                                 "Tumor - Ductal Epithelium" = "#8B0000",
                                 "Tumor - Unspecified" = "#4F0000")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    geom_hline(yintercept = 0, alpha = 0.1) +
    xlab("Epithelial Cell Type") + ylab(paste0(pattern, " Weight")) +
    scale_x_discrete(labels=c("Cancer", "Ductal Epithelium", "Unspecified",
                              "Cancer", "Ductal Epithelium", "Unspecified"))
  ggsave(paste0("results/figures/Pattern_Weights/", pattern, "_TN_Epithelial_Weights.pdf"),
         width = unit(8, "in"),
         height = unit(5, "in"))
  
  # Jittered violin plots of pattern weights control/tumor ductal epithelium ##
  plot2 <- ggplot(tnde, 
                  aes(factor(tnde[,"TumorProx"]), tnde[,pattern]))
  plot2 + 
    geom_boxplot(aes(fill = Tissue2), outlier.shape = NA) +
    geom_jitter(shape = 21, aes(fill = Tissue2),
                width = 0.1,
                show.legend = FALSE) +
    stat_compare_means(comparisons = list(c("T_Epithelial_normal","N_Epithelial_normal")),
                       method = "wilcox") +
    labs(fill = "Pancreas Tissue") +
    scale_fill_manual(values = c("Control" = "#FFA500",
                                 "Tumor" = "#8B0000")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    geom_hline(yintercept = 0, alpha = 0.1) +
    xlab("Ductal Epithelium") + ylab(paste0(pattern, " Weight")) +
    ylim(0,1.1) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_x_discrete(labels=NULL)
  ggsave(paste0("results/figures/Pattern_Weights/", pattern, "_TN_Ductal_Weights.pdf"),
         width = unit(4, "in"),
         height = unit(5, "in"))
}

