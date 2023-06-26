# Author: Jacob Mitchell
# Date: 3/28/2022

library(ggplot2)
library(ggpubr)
library(monocle3)
library(CoGAPS)
library(dplyr)
library(forcats)

# Increase memory limit for R on Windows
if (memory.limit() < 15000){
  memory.limit(size = 15000)
}

# create sub-directory of figures for correlation plots
if (!dir.exists("results/figures/CAF_correlations")){
  dir.create("results/figures/CAF_correlations")
}

# save sessionInfo as a text document with the figures
writeLines(capture.output(sessionInfo()),
           "results/figures/CAF_correlations/sessionInfo.txt")

# CoGAPS Result
cogaps <- readRDS("CoGAPS_output/epiMat-e74511f6-8a1d-4929-ab42-1cb67c5a4fb6-result-8pattern.rds")

# full PDAC Atlas cds
cds <- readRDS("data/cds_combined_8_SHARED_GENES_QC_filtered_harmonized_preprocessed_aligned_manuscript_UMAP_learn_graph(6).rds")
# Subset to Peng and Steele manuscripts
cds_PS <- cds[,(pData(cds)$manuscript == "Peng" | 
                  pData(cds)$manuscript == "Steele")]

# cds containing epithelial cells from Peng and Steele Manuscripts
cds_epi <- readRDS("CoGAPS/cds_combined_epithelial_Peng&Steele_8_SHARED_GENES_QC_filtered_harmonized_preprocessed_aligned_manuscript_UMAP_learn_graph(6).rds")

# Incorporate pattern weights from the CoGAPS result into the cell meta data
# of the Epithelial cells from Peng and Steele

# Extract the pattern matrix
patMat <- cogaps@sampleFactors
# Merge with the column data in the cds
cds_epi@colData <- cbind(colData(cds_epi), patMat[colnames(cds_epi), ])

# create a data frame of mean pattern weights for each subject
epi_cell_meta <- as.data.frame(pData(cds_epi))

subject_pattern_weights <- epi_cell_meta %>% 
  group_by(sample_ID, TN, manuscript) %>% 
  summarize(Pattern_1_mean = mean(Pattern_1),
            Pattern_2_mean = mean(Pattern_2),
            Pattern_3_mean = mean(Pattern_3),
            Pattern_4_mean = mean(Pattern_4),
            Pattern_5_mean = mean(Pattern_5),
            Pattern_6_mean = mean(Pattern_6),
            Pattern_7_mean = mean(Pattern_7),
            Pattern_8_mean = mean(Pattern_8))

# Calculate proportions of CAFs to all cells for the subjects in Peng and Steele
PS_cell_meta <- as.data.frame(pData(cds_PS))

subject_fibroblast_proportion <- PS_cell_meta %>% 
  group_by(sample_ID, TN, manuscript) %>% 
  summarise(
    fibroblast_count = sum(TN_assigned_cell_type_immune_broad %in% c("Fibroblast")),
    
    all_count = length(TN_assigned_cell_type_immune_broad),
    
    # proportion of fibroblasts / all cells from the subject
    fibroblast_proportion = 
      sum(TN_assigned_cell_type_immune_broad %in% c("Fibroblast")) / 
      length(TN_assigned_cell_type_immune_broad),
    
    # iCAF and myCAF counts only consider CAFs also annotated as Fibroblasts
    # within the "TN_assigned_cell_type_immune_broad" classifier
    
    iCAF_count = sum(TN_assigned_cell_type_immune_broad %in% c("Fibroblast") &
                       Classifier_T_Fibroblast_only %in% c("iCAF")),
    myCAF_count = sum(TN_assigned_cell_type_immune_broad %in% c("Fibroblast") &
                        Classifier_T_Fibroblast_only %in% c("myCAF")),
    
    # proportion of iCAFs / fibroblasts from the subject
    iCAF_proportion = 
      if (sum(TN_assigned_cell_type_immune_broad %in% c("Fibroblast") != 0)) {
        sum(sum(TN_assigned_cell_type_immune_broad %in% c("Fibroblast") &
                  Classifier_T_Fibroblast_only %in% c("iCAF"))) / 
          sum(TN_assigned_cell_type_immune_broad %in% c("Fibroblast"))
      } else {NA},
    
    # proportion of myCAFs / fibroblasts from the subject
    myCAF_proportion = 
      if (sum(TN_assigned_cell_type_immune_broad %in% c("Fibroblast") != 0)) {
        sum(sum(TN_assigned_cell_type_immune_broad %in% c("Fibroblast") &
                  Classifier_T_Fibroblast_only %in% c("myCAF"))) / 
          sum(TN_assigned_cell_type_immune_broad %in% c("Fibroblast"))
      } else {NA}
  )

# merge data frames of mean pattern weights with CAF proportions
subject_pat_caf <- right_join( x = subject_pattern_weights,
                               y = subject_fibroblast_proportion,
                               by = c("sample_ID", "TN", "manuscript"))

# save mean pattern weights and CAF proportions for all subjects
write.csv(subject_pat_caf,
          file = "results/figures/CAF_correlations/subject_pattern_weights_CAF_proportions.csv")

## Plotting ###################################################################

# subsetting the data frame of mean pattern weights and cell proportions
# for plotting

# Tumor samples, Peng & Steele
tumor_pat_caf <- subject_pat_caf[subject_pat_caf$TN == "T",]
# Tumor samples, Peng
peng_tumor <- subject_pat_caf[subject_pat_caf$TN == "T" &
                                subject_pat_caf$manuscript == "Peng",]
# Tumor samples, Steele
steele_tumor <- subject_pat_caf[subject_pat_caf$TN == "T" &
                                  subject_pat_caf$manuscript == "Steele",]

# list of cellular proportions to plot
caf_types <- c("fibroblast_proportion", "iCAF_proportion", "myCAF_proportion")

# plot regressions for Peng and Steele Tumor Samples with Pattern 7 weight
for (caf in caf_types){
  plot <- ggplot(tumor_pat_caf, aes(x = tumor_pat_caf[[caf]],
                                    y = Pattern_7_mean,
                                    shape = TN)) +
    geom_point(size = 8, aes(col = manuscript)) +
    geom_smooth(size = 1, method = "lm", color = "red") +
    stat_cor(method = "pearson",
             label.x = 0,
             label.y = 0.6) +
    xlim(0, 1) +
    ylim(0, max(tumor_pat_caf$Pattern_7_mean)) +
    xlab(gsub("_", " ", caf)) +
    ylab("Mean Pattern 7 Weight in Epithelium") +
    scale_color_manual(values = c("#009ddc", "#963d97")) +
    theme_minimal() +
    guides(shape = "none")
  ggsave(paste0("results/figures/CAF_correlations/",
                "Peng+Steele_Tumor_",
                caf,
                "_Pattern7_Regression.pdf"),
         plot = plot,
         width = unit(8, "in"),
         height = unit(6, "in"))
}

# plot regressions for Peng Tumor Samples with Pattern 7 weight
for (caf in caf_types){
  plot <- ggplot(peng_tumor, aes(x = peng_tumor[[caf]],
                                 y = Pattern_7_mean,
                                 shape = TN)) +
    geom_point(size = 8, aes(col = manuscript)) +
    geom_smooth(size = 1, method = "lm", color = "red") +
    stat_cor(method = "pearson",
             label.x = 0,
             label.y = 0.6) +
    xlim(0, 1) +
    ylim(0, max(tumor_pat_caf$Pattern_7_mean)) +
    xlab(gsub("_", " ", caf)) +
    ylab("Mean Pattern 7 Weight in Epithelium") +
    scale_color_manual(values = c("#009ddc")) +
    theme_minimal() +
    guides(shape = "none")
  ggsave(paste0("results/figures/CAF_correlations/",
                "Peng_Tumor_",
                caf,
                "_Pattern7_Regression.pdf"),
         plot = plot,
         width = unit(8, "in"),
         height = unit(6, "in"))
}

# plot regressions for Steele Tumor Samples with Pattern 7 weight
for (caf in caf_types){
  plot <- ggplot(steele_tumor, aes(x = steele_tumor[[caf]],
                                   y = Pattern_7_mean,
                                   shape = TN)) +
    geom_point(size = 8, aes(col = manuscript)) +
    geom_smooth(size = 1, method = "lm", color = "red") +
    stat_cor(method = "pearson",
             label.x = 0,
             label.y = 0.6) +
    xlim(0, 1) +
    ylim(0, max(tumor_pat_caf$Pattern_7_mean)) +
    xlab(gsub("_", " ", caf)) +
    ylab("Mean Pattern 7 Weight in Epithelium") +
    scale_color_manual(values = c("#963d97")) +
    theme_minimal() +
    guides(shape = "none")
  ggsave(paste0("results/figures/CAF_correlations/",
                "Steele_Tumor_",
                caf,
                "_Pattern7_Regression.pdf"),
         plot = plot,
         width = unit(8, "in"),
         height = unit(6, "in"))
}

# loop through all patterns
patterns <- paste0("Pattern_", 1:8)
for(pat in patterns){
  for (caf in caf_types){
    plot <- ggplot(tumor_pat_caf, aes(x = tumor_pat_caf[[caf]],
                                      y = .data[[paste0(pat, "_mean")]],
                                      shape = TN)) +
      geom_point(size = 8, aes(col = manuscript)) +
      geom_smooth(size = 1, method = "lm", color = "red") +
      stat_cor(method = "pearson",
               label.x = 0,
               label.y = 0.6) +
      xlim(0, 1) +
      ylim(0, max(tumor_pat_caf[[paste0(pat, "_mean")]])) +
      xlab(gsub("_", " ", caf)) +
      ylab(paste0("Mean ", pat, " Weight in Epithelium")) +
      scale_color_manual(values = c("#009ddc", "#963d97")) +
      theme_minimal() +
      guides(shape = "none")
    ggsave(paste0("results/figures/CAF_correlations/",
                  "Peng+Steele_Tumor_",
                  caf,
                  "_", pat, "_Regression.pdf"),
           plot = plot,
           width = unit(8, "in"),
           height = unit(6, "in"))
  }
}

