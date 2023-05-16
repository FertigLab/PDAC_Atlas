# Author: Jacob Mitchell

# Creation of the CoGAPS parameters file for running CoGAPS on AWS

library(CoGAPS)
library(Matrix)
geneNames <- readRDS("CoGAPS/geneNames.rds")
sampleNames <- readRDS("CoGAPS/sampleNames.rds")
EMT_Params_MEL <- readRDS("data/emtParams (2).rds")
EMT_Params_MEL@geneNames <- geneNames
EMT_Params_MEL@sampleNames <- sampleNames
EMT_Params_MEL@distributed <- "Single-Cell"
EMT_Params_MEL@nSets <- 15
EMT_Params_MEL@nPatterns <- 8
EMT_Params_MEL@nIterations <- 50000
EMT_Params_MEL@cut <- 10
EMT_Params_MEL@minNS <- 8
EMT_Params_MEL@maxNS <- 23
EMT_Params_MEL
saveRDS(EMT_Params_MEL, file = "CoGAPS/08_nPattern_epiParams.rds")
