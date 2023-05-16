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
  message("Usage: Rscript deMULTIplex.R -i dge_path")
  q()
}


###################################################
### Load libraries
library(deMULTIplex)
library(ggplot2)
library(Seurat)


###################################################
### read MULTI-seq sample barcode alignments
bar.table <- read.table("barMatrix.csv", header=T, sep=",", stringsAsFactors=F)
rownames(bar.table) <- bar.table[,1]
bar.table <- bar.table[,-c(1, ncol(bar.table)-1, ncol(bar.table))]
apply(bar.table, 2, function(x) {length(x[x>0])})


###################################################
### Step 2: Visually inspect sample barcode quality
## Visualize barcode space
bar.tsne <- barTSNE(bar.table) 

pdf("bc.check.pdf")
for (i in 3:ncol(bar.tsne)) {
    g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
    print(g)
}
dev.off()


###################################################
### Step 3: Sample Classification
## Round 1 -----------------------------------------------------------------------------------------------------
## Perform Quantile Sweep
bar.table.full <- bar.table
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

## Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results1 <- findThresh(call.list=bar.table_sweep.list)
pdf("inter_maxima_quantile_round1.pdf", width=12, height=6)
ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + labs(x="Inter-maxima qunatile", y="Proportion") + 
  geom_vline(xintercept=threshold.results1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
dev.off()

## Finalize round 1 classifications, remove negative cells
round1.calls <- classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema))
neg.cells <- names(round1.calls)[which(round1.calls == "Negative")]
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

## Round 2 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

threshold.results2 <- findThresh(call.list=bar.table_sweep.list)
pdf("inter_maxima_quantile_round2.pdf", width=12, height=6)
ggplot(data=threshold.results2$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + labs(x="Inter-maxima qunatile", y="Proportion") +
  geom_vline(xintercept=threshold.results2$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
dev.off()
round2.calls <- classifyCells(bar.table, q=findQ(threshold.results2$res, threshold.results2$extrema))
neg.cells <- c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

## Round 3 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

threshold.results3 <- findThresh(call.list=bar.table_sweep.list)
pdf("inter_maxima_quantile_round3.pdf", width=12, height=6)
ggplot(data=threshold.results3$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + labs(x="Inter-maxima qunatile", y="Proportion") +
  geom_vline(xintercept=threshold.results3$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
dev.off()
round3.calls <- classifyCells(bar.table, q=findQ(threshold.results3$res, threshold.results3$extrema))
neg.cells <- c(neg.cells, names(round3.calls)[which(round3.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

## Round 4 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

threshold.results4 <- findThresh(call.list=bar.table_sweep.list)
pdf("inter_maxima_quantile_round4.pdf", width=12, height=6)
ggplot(data=threshold.results4$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + labs(x="Inter-maxima qunatile", y="Proportion") +
  geom_vline(xintercept=threshold.results4$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
dev.off()
round4.calls <- classifyCells(bar.table, q=findQ(threshold.results4$res, threshold.results4$extrema))
neg.cells <- c(neg.cells, names(round4.calls)[which(round4.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

## Round 5 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
    print(q) 
    n <- n + 1
    bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
    names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

threshold.results5 <- findThresh(call.list=bar.table_sweep.list)
pdf("inter_maxima_quantile_round5.pdf", width=12, height=6)
ggplot(data=threshold.results5$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + labs(x="Inter-maxima qunatile", y="Proportion") +
    geom_vline(xintercept=threshold.results5$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
dev.off()
round5.calls <- classifyCells(bar.table, q=findQ(threshold.results5$res, threshold.results5$extrema))
neg.cells <- c(neg.cells, names(round5.calls)[which(round5.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

## Round 6 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
    print(q)
    n <- n + 1
    bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
    names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

threshold.results6 <- findThresh(call.list=bar.table_sweep.list)
pdf("inter_maxima_quantile_round6.pdf", width=12, height=6)
ggplot(data=threshold.results6$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + labs(x="Inter-maxima qunatile", y="Proportion") +
    geom_vline(xintercept=threshold.results6$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
dev.off()
round6.calls <- classifyCells(bar.table, q=findQ(threshold.results6$res, threshold.results6$extrema))
neg.cells <- c(neg.cells, names(round6.calls)[which(round6.calls == "Negative")])
paste0("NUMBER OF NEGATIVE CELLS LEFT: ", length(names(round6.calls)[which(round6.calls == "Negative")]))

## Repeat until all no negative cells remain (usually 3 rounds)...
final.calls <- c(round6.calls, rep("Negative",length(neg.cells)))
names(final.calls) <- c(names(round6.calls),neg.cells)


###########################
### gete demultiplexing summary
barcodes <- read.table("ssSeq_Barcode.csv", sep=",")
colnames(barcodes) <- c("Samplename", "SampleIndex")

demultiplexSummary <- data.frame(Number_of_Cells=as.vector(table(final.calls)), SampleIndex=names(table(final.calls)))
demultiplexSummary <- merge(demultiplexSummary, barcodes, all=T)
write.table(demultiplexSummary, file="demultiplexSummary.txt", row.names=F, sep="\t", quote=F)


###########################
### read count data
umis <- Read10X(data.dir=dge_path, strip.suffix=TRUE)
umis.barcodes <- data.frame(SNo=c(1:length(umis@Dimnames[[2]])), FeatureIndex=umis@Dimnames[[2]])
multiplex.barcodes <- data.frame(FeatureIndex=names(final.calls), SampleIndex=as.vector(final.calls))


###########################
### check for missing cell barcodes
missing.barcodes <- setdiff(umis.barcodes$FeatureIndex, multiplex.barcodes$FeatureIndex)
if(length(missing.barcodes) > 0) 
{
	for (i in 1:length(missing.barcodes))
	{
		multiplex.barcodes[nrow(multiplex.barcodes)+i,] <- c(missing.barcodes[i], "Negative")
	}
}

multiplex.barcodes.merged <- merge(multiplex.barcodes, barcodes, all.x=T)
multiplex.barcodes.merged$Samplename[multiplex.barcodes.merged$SampleIndex=="Doublet"] <- "Doublet"
multiplex.barcodes.merged$Samplename[multiplex.barcodes.merged$SampleIndex=="Negative"] <- "Negative"

umis.barcodes.merged <- merge(umis.barcodes, multiplex.barcodes.merged)
umis.barcodes.merged <- umis.barcodes.merged[order(umis.barcodes.merged$SNo),]

## Check order
if(all(umis.barcodes.merged$FeatureIndex==umis@Dimnames[[2]])) 
{
	print("Order looks good")
}


###################################################
### create Seurat Object containing demultiplexed UMI count
myList <- list()
for (i in 1:nrow(barcodes))
{
	tab <- CreateSeuratObject(counts = as.matrix(umis)[,umis.barcodes.merged$FeatureIndex[umis.barcodes.merged$Samplename==barcodes$Samplename[i]]], project = barcodes$Samplename[i])
	myList[[i]] <- tab
}

aggrData <- merge(myList[[1]], c(myList[[2]], myList[[3]], myList[[4]]), add.cell.ids = barcodes$Samplename)


###########################
### save data
save(bar.table, bar.tsne, bar.table.full, final.calls, umis.barcodes.merged, aggrData, file="deMULTIplex_objs.rda")


###########################
### sessionInfo and clean and quit
sessionInfo()
date()

rm(list=ls())


###########################
####quit
q("yes")
