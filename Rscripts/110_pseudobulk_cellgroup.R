# load packages ----
library(tidyverse)
library(Seurat)
library(DESeq2)

# create directories ----
title <- "110_pseudobulk_cellgroup"
dir.create("results")
path <- file.path("results", title)
dir.create(path)


# load data ----
seu <- readRDS(file = "RDSfiles/seu_101_combine_npc30.RDS")
seu <- JoinLayers(seu)

# pseudobulk ----
bulk <- AggregateExpression(seu, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("orig.ident", "cellgroup"))
tail(Cells(bulk))
bulk$sample_ID <- sapply(strsplit(Cells(bulk), split = "_"), "[", 1)
bulk$cellgroup <- sapply(strsplit(Cells(bulk), split = "_"), "[", 2)
bulk$exp_group <- substring(bulk$sample_ID, 1, 2)
bulk_all <- bulk

# subset ----
subset <- "epithelial"
bulk <- subset(bulk_all, subset = cellgroup == subset)

cts <- as.matrix(bulk[["RNA"]]$counts)
coldata <- cbind(colnames(cts), conditions = bulk$exp_group) %>% as.data.frame()
coldata$conditions <- factor(coldata$conditions, levels = c("HC", "UC", "CD"))
rownames(coldata) <- coldata$V1

dds <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata, design = ~conditions)

smallestGroupSize <- 6
keep <- rowSums(counts(dds) >= 3) >= smallestGroupSize
dds <- dds[keep,]

dds <- DESeq(dds)
resultsNames(dds)

vst = varianceStabilizingTransformation(dds)
pcaData <- plotPCA(vst, intgroup = "conditions", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=conditions)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()

exp_group <- "UC"
cont_group <- "HC"
description <- paste(subset, exp_group, "vs", cont_group, sep = "_")
res <- results(dds, contrast = c("conditions", exp_group, cont_group))
res <- data.frame(res) %>% rownames_to_column(var = "gene_name")
write_tsv(res, paste0(path, "//", description, ".tsv"))
write.csv(res$gene_name, paste0(path, "//", subset, "_gene_name.csv"))
