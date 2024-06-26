####
# load packages ----
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(Seurat)


####
# set conditions ----
description = "020_QC_filter"
# npcs = 30
path = paste0("plots/", description)
dir.create(path)
RDSfile = paste0("RDSfiles/seu_", description, ".RDS")


####
# check the annotation from the paper ----
anno <- read.delim(file = "data/GSE214695_cell_annotation.csv.gz", sep = ",")
anno$cellnames <- paste(anno$sample, anno$cell_id, sep = "_")
anno <- distinct(anno, cellnames, .keep_all = TRUE) %>% 
  select(cellnames, annotation, nanostring_reference)


####
# load data and filter by doublet calls, nFeature_RNA and percent.mt ----
seu <- readRDS(file = "RDSfiles/seu_010_unfiltered_DoubletFinder.RDS")
View(seu[[]])
seu@meta.data <- seu@meta.data[,c(1:7,11)]

# add annotation to metadata
meta <- seu[[]] %>% rownames_to_column(var = "cellnames")
meta_anno <- left_join(meta, anno, by = "cellnames") %>% column_to_rownames(var = "cellnames")
seu[[]] <- meta_anno

Idents(seu) <- "orig.ident"
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("vln_RNA_mt_unfiltered.png", path = path, width = 12, height = 3, units = "in", dpi = 150)

seu <- subset(seu, subset = doublet_finder == "Singlet")
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("vln_RNA_mt_singlet.png", path = path, width = 12, height = 3, units = "in", dpi = 150)

# filter mt<50 etc.
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 50)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("vln_RNA_filtered.png", path = path, width = 12, height = 3, units = "in", dpi = 150)

# select only annotated cells
# seu <- seu[ , rownames(anno)]
# anno <- anno[colnames(seu), ]
# seu$annotation <- anno$annotation
# seu$nanostring_reference <- anno$nanostring_reference
# VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
# ggsave("vln_RNA_annotated.png", path = path, width = 12, height = 3, units = "in", dpi = 150)


saveRDS(seu, file = "RDSfiles/seu_020_filtered.RDS")
