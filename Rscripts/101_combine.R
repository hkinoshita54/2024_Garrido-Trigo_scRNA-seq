####
# load packages ----
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(RColorBrewer)
# library(CytoTRACE)


####
# set conditions ----
description = "101_combine"
npcs = 30
path = paste0("plots/", description, "_npc", as.character(npcs))
dir.create(path)
RDSfile = paste0("RDSfiles/seu_", description, "_npc", as.character(npcs), ".RDS")


####
# load data ----
epi <- readRDS(file = "RDSfiles/seu_050_epithelial_npc30.RDS")
str <- readRDS(file = "RDSfiles/seu_060_stromal_npc20.RDS")
Bcell <- readRDS(file = "RDSfiles/seu_070_Bcell_npc30.RDS")
Tcell <- readRDS(file = "RDSfiles/seu_080_Tcell_npc30_removed31_2.RDS")
mye <- readRDS(file = "RDSfiles/seu_090_myeloid_npc20.RDS")

# cells in both Bcell and Tcell > already removed (in 101_combine.R)
# ambiguous <- intersect(colnames(Bcell), colnames(Tcell))
# Bcell$celltype[ambiguous]    # most of them are cycling
# Tcell$celltype[ambiguous]    # most of them are cycling or MT-
# Tcell <- subset(Tcell, cells = ambiguous, invert = TRUE)    # remove ambiguous cells from Tcells (this is argitrary)
# saveRDS(Tcell, "RDSfiles/seu_080_Tcell_npc30_removed31.RDS")

# merge
seu <- merge(x = epi, y = c(str, Bcell, Tcell, mye))
rm(epi, str, Bcell, Tcell, mye)
table(seu$cellgroup)
seu$cellgroup <- factor(seu$cellgroup, levels = c("epithelial", "stromal", "Bcell", "Tcell", "myeloid"))

####
# clustering with harmony integration ----
seu <-DietSeurat(seu) 
seu <- JoinLayers(seu)
seu[["RNA"]]$data <- NULL
seu[["RNA"]]$scale.data <- NULL
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)
seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
seu <- RunPCA(seu, npcs = npcs)
seu <- IntegrateLayers(
  object = seu, method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, reduction = "harmony", resolution = 2, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("cluster.png", path = path, width = 5, height = 5, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("id.png", path = path, width = 5, height = 5, units = "in", dpi = 150)
DimPlot(seu, group.by = "exp_group") + NoAxes()
ggsave("exp_group.png", path = path, width = 5, height = 5, units = "in", dpi = 150)


####
# feature plots ----
# files <- list.files(path = "gene_set/", full.names = TRUE)
# features <- lapply(files, FUN = readLines) %>% unlist()
features <- readLines("gene_set/cellgroup.txt")
for(i in 1:length(features)){
  tryCatch({
    p <- FeaturePlot(seu, features = features[i], cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) +
      NoAxes() + NoLegend()
    ggsave(paste0(as.character(features[i]), ".png"), plot = p, 
           path = path, 
           width = 5, height = 5, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# FeaturePlot(seu, features = "CRYAB", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()


Idents(seu) <- "nanostring_reference"
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("nanostring_reference.png", path = path, width = 12, height = 5, units = "in", dpi = 150)
Idents(seu) <- "annotation"
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("annotation.png", path = path, width = 18, height = 6, units = "in", dpi = 150)

Idents(seu) <- "cellgroup"
DimPlot(seu, label = TRUE, repel = TRUE, cols = brewer.pal(5, "Set1")) + NoAxes()
ggsave("cellgroup.png", path = path, width = 6, height = 5, units = "in", dpi = 150)
Idents(seu) <- "celltype"
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("celltype.png", path = path, width = 10, height = 5, units = "in", dpi = 150)
Idents(seu) <- "celltype_cr"
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()
ggsave("celltype_cr.png", path = path, width = 10, height = 5, units = "in", dpi = 150)
# DimPlot(seu, split.by = "exp_group") + NoAxes()
# ggsave("split.png", path = path, width = 12, height = 4, units = "in", dpi = 150)
# VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
# ggsave("vln_RNA_mt_epi.png", path = path, width = 12, height = 3, units = "in", dpi = 150)


####
# save ----
saveRDS(seu, file = RDSfile)
