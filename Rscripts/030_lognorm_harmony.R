####
# load packages ----
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(Seurat)


####
# set conditions ----
description = "030_lognorm_harmony"
npcs = 30
path = paste0("plots/", description, "_npc", as.character(npcs))
dir.create(path)
RDSfile = paste0("RDSfiles/seu_", description, "_npc", as.character(npcs), ".RDS")


####
# load data ----
seu <- readRDS(file = "RDSfiles/seu_020_filtered.RDS")


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
seu <- FindClusters(seu, reduction = "harmony", resolution = 0.6, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("cluster.png", path = path, width = 5, height = 5, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("id.png", path = path, width = 5, height = 5, units = "in", dpi = 150)
DimPlot(seu, group.by = "exp_group") + NoAxes()
ggsave("exp_group.png", path = path, width = 5, height = 5, units = "in", dpi = 150)

Idents(seu) <- "nanostring_reference"
DimPlot(seu) + NoAxes()

####
# save ----
saveRDS(seu, file = RDSfile)

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



