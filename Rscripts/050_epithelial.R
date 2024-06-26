####
# load packages ----
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(Seurat)
# library(CytoTRACE)


####
# set conditions ----
description = "050_epithelial"
npcs = 30
path = paste0("plots/", description, "_npc", as.character(npcs))
dir.create(path)
RDSfile = paste0("RDSfiles/seu_", description, "_npc", as.character(npcs), ".RDS")


####
# load data ----
seu <- readRDS(file = "RDSfiles/seu_040_cellgroup_annotation.RDS")

# subset
seu <- subset(seu, subset = cellgroup == "epithelial")
seu <- subset(seu, subset = CD3D > 0 | CD3E > 0 | CD3G > 0 | C1QA > 0 | DERL3 > 0 | MS4A1 >0, invert = TRUE)



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


####
# feature plots ----
# files <- list.files(path = "gene_set/", full.names = TRUE)
# features <- lapply(files, FUN = readLines) %>% unlist()
features <- readLines("gene_set/epi.txt")
for(i in 1:length(features)){
  tryCatch({
    p <- FeaturePlot(seu, features = features[i], cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) +
      NoAxes() + NoLegend()
    ggsave(paste0(as.character(features[i]), ".png"), plot = p, 
           path = path, 
           width = 5, height = 5, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

FeaturePlot(seu, features = "CD3G", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()


####
# CytoTRACE ----
# dir.create(paste0(path, "/CytoTRACE"))
# counts_matrix <- LayerData(seu, assay='RNA', layer='counts') %>% as.data.frame()
# obj_cell_type_anno <- as.data.frame(seu$seurat_clusters)
# results <- CytoTRACE(counts_matrix, ncores = 4)
# pheno <- as.character(seu$seurat_clusters)
# names(pheno) <- colnames(seu)
# plotCytoTRACE(results, phenotype = pheno, outputDir = paste0(path, "/CytoTRACE/"))
# Get umap embeddings from Seurat...
# plotCytoTRACE(results, phenotype = pheno, gene = "TIMP1", emb = umapEmb)


####
# check markers ----
seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25)
markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 10) %>% ungroup() -> top10
# c9markers <- FindMarkers(seu, ident.1 = 9, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
# c9vsc0c12c16markers <- FindMarkers(seu, ident.1 = 9, ident.2 = c(0, 12, 16), logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)


Idents(seu) <- "seurat_clusters"
# seu <- subset(seu, subset = seurat_clusters %in% c(13, 14), invert = TRUE)
seu <- RenameIdents(
  seu,
  `0` = "Colono", `7` = "Colono", `10` = "Colono", 
  `3` = "PLGC-C",
  `8` = "LAM-C",
  `13` = "BEST4-C",
  `2` = "Rib-C", `4` = "Rib-C", 
  `6` = "TA",
  `5` = "Sec-Pro", `12` = "Sec-Pro",
  `1` = "Goblet", `11` = "Mature-G",
  `9` = "Tuft", `14` = "Paneth", `15` = "EE"
)
seu$celltype <- Idents(seu)
seu$celltype <- factor(seu$celltype, 
                       levels = c("TA", "Rib-C", "LAM-C", "PLGC-C", "Colono", "BEST4-C", 
                                  "Sec-Pro", "Goblet", "Mature-G", "Paneth", "EE", "Tuft"))
Idents(seu) <- "celltype"
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("celltype.png", path = path, width = 7, height = 5, units = "in", dpi = 150)
DimPlot(seu, split.by = "exp_group") + NoAxes()
ggsave("split.png", path = path, width = 12, height = 4, units = "in", dpi = 150)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("vln_RNA_mt_epi.png", path = path, width = 12, height = 3, units = "in", dpi = 150)


####
# save ----
saveRDS(seu, file = RDSfile)


Idents(seu) <- "annotation"
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("annotation.png", path = path, width = 7, height = 5, units = "in", dpi = 150)


####
# add crude annotation ----
seu <- readRDS("RDSfiles/seu_050_epithelial_npc30.RDS")
levels <- levels(seu$celltype)
levels[levels %in% c("Rib-C", "LAM-C", "PLGC-C")] <- "Colono"
levels[levels == "Mature-G"] <- "Goblet"
seu$celltype_cr <- seu$celltype
levels(seu$celltype_cr) <- levels
Idents(seu) <- "celltype_cr"
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("celltype_cr.png", path = path, width = 7, height = 5, units = "in", dpi = 150)
saveRDS(seu, file = RDSfile)
