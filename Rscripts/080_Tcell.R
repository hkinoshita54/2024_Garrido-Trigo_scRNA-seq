####
# load packages ----
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(Seurat)
# library(CytoTRACE)


####
# set conditions ----
description = "080_Tcell"
npcs = 30
path = paste0("plots/", description, "_npc", as.character(npcs))
dir.create(path)
RDSfile = paste0("RDSfiles/seu_", description, "_npc", as.character(npcs), ".RDS")
pal = DiscretePalette(n = 36, palette = "polychrome", shuffle = FALSE)


####
# load data ----
seu <- readRDS(file = "RDSfiles/seu_040_cellgroup_annotation.RDS")

# subset
seu <- subset(seu, subset = cellgroup %in% c("Tcell", "TandB"))
seu <- subset(seu, subset = DERL3 > 0 | MS4A1 >0 | C1QA > 0 | EPCAM >0 | CD79A > 0, invert = TRUE)


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
features <- readLines("gene_set/Tcell.txt")
for(i in 1:length(features)){
  tryCatch({
    p <- FeaturePlot(seu, features = features[i], cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) +
      NoAxes() + NoLegend()
    ggsave(paste0(as.character(features[i]), ".png"), plot = p, 
           path = path, 
           width = 5, height = 5, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

FeaturePlot(seu, features = "IFNG", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()


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

Idents(seu) <- "nanostring_reference"
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("nanostring_reference.png", path = path, width = 9, height = 5, units = "in", dpi = 150)
Idents(seu) <- "annotation"
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("annotation.png", path = path, width = 12, height = 5, units = "in", dpi = 150)


####
# rename idents ----
Idents(seu) <- "seurat_clusters"
seu <- subset(seu, subset = seurat_clusters %in% c(18, 24), invert = TRUE)
seu <- RenameIdents(
  seu,
  `0` = "CD8-CTL", `15` = "CD8-CTL",
  `4` = "CD8-TRM", `17` = "CD8-TRM", 
  `21` = "CD8-FGFBP2",
  `9` = "CD4-naive", `13` = "CD4-naive", 
  `1` = "CD4-ANXA1", `3` = "CD4-ANXA1", `16` = "CD4-ANXA1",
  `2` = "CD4-CCL20",
  `8` = "Treg", `10` = "Treg",
  `11` = "Tfh", `25` = "Tfh", 
  `12` = "DN-TNF", 
  `26` = "DN-EOMES",
  `19` = "NK", 
  `6` = "gd-IEL",
  `22` = "ILC3", 
  `20` = "Cycl.-T",
  `7` = "Rib-T", `14` = "Rib-T",
  `5` = "MT-T", `23` = "MT-T"
)
seu$celltype <- Idents(seu)
seu$celltype <- factor(seu$celltype, 
                       levels = c("CD8-CTL", "CD8-TRM", "CD8-FGFBP2", 
                                  "CD4-naive", "CD4-ANXA1", "CD4-CCL20", "Treg", "Tfh", 
                                  "DN-TNF", "DN-EOMES", "NK", "gd-IEL", "ILC3", "Cycl.-T", "Rib-T", "MT-T"))
Idents(seu) <- "celltype"
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("celltype.png", path = path, width = 7, height = 5, units = "in", dpi = 150)
DimPlot(seu, split.by = "exp_group") + NoAxes()
ggsave("split.png", path = path, width = 12, height = 4, units = "in", dpi = 150)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("vln_RNA_mt_epi.png", path = path, width = 12, height = 3, units = "in", dpi = 150)

seu$cellgroup <- "Tcell"

####
# save ----
saveRDS(seu, file = RDSfile)


####
# add crude annotation ----
seu <- readRDS("RDSfiles/seu_080_Tcell_npc30_removed31.RDS")
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
seu <- subset(seu, subset = celltype %in% c("Rib-T", "MT-T"), invert = TRUE) # removed undefined cluster (Rib and MT)
levels <- levels(seu$celltype)
levels[levels %in% c("CD8-CTL", "CD8-TRM", "CD8-FGFBP2")] <- "CD8"
levels[levels %in% c("CD4-naive", "CD4-ANXA1", "CD4-CCL20")] <- "CD4"
levels[levels %in% c("DN-TNF", "DN-EOMES")] <- "DN"
# levels <- levels[! levels %in% c("Rib-T", "MT-T")]
seu$celltype_cr <- seu$celltype
levels(seu$celltype_cr) <- levels
Idents(seu) <- "celltype_cr"
DimPlot(seu, label = TRUE, repel = TRUE, cols = pal[21:29]) + NoAxes()
ggsave("celltype_cr.png", path = path, width = 7, height = 5, units = "in", dpi = 150)
saveRDS(seu, file = "RDSfiles/seu_080_Tcell_npc30_removed31_2.RDS")

