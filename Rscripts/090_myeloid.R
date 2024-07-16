####
# load packages ----
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(Seurat)
# library(RColorBrewer)
# library(CytoTRACE)


####
# set conditions ----
description = "090_myeloid"
npcs = 20
path = paste0("plots/", description, "_npc", as.character(npcs))
dir.create(path)
RDSfile = paste0("RDSfiles/seu_", description, "_npc", as.character(npcs), ".RDS")
pal = DiscretePalette(n = 36, palette = "polychrome", shuffle = FALSE)

####
# load data ----
seu <- readRDS(file = "RDSfiles/seu_040_cellgroup_annotation.RDS")

# subset
seu <- subset(seu, subset = cellgroup == "myeloid")
seu <- subset(seu, 
              subset = CD3D > 0 | CD3E > 0 | CD3G > 0 | THY1 > 0 | DERL3 > 0 | MS4A1 >0 | EPCAM >0, 
              invert = TRUE)


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
features <- readLines("gene_set/myeloid.txt")
for(i in 1:length(features)){
  tryCatch({
    p <- FeaturePlot(seu, features = features[i], cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) +
      NoAxes() + NoLegend()
    ggsave(paste0(as.character(features[i]), ".png"), plot = p, 
           path = path, 
           width = 5, height = 5, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

FeaturePlot(seu, features = "IL3RA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()


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
ggsave("nanostring_reference.png", path = path, width = 7, height = 5, units = "in", dpi = 150)
Idents(seu) <- "annotation"
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("annotation.png", path = path, width = 9, height = 5, units = "in", dpi = 150)


####
# rename idents ----
Idents(seu) <- "seurat_clusters"
seu <- subset(seu, subset = seurat_clusters %in% c(16, 20, 22), invert = TRUE)
seu <- RenameIdents(
  seu,
  `0` = "M0", `8` = "M0",
  `2` = "Rib-M0", `17` = "Rib-M0",
  `6` = "M1", 
  `5` = "M2", `21` = "M2",
  `4` = "IDA-Mac", `14` = "IDA-Mac", `19` = "IDA-Mac", 
  `11` = "Infl.-Mono", 
  `3` = "Mast", `9` = "Mast", `13` = "Mast", 
  `1` = "Neutro", `7` = "Neutro", `10` = "Neutro",
  `15` = "Eosino", 
  `12` = "DC-CD1C", 
  `18` = "DC-CCL22"
)
seu$celltype <- Idents(seu)
seu$celltype <- factor(seu$celltype, 
                       levels = c("M0", "Rib-M0", "M1", "M2", "IDA-Mac", "Infl.-Mono",
                                  "Mast", "Neutro", "Eosino", "DC-CD1C", "DC-CCL22"))
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


####
# add crude annotation ----
seu <- readRDS("RDSfiles/seu_090_myeloid_npc20.RDS")
levels <- levels(seu$celltype)
levels[levels %in% c("Rib-M0")] <- "M0"
levels[levels %in% c("Infl.-Mono", "IDA-Mac")] <- "M1"
levels[levels %in% c("DC-CD1C", "DC-CCL22")] <- "DC"
seu$celltype_cr <- seu$celltype
levels(seu$celltype_cr) <- levels
Idents(seu) <- "celltype_cr"
DimPlot(seu, label = TRUE, repel = TRUE, cols = pal[30:36]) + NoAxes()
ggsave("celltype_cr.png", path = path, width = 7, height = 5, units = "in", dpi = 150)
saveRDS(seu, file = RDSfile)

DimPlot(seu, 
        # group.by = "celltype",
        cols = "polychrome",
        label = TRUE, repel = TRUE
) + NoAxes()
