####
# description ----

# proceed from lognorm_harmony_npc30, res 0.6


####
# load packages ----
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(RColorBrewer)


####
# set conditions ----
description = "040_cellgroup_annotation"
path = paste0("plots/", description)
dir.create(path)
RDSfile = paste0("RDSfiles/seu_", description, ".RDS")


####
# load data ----
seu <- readRDS(file = "RDSfiles/seu_030_lognorm_harmony_npc30.RDS")
DimPlot(seu)
View(seu[[]])


####
# celltype annotation ----
Idents(seu) <- "seurat_clusters"
seu <- RenameIdents(
  seu,
  `2` = "epithelial", `7` = "epithelial", `9` = "epithelial", `12` = "epithelial", `17` = "epithelial", `20` = "epithelial",
  `1` = "Bcell", `3` = "Bcell", `5` = "Bcell", `11` = "Bcell",
  `0` = "Tcell", `4` = "Tcell", `10` = "Tcell", `13` = "Tcell", `22` = "Tcell",
  `8` = "stromal", `18` = "stromal", `19` = "stromal", `21` = "stromal", 
  `6` = "myeloid", `14` = "myeloid", `15` = "myeloid", 
  `16` = "TandB"
)
seu$cellgroup <- Idents(seu)
seu$cellgroup <- factor(seu$cellgroup, levels = c("epithelial", "stromal", "Tcell", "Bcell", "myeloid", "TandB"))
Idents(seu) <- "cellgroup"
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes() 
ggsave("cellgroup.png", path = path, width = 7, height = 5, units = "in", dpi = 150)
DimPlot(seu, split.by = "exp_group") + NoAxes()
ggsave("split.png", path = path, width = 12, height = 4, units = "in", dpi = 150)


####
# save ----
saveRDS(seu, file = RDSfile)
