####
# load packages ----
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(DoubletFinder)


####
# load data ----
sra_table <- read.delim(file = "data/SraRunTable.txt", sep = ",")
GSMids <- sra_table$Sample.Name %>% unique() %>% sort()

ids <- c(paste0("HC", as.character(1:6)), paste0("UC", as.character(1:6)), paste0("CD", as.character(1:6)))    # make character vector of sample ids
exp_group <- c(rep(c("HC", "UC", "CD"), each = 6))    # to add to the meta data
pts_data <- read.delim(file = "data/metadata_scrnaseq.txt", sep = "\t")    # to add to the meta data

seu_list <- list()
for (i in 1:length(ids)){
  mtx <- dir(path = "data/GSE214695_RAW/", pattern = paste0(GSMids[i], ".*", "matrix"), full.names = TRUE)
  features <- dir(path = "data/GSE214695_RAW/", pattern = paste0(GSMids[i], ".*", "features"), full.names = TRUE)
  barcodes <- dir(path = "data/GSE214695_RAW/", pattern = paste0(GSMids[i], ".*", "barcodes"), full.names = TRUE)
  cts <- ReadMtx(mtx = mtx, features = features, cells = barcodes)
  seu <- CreateSeuratObject(counts = cts, project = ids[i], min.cells = 3, min.features = 200)
  seu$exp_group <- exp_group[i]
  seu$age <- pts_data$Age[i]
  seu$sex <- pts_data$Sex[i]
  seu$tissue <- pts_data$Biopsy.surgical.resection.area[i]
  seu_list <- append(seu_list, seu)
  rm(cts, features, barcodes)
}


####
# DoubletFinder ----
for (i in 1:length(seu_list)) {
  # Pre-process seurat object with standard seurat workflow
  seu <- NormalizeData(seu_list[[i]])
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu)
  seu <- RunUMAP(seu, dims = 1:10)
  seu <- FindNeighbors(object = seu, dims = 1:10)              
  seu <- FindClusters(object = seu, resolution = 0.1)
  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep(seu, PCs = 1:10)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  annotations <- seu@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round((0.01 * ncol(seu) / 1000) * nrow(seu@meta.data)) ## doublet formation rate - tailor for your dataset, see below
  # https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  seu <- doubletFinder(seu = seu, 
                          PCs = 1:10, 
                          pK = optimal.pk,
                          nExp = nExp.poi.adj)
  metadata <- seu@meta.data
  colnames(metadata)[ncol(metadata)] <- "doublet_finder"
  seu@meta.data <- metadata 
  
  # save
  seu_list[[i]] <- seu
}


####
# visualize doublet ----
path = "plots/010_load_data_DoubletFinder"
dir.create(path, recursive = TRUE)
for(i in 1:length(seu_list)){
  p <- DimPlot(seu_list[[i]], group.by = "doublet_finder") + NoAxes()
  ggsave(paste0(as.character(ids[i]), "_DoubletFinder.png"), plot = p, 
         path = path, 
         width = 5, height = 5, units = "in", dpi = 150)
  v <- VlnPlot(seu_list[[i]], features = c("nFeature_RNA", "nCount_RNA"), group.by = "doublet_finder")
  ggsave(paste0("vln_", as.character(ids[i]), "_DoubletFinder.png"), plot = v, 
         path = path, 
         width = 5, height = 3, units = "in", dpi = 150)
}


####
# merge the list and save ----
seu <- merge(x = seu_list[[1]], y = seu_list[2:length(ids)], add.cell.ids = ids)
saveRDS(seu, file = "RDSfiles/seu_010_unfiltered_DoubletFinder.RDS")
