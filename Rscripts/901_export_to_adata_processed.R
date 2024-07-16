####
# load packages ----
library(tidyverse)
library(Seurat)


####################
# convert Seurat object to anndata manually following the tutorial below ----
# https://smorabit.github.io/blog/2021/velocyto/

# save metadata table
dir.create("out/processed")
seu <- readRDS("RDSfiles/seu_101_combine_npc30.RDS")
seu <- JoinLayers(seu)
seu@meta.data <- seu@meta.data[,c(1:7,14:17)]
seu$barcode <- colnames(seu)
seu$UMAP_1 <- seu@reductions$umap@cell.embeddings[,1]
seu$UMAP_2 <- seu@reductions$umap@cell.embeddings[,2]
write.csv(seu@meta.data, file='out/processed/seu_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- LayerData(seu, assay = 'RNA', layer = 'counts')
writeMM(counts_matrix, file = 'out/processed/seu_counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seu@reductions$pca@cell.embeddings, file = 'out/processed/seu_pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene' = rownames(counts_matrix)), file = 'out/processed/seu_gene_names.csv',
  quote = F, row.names = F, col.names = F
)
