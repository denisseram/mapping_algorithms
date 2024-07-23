
source('utils_seurat.R')

suppressPackageStartupMessages({
  library(symphony)
  library(Seurat)
  suppressWarnings({library(SeuratData)})
  library(ggplot2)
  library(dplyr)
  library(magrittr)
  library(Matrix)
  library(sctransform)
})


suppressWarnings({
  #SeuratData::InstallData('hcabm40k')
  SeuratData::LoadData('hcabm40k')    
})

cells_ref <- hcabm40k@meta.data %>% subset(orig.ident %in% paste0('MantonBM', 1:4)) %>% rownames()
cells_query <- hcabm40k@meta.data %>% subset(orig.ident %in% paste0('MantonBM', 5:8)) %>% rownames()
hcabm40k@meta.data %>% head(3)

.verbose <- FALSE
# Run standard Seurat pipeline with log normalization
obj <- Seurat::CreateSeuratObject(hcabm40k@assays$RNA@counts[, cells_ref]) %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = .verbose) %>% 
  RunPCA(verbose = .verbose)  %>% 
  RunHarmony.Seurat('orig.ident', verbose = .verbose) 
  #FindNeighbors(dims = 1:20, reduction = 'harmony', verbose = .verbose) %>%  # previous version of this tutorial was missing reduction argument
  #FindClusters(resolution = 0.5, verbose = .verbose)
  

obj[['umap']] <- RunUMAP2(Embeddings(obj, 'harmony')[, 1:20], 
                          assay='RNA', verbose=FALSE, umap.method='uwot', return.model=TRUE)

options(repr.plot.height = 4, repr.plot.width = 6)
DimPlot(obj, reduction = 'umap', group.by = 'seurat_clusters', shuffle = TRUE)

#make symphony ref object