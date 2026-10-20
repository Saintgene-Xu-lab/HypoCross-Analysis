setwd("")
if(T){
  library(SeuratWrappers)
  library(Seurat)
  library(BPCells)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(Seurat)
  library(SeuratData)
  library(Azimuth)
  library(ggplot2)
  library(patchwork)
  library(COSG)
  # set this option when analyzing large datasets
  options(future.globals.maxSize = 3e+09)
}

### function
seurat_int_sct <- function(obj,int_method,group){
  #obj[["RNA"]] <- as(obj[["RNA"]],Class = "Assay")
  obj$int <- obj@meta.data[,group]
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$int)
  obj <- SCTransform(obj,variable.features.n = 3000)
  obj <- RunPCA(obj,npcs = 50)
  obj <- RunUMAP(obj, dims = 1:50)
  if (int_method == "CCA"){
    obj <- IntegrateLayers(
      object = obj, method = CCAIntegration, normalization.method = "SCT",
      orig.reduction = "pca", new.reduction = "integrated.cca",
      verbose = FALSE,dims = 1:50) 
    obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:50)
    obj <- FindClusters(obj, resolution = 1, cluster.name = "cca_clusters")
    obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:50, reduction.name = "umap")
    
  }
  if (int_method == "RPCA"){
    obj <- IntegrateLayers(
      object = obj, method = RPCAIntegration, normalization.method = "SCT",
      orig.reduction = "pca", new.reduction = "integrated.rpca",
      verbose = FALSE,dims = 1:50
    )
    obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:50)
    obj <- FindClusters(obj, resolution = 1, cluster.name = "prca_clusters")
    obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:50, reduction.name = "umap")
    
  }
  if (int_method == "harmony"){
    obj <- IntegrateLayers(
      object = obj, method = HarmonyIntegration, normalization.method = "SCT",
      orig.reduction = "pca", new.reduction = "harmony",
      verbose = FALSE,dims = 1:50
      
    )
    obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:50)
    obj <- FindClusters(obj, resolution = 1, cluster.name = "harmony_clusters")
    obj <- RunUMAP(obj, reduction = "harmony", dims = 1:50, reduction.name = "umap")
    
    
  }
  if (int_method == "mnn"){
    obj <- IntegrateLayers(
      object = obj, method = FastMNNIntegration,
      new.reduction = "integrated.mnn", normalization.method = "SCT",
      verbose = FALSE,dims = 1:50)
    obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:50)
    obj <- FindClusters(obj, resolution = 1, cluster.name = "mnn_clusters")
    obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:50, reduction.name = "umap")
    
  }
  return(obj)
}

### int
neu <- readRDS("neu_subset.rds")
neu_int_cca <- seurat_int_sct(neu,int_method = "CCA",group = "sub_batch")
neu_int_rpca <- seurat_int_sct(neu,int_method = "RPCA",group = "sub_batch")
neu_int_har <- seurat_int_sct(neu,int_method = "harmony",group = "sub_batch")
neu_int_mmn <- seurat_int_sct(neu,int_method = "mnn",group = "sub_batch")

saveRDS(neu_int_cca,file = paste0(path,"/neu_int_cca.rds"))
saveRDS(neu_int_rpca,file = paste0(path,"/neu_int_rpca.rds"))
saveRDS(neu_int_har,file = paste0(path,"/neu_int_har.rds"))
saveRDS(neu_int_mmn,file = paste0(path,"/neu_int_mmn.rds"))
rm(list=ls())


