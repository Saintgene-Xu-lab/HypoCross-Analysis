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
  options(Seurat.object.assay.version = "v5")
}

### function
seuratv5_int_sct <- function(obj,int_method,group){
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
data <- readRDS("raw.rds")
dat_int_cca <- seuratv5_int_sct(data,int_method = "CCA",group = "sub_batch")
dat_int_rpca <- seuratv5_int_sct(data,int_method = "RPCA",group = "sub_batch")
dat_int_har <- seuratv5_int_sct(data,int_method = "harmony",group = "sub_batch")
dat_int_mmn <- seuratv5_int_sct(data,int_method = "mnn",group = "sub_batch")

saveRDS(dat_int_cca,file = paste0(path,"/all_int_cca.rds"))
saveRDS(dat_int_rpca,file = paste0(path,"/all_int_rpca.rds"))
saveRDS(dat_int_har,file = paste0(path,"/all_int_har.rds"))
saveRDS(dat_int_mmn,file = paste0(path,"/all_int_mmn.rds"))

### screan cluster
re <- seq(0.1,2,0.1)
for (i in re) {
  dat_int_rpca <- FindClusters(dat_int_rpca, resolution = re[i], cluster.name = paste0("re_",re[i]))
}
saveRDS(dat_int_rpca,file = paste0(path,"/all_int_rpca.rds"))


rm(list=ls())


