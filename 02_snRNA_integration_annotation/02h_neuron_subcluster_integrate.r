if(T){
    library(Seurat)
    library(tidyverse)
    library(reshape2)
    library(stringr)
}

### loading data
meta <- readRDS("meta.rds")
cells <- names(table(meta$anno6))
### 

tre_cell <- c("EN26","EN30","IN28","IN30","IN31")
cells <- setdiff(cells,tre_cell)

sct_int_pipline <- function(data){
    data <- SplitObject(data,split.by = "specie")
    data <- lapply(data,function(x){
    x <- SCTransform(x,verbose = FALSE,vst.flavor = "v2") %>%
        RunPCA(npcs = 30, verbose = FALSE)
    })

    features <- SelectIntegrationFeatures(object.list = data, nfeatures = 3000)
    data <- PrepSCTIntegration(object.list = data, anchor.features = features)

    anchors <- FindIntegrationAnchors(object.list = data, normalization.method = "SCT",
                                            anchor.features = features)
    data.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

    DefaultAssay(data.sct) <- "integrated"
    data.sct <- RunPCA(data.sct, verbose = FALSE)
    data.sct <- RunUMAP(data.sct, reduction = "pca", dims = 1:20, verbose = FALSE)
    data.sct <- FindNeighbors(data.sct, reduction = "pca", dims = 1:30)
    data.sct <- FindClusters(data.sct, resolution = 0.1)
    data.sct <- FindClusters(data.sct, resolution = 0.2)
    data.sct <- FindClusters(data.sct, resolution = 0.3)
    data.sct <- FindClusters(data.sct, resolution = 0.4)
    data.sct <- FindClusters(data.sct, resolution = 0.5)
    data.sct <- FindClusters(data.sct, resolution = 0.6)
    data.sct <- FindClusters(data.sct, resolution = 0.7)
    data.sct <- FindClusters(data.sct, resolution = 0.8)
    data.sct <- FindClusters(data.sct, resolution = 0.9)
    data.sct <- FindClusters(data.sct, resolution = 1)
    return(data.sct)
}

path <- "sub_cluster"
path_save <- "sub_cluster_int"
cells
data <- c()
meta <- c()
for(i in 1:length(cells)){
    print(paste0(path,"/",cells[i],".rds"))
    data <- readRDS(paste0(path,"/",cells[i],".rds"))
    DefaultAssay(data) <- "RNA"
    Idents(data) <- data$specie
    print(paste0("runing : ",cells[i]))
    data <- sct_int_pipline(data)
    saveRDS(data,paste0(path_save,"/",cells[i],"_int.rds"))
    meta <- data@meta.data
    saveRDS(meta,paste0(path_save,"/",cells[i],"_meta.rds"))
    print(paste0("Done: ",cells[i]," !!!"))
    data <- c()
    gc()
}

rm(list=ls())
