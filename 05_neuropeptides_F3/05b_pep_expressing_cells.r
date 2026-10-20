
if(T){
  library(SeuratWrappers)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ggplot2)
  library(patchwork)
  library(COSG)
  library(qs)
  options(future.globals.maxSize = 3e+09)
  options(Seurat.object.assay.version = "v5")
}


load("regulatory_database.Rdata")
neu <- qread("")
####　
peti_use2 <- qread("p_pep_stat_use_all_spe.qs")
peti_use <- unique(as.character(peti_use2$features.plot))
DefaultAssay(neu) <- "RNA"
neu <- NormalizeData(neu)

Idents(neu) <- neu$anno6
neu <- SplitObject(neu,split.by = "specie")
neu[["mouse"]] <- subset(neu[["mouse"]],cells = rownames(neu[["mouse"]]@meta.data[neu[["mouse"]]$sample_area %in% "HY",]))

path <- ""

#### mouse id 
mou_id <- list()
for (i in 1:length(peti_use)) {

  result <- tryCatch({
    
    p <- VlnPlot(neu[["mouse"]],features = peti_use[i],idents = as.character(peti_use2[peti_use2$features.plot %in% peti_use[i],]$id) )
    p <- p$data
    colnames(p) <- c("value","idents")
    cells <- rownames(p[p$value > 0.1,])
    mou_id[[peti_use[i]]] <- cells
    print(paste0("Done",peti_use[i]))
    
  }, error = function(e) {
    
    i= 1+1
    p <- VlnPlot(neu[["mouse"]],features = peti_use[i],idents = as.character(peti_use2[peti_use2$features.plot %in% peti_use[i],]$id) )
    p <- p$data
    colnames(p) <- c("value","idents")
    cells <- rownames(p[p$value > 0.1,])
    mou_id[[peti_use[i]]] <- cells
    
  })
  
  print(result)


  
}
qsave(mou_id,file = paste0(path,"mou_pet_id.qs"))
print(paste0("Done  mou !!!!!!!!!"))

### mar id 
mar_id <- list()
for (i in 1:length(peti_use)) {
  
  result <- tryCatch({
    
    p <- VlnPlot(neu[["marmoset"]],features = peti_use[i],idents = as.character(peti_use2[peti_use2$features.plot %in% peti_use[i],]$id) )
    p <- p$data
    colnames(p) <- c("value","idents")
    cells <- rownames(p[p$value > 0.1,])
    mar_id[[peti_use[i]]] <- cells
    print(paste0("Done",peti_use[i]))
    
  }, error = function(e) {
    
    i= 1+1
    p <- VlnPlot(neu[["marmoset"]],features = peti_use[i],idents = as.character(peti_use2[peti_use2$features.plot %in% peti_use[i],]$id) )
    p <- p$data
    colnames(p) <- c("value","idents")
    cells <- rownames(p[p$value > 0.1,])
    mar_id[[peti_use[i]]] <- cells
    print(paste0("Done",peti_use[i]))
    
  })
  
  print(result)
  
  
  
}
qsave(mar_id,file = paste0(path,"mar_pet_id.qs"))
print(paste0("Done  mar !!!!!!!!!"))


### mac id 
mac_id <- list()
for (i in 1:length(peti_use)) {
  
  result <- tryCatch({
    
    p <- VlnPlot(neu[["macaque"]],features = peti_use[i],idents = as.character(peti_use2[peti_use2$features.plot %in% peti_use[i],]$id) )
    p <- p$data
    colnames(p) <- c("value","idents")
    cells <- rownames(p[p$value > 0.1,])
    mac_id[[peti_use[i]]] <- cells
    print(paste0("Done",peti_use[i]))
    
  }, error = function(e) {
    
    i= 1+1
    p <- VlnPlot(neu[["macaque"]],features = peti_use[i],idents = as.character(peti_use2[peti_use2$features.plot %in% peti_use[i],]$id) )
    p <- p$data
    colnames(p) <- c("value","idents")
    cells <- rownames(p[p$value > 0.1,])
    mac_id[[peti_use[i]]] <- cells
    print(paste0("Done",peti_use[i]))
    
  })
  
  print(result)
  
  
  
}
qsave(mac_id,file = paste0(path,"mac_pet_id.qs"))
print(paste0("Done  mc !!!!!!!!!"))
rm(list = ls())

