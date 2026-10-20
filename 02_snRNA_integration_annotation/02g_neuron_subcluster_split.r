library(Seurat)
library(tidyverse)
library(reshape2)
library(stringr)


##### load data 
neu <- readRDS("")
neu2 <- CreateSeuratObject(counts = neu@assays$RNA@counts,meta = neu@meta.data)
rm(neu)
######
path <- ""

Idents(neu2) <- neu2$anno6
cells <- levels(neu2)
data_sot <- c()
for (i in 1:length(cells)){
  print(paste0("get:",cells[i]))
  data_sot <- subset(neu2,idents = cells[i])
  saveRDS(data_sot,file = paste0(path,"/",cells[i],".rds"))
  print(paste0("save:",cells[i],"done !!"))
  data_sot <- c()
}
rm(list=ls())

