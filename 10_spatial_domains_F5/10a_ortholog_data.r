### conda activate R_4.3
if(T){
  library(Seurat)
  library(qs)
  library(Seurat)
  library(qs)
  library(tidyverse)
  library(reshape2)
  library(ggplot2)
}
path <- ""
mac <- readRDS(paste0(path,"macaque.rds"))
DefaultAssay(mac) <- "RNA"
mac$specie <- "mac"
mar <- readRDS(paste0(path,"marmoset.rds"))
DefaultAssay(mar) <- "RNA"
mar$specie <- "mar"
mou <- readRDS(paste0(path,"mouse.rds"))
DefaultAssay(mou) <- "RNA"
mou$specie <- "mou"

mac[["RNA"]] <- as(mac[["RNA"]],Class = "Assay")
mar[["RNA"]] <- as(mar[["RNA"]],Class = "Assay")
mou[["RNA"]] <- as(mou[["RNA"]],Class = "Assay")

path <- ""
gene <- read.csv("conv_tb.gprofiler.20230729.csv")
print(paste0("Doing orh data......"))
count <- mar@assays$RNA@counts[intersect(rownames(mar),gene$marmosetGene),]
mar_orh <- CreateSeuratObject(counts = count,meta.data = mar@meta.data)
mar_orh$specie <- "marmoset"
mar_orh[["RNA"]] <- as(mar_orh[["RNA"]],Class = "Assay")
a <- colSums(mar_orh@assays$RNA@counts)
mar_orh <- subset(mar_orh,cells = names(which(a == 0)),invert =T)
qsave(mar_orh,file = paste0(path,"mar_orh.qs"))

count <- mac@assays$RNA@counts[intersect(rownames(mac),gene$macaqueGene),]
mac_orh <- CreateSeuratObject(counts = count,meta.data = mac@meta.data)
mac_orh$specie <- "macaque"
mac_orh[["RNA"]] <- as(mac_orh[["RNA"]],Class = "Assay")
a <- colSums(mac_orh@assays$RNA@counts)
mac_orh <- subset(mac_orh,cells = names(which(a == 0)),invert =T)
qsave(mac_orh,file = paste0(path,"mac_orh.qs"))
print(paste0("Done orh  data......"))

gene_mou <- gene
rownames(gene_mou) <- gene$mouseGene
count <- mou@assays$RNA@counts[intersect(rownames(mou),gene$mouseGene),]
gene_mou <- gene_mou[rownames(count),]
rownames(count) <- gene_mou$marmosetGene
mou_orh <- CreateSeuratObject(counts = count,meta.data = mou@meta.data)
mou_orh$specie <- "mouse"
mou_orh[["RNA"]] <- as(mou_orh[["RNA"]],Class = "Assay")
a <- colSums(mou_orh@assays$RNA@counts)
mou_orh <- subset(mou_orh,cells = names(which(a == 0)),invert =T)
qsave(mou_orh,file = paste0(path,"mou_orh.qs"))

print(paste0("Merge specie data......"))
mac_mar_mou_orh <-  merge(mac_orh, y = c(mar_orh, mou_orh),add.cell.ids = c("mac", "mar", "mou")) # 
qsave(mac_mar_mou_orh,file = paste0(path,"three_specie_merge_data.qs"))

rm(list=ls())


