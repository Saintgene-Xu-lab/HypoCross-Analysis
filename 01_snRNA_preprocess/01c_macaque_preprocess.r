#! /usr/bin/env R

# macaque SoupX data preprocessing

suppressMessages({
  # basic
  library(magrittr)
  library(Seurat)
  library(Matrix)
  library(future)
  library(furrr)
  library(patchwork)
  library(pals)
  library(future.apply)
  library(svglite)
  library(tidyverse)
  
  # specific lib
  library(DoubletFinder)
})
###########################################
setwd("/public/group/weiwu_group/majianyun/data/Hypothalamus/raw.macaque.snRNA")
source("/public/group/weiwu_group/majianyun/data/Hypothalamus/preprocess_lib.r")
options(future.globals.maxSize=+Inf)
# theme set
theme_set(theme_light())
###########################################
plan(multicore, workers=8)

# data path
dpath = list.dirs("./20221214_all/macaque_all", recursive = F)
sample_name = basename(dpath)
data_folder = "./20221214_all"
batch_l1 = sample_name

mt_genes = c("ND6","COX3","COX1","ND5","ND4","ND2","ND4L","ATP8","CYTB","COX2","ND3","ATP6","ND1")

# load mj data
mjrds %<-% readRDS("macaque.snRNA.mj.rds")
  
# generate Seurat list
seul = future_lapply(1:length(dpath), function(p) {
  message("##### ", sample_name[p], " #####")
  seu = Read10X(dpath[p], gene.column=1) %>%
          CreateSeuratObject(project=sample_name[p]) 
  mt_genes_sub = mt_genes[mt_genes %in% rownames(seu)]
  seu = PercentageFeatureSet(seu, features=mt_genes_sub, col.name="mt")
  seu$batch = batch_l1[p]
  seu$batch2 = "HY"
  seu$specie = "macaque"
  seu$cf_ratio = exp(log(seu$nCount_RNA) - log(seu$nFeature_RNA))
  
  # standard processing for macaque 
  seu = subset(seu, subset=nFeature_RNA > 800 & mt < 5 & cf_ratio > 1.2)
  seu %<>% 
        SCTransform(vars.to.regress="mt", 
                    method="glmGamPoi"
                    ) %>%
        RunPCA(npcs=50) %>%
        RunUMAP(dims=1:50) %>%
        FindNeighbors(dims=1:50) %>%
        FindClusters()
  
  # add label
  seu$sex = "M"
  seu$sample_area = "HY"
  nseu = doublet_pipe(seu)
  # new name
  nseu = RenameCells(nseu, add.cell.id=sample_name[p])

  # plot
  insp_seu(nseu, data_folder, sample_name[p])
  
  # return
  nseu
}, future.seed=42)
names(seul) = sample_name

# add & remove label
mjrds$brain_area = NULL
mjrds$anno_MJ = mjrds$anno
mjrds$anno = NULL
mjrds$batch = "MJ"
mjrds$batch2 = "HY"
mjrds$specie = "macaque"
mjrds$cf_ratio = exp(log(mjrds$nCount_RNA) - log(mjrds$nFeature_RNA))
mjrds$sex = "M"
mjrds$sample_area = "HY"
mt_genes_sub = mt_genes[mt_genes %in% rownames(mjrds)]
DefaultAssay(mjrds) = "RNA"
mjrds = PercentageFeatureSet(mjrds, features=mt_genes_sub, col.name="mt")
# deal with data
# standard processing for macaque 
mjrds = subset(mjrds, subset=nFeature_RNA > 800 & mt < 5 & cf_ratio > 1.2)
mjrds %<>% 
      SCTransform(vars.to.regress="mt", 
                  method="glmGamPoi"
                  ) %>%
      RunPCA(npcs=50) %>%
      RunUMAP(dims=1:50) %>%
      FindNeighbors(dims=1:50) %>%
      FindClusters() 
mjrds = RenameCells(mjrds, add.cell.id="MJ")
mjrds = doublet_pipe(mjrds)

# append
seul = append(seul, mjrds)
batch_l1 = names(seul) = c(sample_name, "MJ")
 

cseu = merge(seul[[1]], y = as.vector(seul[-1]),
          project = "macaque.hypothalamus")

gc()
plan(sequential)
cseu %<>% 
  SCTransform(vars.to.regress="mt", method="glmGamPoi", variable.features.n=1e4) %>%
  RunPCA(npcs=50) %>%
  RunUMAP(dims=1:50) %>%
  FindNeighbors(dims=1:50) %>%
  FindClusters()
gc()
insp_seu_all(cseu, "20221214_all", "macaque_all", batch_l1)

ccseu = cseu[, cseu$DFRes != "Doublet"]
insp_seu_all(ccseu, "20221214_all", "macaque_all_doublet_rm", batch_l1)

ncseu = ccseu %>% 
  RunPCA(npcs=50) %>%
  RunUMAP(dims=1:50) %>%
  FindNeighbors(dims=1:50) %>%
  FindClusters()
saveRDS(ncseu, "20221214_all/rds/macaque_all_doublet_rm.rf1.rds")


