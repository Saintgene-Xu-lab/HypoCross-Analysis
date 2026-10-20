#! /usr/bin/env R

# marmoset SoupX data preprocessing

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
setwd("/public/group/weiwu_group/majianyun/data/Hypothalamus/raw.marmoset.snRNA")
source("/public/group/weiwu_group/majianyun/data/Hypothalamus/preprocess_lib.r")
options(future.globals.maxSize=+Inf)
# theme set
theme_set(theme_light())
###########################################
plan(multicore, workers=4)

# data path
dpath = list.dirs("./20221201_all/marmoset_all/", recursive = F)
sample_name = basename(dpath)
data_folder = "./20221201_all/"
batch_l1 = sample_name

# generate Seurat list
seul = future_lapply(1:length(dpath), function(p) {
  message("##### ", sample_name[p], " #####")
  seu = Read10X(dpath[p], gene.column=1) %>%
          CreateSeuratObject(project=sample_name[p]) 
  mt_genes = c("ND6","COX3","COX1","ND5","ND4","ND2","ND4L","ATP8","CYTB","COX2","ND3","ATP6","ND1")
  mt_genes = mt_genes[mt_genes %in% rownames(seu)]
  seu = PercentageFeatureSet(seu, features=mt_genes, col.name="mt")
  seu$batch = batch_l1[p]
  seu$batch2 = "HY"
  seu$specie = "marmoset"
  seu$cf_ratio = exp(log(seu$nCount_RNA) - log(seu$nFeature_RNA))
  
  # standard processing for marmoset 
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
  
  # doublet identification
  filter_rate = 0.075
  t = paramSweep_v3(seu, 1:10, sct=T, num.cores=1)
  tt = summarizeSweep(t, GT = FALSE)
  ttt = find.pK(tt)
  pk = as.numeric(levels(ttt[which(ttt$BCmetric == max(ttt$BCmetric)), "pK"]))[ttt[which(ttt$BCmetric == max(ttt$BCmetric)), "pK"]]
  pn = tt %>% filter(pK == pk)
  pn = as.numeric(levels(pn[which(pn$BCreal == max(pn$BCreal)), "pN"]))[pn[which(pn$BCreal == max(pn$BCreal)), "pN"]]
  mht = modelHomotypic(seu$seurat_clusters)
  nExp_poi = round(filter_rate*ncol(seu))
  nExp_poi.adj = round(nExp_poi*(1-mht))
  DFparam = paste(0.25, pk, nExp_poi, sep='_')
  nseu = doubletFinder_v3(seu, PCs=1:10, pN=0.25, pK=pk, nExp=nExp_poi, reuse.pANN=FALSE, sct=T)
  df_slot_name = paste("DF.classifications", DFparam, sep='_')
  df_slot_name2 = paste('pANN', DFparam, sep='_')
  nseu[["DFRes"]] = nseu[[df_slot_name]]
  nseu[["DFparam"]] = DFparam
  nseu[[df_slot_name]] = NULL
  nseu[[df_slot_name2]] = NULL
  # nseu = nseu[, nseu[[df_slot_name]] == "Singlet"]
  
  # new name
  nseu = RenameCells(nseu, add.cell.id=sample_name[p])

  # plot
  insp_seu(nseu, data_folder, sample_name[p])
  
  # return
  nseu
}, future.seed=42)
names(seul) = sample_name
cseu = merge(seul[[1]], y = as.vector(seul[-1]),
          project = "marmoset.hypothalamus")

gc()
plan(sequential)
cseu %<>% 
  SCTransform(vars.to.regress="mt", method="glmGamPoi", variable.features.n=1e4) %>%
  RunPCA(npcs=50) %>%
  RunUMAP(dims=1:50) %>%
  FindNeighbors(dims=1:50) %>%
  FindClusters()
gc()
insp_seu_all(cseu, "20221201_all", "marmoset_all", batch_l1)

ccseu = cseu[, cseu$DFRes != "Doublet"]
insp_seu_all(ccseu, "20221201_all", "marmoset_all_doublet_rm", batch_l1)

ncseu = ccseu %>% 
  RunPCA(npcs=50) %>%
  RunUMAP(dims=1:50) %>%
  FindNeighbors(dims=1:50) %>%
  FindClusters()
saveRDS(ncseu, "20221201_all/rds/marmoset_all_doublet_rm.rf1.rds")