#! /usr/bin/env R

# combine data and clean

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
  library(SoupX)
  library(harmony)
  # library(SeuratWrappers)
  library(DoubletFinder)
})
###########################################
setwd("/public/group/weiwu_group/majianyun/data/Hypothalamus/raw.mouse.snRNA")
source("/public/group/weiwu_group/majianyun/data/Hypothalamus/preprocess_lib.r")
options(future.globals.maxSize=+Inf)
# theme set
theme_set(theme_light())
###########################################
plan(multicore, workers=4)

# data path
dpath = paste0(list.dirs("./20221128_czh/mouse_all", 
                                recursive = F))
dpath2 = list.dirs("/public/group/weiwu_group/majianyun/data/Hypothalamus/raw.mouse.snRNA/own_mouse/mtx", recursive=F)
data_folder = c(rep("20221128_czh", length(dpath)), rep("own_mouse", length(dpath2)))
sample_name = c(
  str_extract(dpath, pattern="R\\d{4}"),
  basename(dpath2)
)
batch_l1 = sample_name
batch_l2 = c(rep("wqf", length(dpath)), rep("hl", length(dpath2)))
dpath = c(dpath,dpath2)

# generate Seurat list
seul = future_lapply(1:length(dpath), function(p) {
  message("##### ", sample_name[p], " #####")
  seu = Read10X(dpath[p], gene.column=1) %>%
          CreateSeuratObject(project=sample_name[p]) %>%
          PercentageFeatureSet(pattern="^mt-", col.name="mt") 
  seu$batch = batch_l1[p]
  seu$batch2 = batch_l2[p]
  seu$specie = "mouse"
  seu$cf_ratio = exp(log(seu$nCount_RNA) - log(seu$nFeature_RNA))
  
  # standard processing for mouse
  seu = subset(seu, subset=nFeature_RNA > 800 & mt < 5 & cf_ratio > 1.2)
  seu %<>% 
        SCTransform(vars.to.regress="mt", 
                    method="glmGamPoi",
                    variable.features.n=1e4
                    ) %>%
        RunPCA(npcs=100) %>%
        RunUMAP(dims=1:100) %>%
        FindNeighbors(dims=1:100) %>%
        FindClusters()
  
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
  

  # sex
  sex = str_extract(dpath[p], '.$')
  if (sex %in% as.character(3:7)) {
    sex = "M"
  }
  nseu$sex = sex

  # sample area
  sample_area = str_extract(dpath[p], "R\\d{4}_.*-" )
  sample_area = ifelse(is.na(sample_area), "_HY-", sample_area)
  sample_area = str_remove(sample_area, '-')
  sample_area = str_remove(sample_area, '.*_')
  nseu$sample_area = sample_area

  # new name
  nseu = RenameCells(nseu, add.cell.id=sample_name[p])

  # plot
  insp_seu(nseu, data_folder[p], sample_name[p])
  
  # return
  nseu
})
names(seul) = sample_name

cseu = merge(seul[[1]], y = as.vector(seul[-1]),
          project = "mouse.hypothalamus")
gc()
plan(sequential)
cseu %<>% 
  SCTransform(vars.to.regress="mt", method="glmGamPoi", 
                    variable.features.n=1e4) %>%
  RunPCA(npcs=100) %>%
  RunUMAP(dims=1:100) %>%
  FindNeighbors(dims=1:100) %>%
  FindClusters()
gc()
insp_seu_all(cseu, "20221128_czh", "mouse_czh_own", batch_l1)

ccseu = cseu[, cseu$DFRes != "Doublet"]
insp_seu_all(ccseu, "20221128_czh", "mouse_czh_own_doublet_rm", batch_l1)

# Idents(ccseu) = "no"
# cccseu = subset(ccseu, downsample=5e4)
# saveRDS(cccseu, "cccseu.rds", compress="xz")



### based on czh's decisiion
ccseu = readRDS("/public/group/weiwu_group/majianyun/data/Hypothalamus/raw.mouse.snRNA/20221128_czh/rds/mouse_czh_own_doublet_rm.rds")
ncseu = ccseu[, !(ccseu$seurat_clusters %in% c(4, 53:55))]

ncseu  %<>% 
  RunPCA(npcs=100) %>%
  RunUMAP(dims=1:100) %>%
  FindNeighbors(dims=1:100) %>%
  FindClusters()
saveRDS(ncseu, "20221128_czh/rds/mouse_czh_own_doublet_rm.rf1.rds")
