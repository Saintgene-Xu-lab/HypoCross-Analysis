`#! /usr/bin/env R
# refine cluster to improve transfer result

suppressMessages({
  # basic
  library(magrittr)
  library(Matrix)
  library(future)
  library(furrr)
  library(patchwork)
  library(pals)
  library(future.apply)
  # library(svglite)
  library(tidyverse)
  library(reticulate)
  
  # specific lib
  library(Seurat)
  # library(ComplexHeatmap)
})
###########################################
setwd("/public/group/weiwu_group/majianyun/data/Hypothalamus/raw.all.snRNA")
options(future.globals.maxSize=+Inf)
# theme set
theme_set(theme_light())
# plan(multicore, workers=16)



############### 2025 update
#### output object for transfer
dir.create("rds/snRNA_SCT_transfer_2025/")
dir.create("rds/final_seu/")


seu = readRDS("/cluster/share/hypo_group/snRNA_20250306_clena_rm_female/all_cell_male_677320_20250305.rds")
# seu = subset(seu, idents=c("GABA", "GLU"), invert=F)
# seu$clu = seu$anno_neu
# fseu = merge(seu, seun)
## fseu = merge(seu, seuf)
seu$clu = seu$annotation
fseu = seu

mouse = readRDS("/public/group/weiwu_group/majianyun/data/Hypothalamus/raw.mouse.snRNA/20221128_czh/rds/mouse_czh_own_doublet_rm.rf1.rds")
snames = intersect(colnames(mouse), colnames(fseu))
mouse = mouse[, snames]
mouse$clu = fseu[[]][snames,"clu"]
mouse = mouse[, mouse$sex == "M"]
saveRDS(mouse, "rds/final_seu/mouse.rds")


mouse = DietSeurat(mouse, counts=F, data=T, scale.data = F, assays = "SCT", dimreducs = F, graphs = F, misc=F)
saveRDS(mouse, "rds/snRNA_SCT_transfer_2025/mouse.rds")



macaque = readRDS("/public/group/weiwu_group/majianyun/data/Hypothalamus/raw.macaque.snRNA/20221214_all/rds/macaque_all_doublet_rm.rf1.rds")
snames = intersect(colnames(macaque), colnames(fseu))
macaque = macaque[, snames]
macaque$clu = fseu[[]][snames,"clu"]
saveRDS(macaque, "rds/final_seu/macaque.rds")

macaque = DietSeurat(macaque, counts=F, data=T, scale.data = F, assays = "SCT", dimreducs = F, graphs = F, misc=F)
saveRDS(macaque, "rds/snRNA_SCT_transfer_2025/macaque.rds")


marmoset = readRDS("/public/group/weiwu_group/majianyun/data/Hypothalamus/raw.marmoset.snRNA/20221201_all/rds/marmoset_all_doublet_rm.rf1.rds")
snames = intersect(colnames(marmoset), colnames(fseu))
marmoset = marmoset[, snames]
marmoset$clu = fseu[[]][snames,"clu"]
saveRDS(marmoset, "rds/final_seu/marmoset.rds")


marmoset = DietSeurat(marmoset, counts=F, data=T, scale.data = F, assays = "SCT", dimreducs = F, graphs = F, misc=F)
saveRDS(marmoset, "rds/snRNA_SCT_transfer_2025/marmoset.rds")



# #### inspect OD&MG marker
# om = markers |> 
#   filter(str_starts(cluster, "OD|MG")) |>
#   group_by(cluster) |>
#   filter(p_val_adj < 0.01 & avg_log2FC > 1.5)|>
#   ungroup() |>
#   group_by(gene) |>
#   mutate(n=n()) |>
#   ungroup() |>
#   filter(n==1)
# write_csv(om, "om.csv")

