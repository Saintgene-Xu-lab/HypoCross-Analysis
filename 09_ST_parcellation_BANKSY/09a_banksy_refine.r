#! /usr/bin/env R

# Banksy refine 
# R4.3 required

suppressMessages({
  # basic
  library(magrittr)
  library(Matrix)
  library(future)
  library(furrr)
  library(patchwork)
  library(pals)
  library(future.apply)
  library(svglite)
  library(tidyverse)
  library(qs)
  
  # specific lib
  library(Seurat)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(scuttle)
  library(Banksy)
  library(presto)
  library(scater)
  library(ComplexHeatmap)
  library(clusteval)
  library(ClusterR)

})
###########################################
setwd("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial")
plan(multicore, workers=30)
options(future.globals.maxSize=+Inf)

source("/public/group/weiwu_group/majianyun/data/Hypothalamus/draw_spatial_lib.r")
#########################################

###### get all chip res
###### load data start
mous = readRDS("/cluster/share/hypo_group/Spatial/spatial_rds_20250328/mouse.rds")
macs = readRDS("/cluster/share/hypo_group/Spatial/spatial_rds_20250328/macaque.rds")
mars = readRDS("/cluster/share/hypo_group/Spatial/spatial_rds_20250328/marmoset.rds")
mars_exclu = mars[, !(mars$chip %in% c("T001_T972", "T002_T971", "T020_T615"))]
mars = mars_exclu


nseul_v2 = lapply(list(mous, macs, mars), function(seu) {
  sl = SplitObject(seu, "chip")
  k = c(15,30)
  lambda = c(0, 0.8)
  res = 1.2
  sll = lapply(sl, function(s) {
    gcm = s@assays$RNA@counts
    locs = as.matrix(s[[c("rx", "ry")]])
    se = SpatialExperiment(assay = list(counts = gcm), spatialCoords = locs)
    se$chip = s$chip[1]
    print(head(se$chip))

    se <- computeLibraryFactors(se)
    aname <- "normcounts"
    assay(se, aname) <- normalizeCounts(se, log = FALSE)

    lambda <- lambda
    res <- res
    k_geom <- k

    se <- Banksy::computeBanksy(se, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)
  })
  se = Reduce(cbind, sll)
  set.seed(42)
  se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda, group = "chip", seed = 42)
  se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda, group = "chip", seed = 42)
  se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = lambda, group = "chip", seed = 42, resolution = res)
  se <- Banksy::connectClusters(se)
  cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(TRUE), lambda, res)
  
  coln = colnames(colData(se))
  print(coln)
  i = str_detect(coln, "clust_M1_lam0\\.8_k50_res.+")
  if (length(i) > 1) i = which(i==T)[1]
  print(coln)
  coln = coln[i]

  se
})
qsave(nseul_v2, "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/3s.bks.qs")


nseul_v3 = lapply(list(macs, mars), function(seu) {
  sl = SplitObject(seu, "chip")
  k = c(15,60)
  lambda = c(0, 0.8)
  res = 1.2
  sll = lapply(sl, function(s) {
    gcm = s@assays$RNA@counts
    locs = as.matrix(s[[c("rx", "ry")]])
    se = SpatialExperiment(assay = list(counts = gcm), spatialCoords = locs)
    se$chip = s$chip[1]
    print(head(se$chip))

    se <- computeLibraryFactors(se)
    aname <- "normcounts"
    assay(se, aname) <- normalizeCounts(se, log = FALSE)

    lambda <- lambda
    res <- res
    k_geom <- k

    se <- Banksy::computeBanksy(se, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)
  })
  se = Reduce(cbind, sll)
  set.seed(42)
  se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda, group = "chip", seed = 42)
  se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda, group = "chip", seed = 42)
  se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = lambda, group = "chip", seed = 42, resolution = res)
  se <- Banksy::connectClusters(se)
  cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(TRUE), lambda, res)
  
  coln = colnames(colData(se))
  print(coln)
  i = str_detect(coln, "clust_M1_lam0\\.8_k50_res.+")
  if (length(i) > 1) i = which(i==T)[1]
  print(coln)
  coln = coln[i]

  se
})

qsave(nseul_v3, "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/3s.bks.v3.qs")


nseul_v6 = lapply(list(mars_exclu), function(seu) {
  sl = SplitObject(seu, "chip")
  k = c(15,60)
  lambda = c(0, 0.8)
  res = 1.2
  sll = lapply(sl, function(s) {
    gcm = s@assays$RNA@counts
    locs = as.matrix(s[[c("rx", "ry")]])
    se = SpatialExperiment(assay = list(counts = gcm), spatialCoords = locs)
    se$chip = s$chip[1]
    print(head(se$chip))

    se <- computeLibraryFactors(se)
    aname <- "normcounts"
    assay(se, aname) <- normalizeCounts(se, log = FALSE)

    lambda <- lambda
    res <- res
    k_geom <- k

    se <- computeBanksy(se, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)
  })
  se = Reduce(cbind, sll)
  set.seed(42)
  se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda, group = "chip", seed = 42)
  se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda, group = "chip", seed = 42)
  se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = lambda, group = "chip", seed = 42, resolution = res)
  se <- Banksy::connectClusters(se)
  cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(TRUE), lambda, res)
  
  coln = colnames(colData(se))
  print(coln)
  i = str_detect(coln, "clust_M1_lam0\\.8_k50_res.+")
  if (length(i) > 1) i = which(i==T)[1]
  print(coln)
  coln = coln[i]

  se
})

qsave(nseul_v6, "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/3s.bks.v6.qs")

nseul_v2 = qread("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/3s.bks.qs")
nseul_v3 = qread("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/3s.bks.v3.qs")
nseul_v6 = qread("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/3s.bks.v6.qs")

mouse_bks = nseul_v2[[1]]
macaque_bks = nseul_v3[[1]]
marmoset_bks = nseul_v6[[1]]

#################################################
#################################################
########## read and plot banksy cluster #########
#################################################
#################################################
#### one to other and pairwise
oto_pw = function(seu, col, pre_fn, pre_path="/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/") {
  Idents(seu) = col
  seu = NormalizeData(seu)
  mks = FindAllMarkers(seu, only.pos=T)
  write_csv(mks, str_c(pre_path, pre_fn, "one_to_other", ".csv.gz"))

  clusters <- unique(Idents(seu))
  pairwise <- combn(clusters, 2)
  results <- list()
  for(i in 1:ncol(pairwise)) {
      markers <- FindMarkers(seu, ident.1 = pairwise[1, i], ident.2 = pairwise[2, i])
      results[[i]] <- markers |> rownames_to_column("gene") |> mutate(ident.1=pairwise[1, i], ident.2=pairwise[2, i])
  }
  res = bind_rows(results) |> as_tibble()
  write_csv(res, str_c(pre_path, pre_fn, "pairwise", ".csv.gz"))


  rres = res |>
    filter(p_val_adj < 1e-8) |>
    filter(!str_starts(gene, "Rps")) |>
    filter(!str_starts(gene, "Rpl")) |>
    filter(!str_starts(gene, "MT-")) |>  
    group_by(ident.1, ident.2) |>
    summarize(np=sum(avg_log2FC > 5), nn=sum(avg_log2FC < -5)) |>
    ungroup() |>
    mutate(res=ifelse(np>5 & nn >5, "Ok", "filter?"))
  write_csv(rres, str_c(pre_path, pre_fn, "pairwise.judge", ".csv"))
}


##### test label
oto_pw(mous, col="clu", pre_fn="mouse.k30")
oto_pw(macs, col="clu", pre_fn="macaque.k60")
oto_pw(mars, col="clu", pre_fn="marmoset.k60")

oto_pw(mous, col="mainAreaID", pre_fn="mouse.mainAreaID")
oto_pw(macs, col="mainAreaID", pre_fn="macaque.mainAreaID")
oto_pw(mars, col="mainAreaID", pre_fn="marmoset.mainAreaID")


###########################################################################
############ load banksy label
###########################################################################
mous_map = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/mouse.k30.csv.gz")
mous_map = setNames(mous_map$label, mous_map$rownames)
mous$clu = mous_map[colnames(mous)]

macs_map = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/macaque.k60.csv.gz")
macs_map = setNames(macs_map$label, macs_map$rownames)
macs$clu = macs_map[colnames(macs)]

mars_map = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/marmoset.k60.csv.gz")
mars_map = setNames(mars_map$label, mars_map$rownames)
mars$clu = mars_map[colnames(mars)]

mouse_anno = read.csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/Banksy-cluster_mouse250722.csv")
macaque_anno = read.csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/Banksy-cluster_macaque250722.csv")
marmoset_anno = read.csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/Banksy-cluster_marmoset250722.csv")

mous$clu2 = as.character(setNames(mouse_anno$Anno, mouse_anno$bks_id)[mous$clu])
mous$clu3 = as.character(setNames(mouse_anno$Simple_Assign, mouse_anno$bks_id)[mous$clu])

macs$clu2 = as.character(setNames(macaque_anno$Anno, macaque_anno$bks_id)[macs$clu])
macs$clu3 = as.character(setNames(macaque_anno$Simple_Assign, macaque_anno$bks_id)[macs$clu])

mars$clu2 = as.character(setNames(marmoset_anno$Anno, marmoset_anno$bks_id)[mars$clu])
mars$clu3 = as.character(setNames(marmoset_anno$Simple_Assign, marmoset_anno$bks_id)[mars$clu])

##################################### mouse all
mous$`clu3.5` = NA
mous_path = "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/bks_clu3.5_mouse/"
area_mks = list()
dir.create(mous_path)
for (area_name in unique(mous$clu3)) {
  area_subset = which(mous$clu3 == area_name)
  mou_pca = mouse_bks@int_colData@listData$reducedDims[["PCA_M1_lam0.8"]][area_subset,]
  km_res = MiniBatchKmeans(mou_pca, clusters=5, batch_size=40, num_init=200, max_iters=100, 
                            init_fraction=0.25, initializer='kmeans++', early_stop_iter=10,
                            verbose=F, seed=42)
  rknn = tibble(cellname=names(area_subset), knn5=predict(object=km_res, newdata=mou_pca))
  mous$`clu3.5` = ifelse(is.na(mous$`clu3.5`), setNames(rknn$knn5, rknn$cellname)[colnames(mous)], mous$`clu3.5`)
}
area_mks = bind_rows(area_mks)
write_csv(area_mks, "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/bks_clu3.5_mouse/mks.csv")
##################################### macaque all
macs$`clu3.5` = NA
macs_path = "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/bks_clu3.5_macaque/"
area_mks = list()
dir.create(macs_path)
for (area_name in unique(macs$clu3)) {
  area_subset = which(macs$clu3 == area_name)
  mac_pca = macaque_bks@int_colData@listData$reducedDims[["PCA_M1_lam0.8"]][area_subset,]
  km_res = MiniBatchKmeans(mac_pca, clusters=5, batch_size=40, num_init=200, max_iters=100, 
                            init_fraction=0.25, initializer='kmeans++', early_stop_iter=10,
                            verbose=F, seed=42)
  rknn = tibble(cellname=names(area_subset), knn5=predict(object=km_res, newdata=mac_pca))
  macs$`clu3.5` = ifelse(is.na(macs$`clu3.5`), setNames(rknn$knn5, rknn$cellname)[colnames(macs)], macs$`clu3.5`)
}
area_mks = bind_rows(area_mks)
write_csv(area_mks, "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/bks_clu3.5_macaque/mks.csv")
##################################### marmoset all
mars$`clu3.5` = NA
mars_path = "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/bks_clu3.5_marmoset/"
area_mks = list()
dir.create(mars_path)
for (area_name in unique(mars$clu3)) {
  area_subset = which(mars$clu3 == area_name)
  mar_pca = marmoset_bks@int_colData@listData$reducedDims[["PCA_M1_lam0.8"]][area_subset,]
  km_res = MiniBatchKmeans(mar_pca, clusters=5, batch_size=40, num_init=200, max_iters=100, 
                            init_fraction=0.25, initializer='kmeans++', early_stop_iter=10,
                            verbose=F, seed=42)
  rknn = tibble(cellname=names(area_subset), knn5=predict(object=km_res, newdata=mar_pca))
  mars$`clu3.5` = ifelse(is.na(mars$`clu3.5`), setNames(rknn$knn5, rknn$cellname)[colnames(mars)], mars$`clu3.5`)
}
area_mks = bind_rows(area_mks)
write_csv(area_mks, "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/bks_clu3.5_marmoset/mks.csv")








################# clu4 add
clu3p5 = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/banksy_cluster_251019.csv") |>
  dplyr::distinct()
mous$specie = "mouse"
macs$specie = "macaque"
mars$specie = "marmoset"
mous@meta.data = left_join(mous@meta.data |> rownames_to_column("cellname"), clu3p5)  |> 
  column_to_rownames("cellname")
macs@meta.data = left_join(macs@meta.data |> rownames_to_column("cellname"), clu3p5)  |> 
  column_to_rownames("cellname")
mars@meta.data = left_join(mars@meta.data |> rownames_to_column("cellname"), clu3p5)  |> 
  column_to_rownames("cellname")

write_csv(mous[[]] |> rownames_to_column("cellname"), "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/Banksy-cluster_mouse250722.clu4.csv")
write_csv(macs[[]] |> rownames_to_column("cellname"), "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/Banksy-cluster_macaque250722.clu4.csv")
write_csv(mars[[]] |> rownames_to_column("cellname"), "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/Banksy-cluster_marmoset250722.clu4.csv")


################# read clu5/clu6 formal
mouse_clu6 = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/251202_final_clu6_csv/Banksy-cluster_mouse250722.clu6_251201.csv")
mous$clu6 = mouse_clu6$clu6
macaque_clu6 = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/251202_final_clu6_csv/Banksy-cluster_macaque250722.clu6_251201.csv")
macs$clu6 = macaque_clu6$clu6
marmoset_clu6 = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/251202_final_clu6_csv/Banksy-cluster_marmoset250722.clu6_251202.csv")
mars$clu6 = marmoset_clu6$clu6

saveRDS(mous, "/cluster/share/hypo_group/Spatial/spatial_rds_20250328/mouse.rds")
saveRDS(macs, "/cluster/share/hypo_group/Spatial/spatial_rds_20250328/macaque.rds")
saveRDS(mars, "/cluster/share/hypo_group/Spatial/spatial_rds_20250328/marmoset.rds")

oto_pw(mous, col="clu6", pre_fn="mouse.anno.250722.clu6")
oto_pw(macs, col="clu6", pre_fn="macaque.anno.250722.clu6")
oto_pw(mars, col="clu6", pre_fn="marmoset.anno.250722.clu6")

###################### get consensus marker
###### load mapping of 3 specie
conversion_table = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/conv_tb.gprofiler.20230729.csv")
conversion_table = as_tibble(lapply(conversion_table, 2, FUN = gsub, pattern = '_', replacement = '-'))
macaque_map = setNames(conversion_table$mouseGene, conversion_table$macaqueGene)
marmoset_map = setNames(conversion_table$mouseGene, conversion_table$marmosetGene)

macaque_map_rev = setNames(conversion_table$macaqueGene, conversion_table$mouseGene)
marmoset_map_rev = setNames(conversion_table$marmosetGene, conversion_table$mouseGene)
macaque_map_revr = setNames(conversion_table$mouseGene, conversion_table$macaqueGene)
marmoset_map_revr = setNames(conversion_table$mouseGene, conversion_table$marmosetGene)

mouse_clu4_mks = read_csv(Sys.glob("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/mouse.anno.250722.clu4one_to_other*.csv.gz")) |>
  mutate(gene_mouse = gene, specie="mouse")
macaque_clu4_mks = read_csv(Sys.glob("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/macaque.anno.250722.clu4one_to_other*.csv.gz")) |>
  mutate(gene_mouse = macaque_map[gene], specie="macaque")
marmoset_clu4_mks = read_csv(Sys.glob("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/marmoset.anno.250722.clu4one_to_other*.csv.gz")) |>
  mutate(gene_mouse = marmoset_map[gene], specie="marmoset")


####### not merge result
csmks_not_merge = bind_rows(mouse_clu4_mks, macaque_clu4_mks, marmoset_clu4_mks) |> # consensus markers
  group_by(cluster) |>
  mutate(cluster_nsp = length(unique(specie))) |>
  ungroup() |>
  filter(!is.na(gene_mouse)) |>
  filter(p_val_adj<0.05 & avg_log2FC>1) |>
  arrange(cluster,gene_mouse) |>
  group_by(cluster, gene_mouse, cluster_nsp) |>
  summarize(n=n(), p_val_adj=mean(p_val_adj), avg_log2FC=mean(avg_log2FC)) |>
  ungroup()

csmks = bind_rows(mouse_clu4_mks, macaque_clu4_mks, marmoset_clu4_mks) |> # consensus markers
  group_by(cluster) |>
  mutate(cluster_nsp = length(unique(specie))) |>
  ungroup() |>
  filter(!is.na(gene_mouse)) |>
  filter(p_val_adj<0.05 & avg_log2FC>1) |>
  arrange(cluster,gene_mouse) |>
  group_by(cluster, gene_mouse, cluster_nsp) |>
  summarize(n=n(), p_val_adj=mean(p_val_adj), avg_log2FC=mean(avg_log2FC)) |>
  ungroup() |>
  filter(n == cluster_nsp) 
write_csv(csmks, "parcellation/Banksy/consensus.mks.csv")



csmks = csmks_not_merge
csmks_top = csmks |>
    group_by(cluster) |>
    slice_max(avg_log2FC, n=10) |>
    ungroup() |>
    mutate(BA_order=cluster, Gene_order=gene_mouse)
write_csv(csmks_top, "parcellation/Banksy/consensus.mks.top.csv")




