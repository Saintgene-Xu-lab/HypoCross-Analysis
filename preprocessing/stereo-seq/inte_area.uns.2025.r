#! /usr/bin/env R

# integration of area (without transfer)

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
  library(ComplexHeatmap)
  library(reshape2)
  library(ClusterR)
  library(Ternary)
  library(uwot)
  library(aricode)
  library(ggrepel)
  library(ggridges)
  library(SpacoR)
  library(aricode)
  library(presto)
  library(qs)
  
})
###########################################
setwd("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial")
theme_set(theme_classic())
# plan(multicore, workers=50)
options(future.globals.maxSize=+Inf)

##### lib
source("/public/group/weiwu_group/majianyun/data/Hypothalamus/cross_inte_lib.r")
source("/public/group/weiwu_group/majianyun/data/Hypothalamus/draw_spatial_lib.r")



##### rename function
RenameGenesSeurat <- function(obj, conversion_table, old_name, new_name) { 
  # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  name_set = setNames(conversion_table[[new_name]],
                      conversion_table[[old_name]])
  obj = obj[rownames(obj) %in% conversion_table[[old_name]],]
  
  RNA <- obj@assays$RNA
  if (!is.null(RNA)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- name_set[RNA@counts@Dimnames[[1]]]
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- name_set[RNA@data@Dimnames[[1]]]
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- name_set[RNA@scale.data@Dimnames[[1]]]
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA

  SCT <- obj@assays$SCT
  if ((!is.null(SCT))) {
    if (length(SCT@counts)) SCT@counts@Dimnames[[1]]            <- name_set[SCT@counts@Dimnames[[1]]]
    if (length(SCT@data)) SCT@data@Dimnames[[1]]                <- name_set[SCT@data@Dimnames[[1]]]
    if (length(SCT@scale.data)) dimnames(SCT@scale.data)[[1]]    <- name_set[dimnames(SCT@scale.data)[[1]]]
    SCT@var.features = name_set[SCT@var.features]
  } else {"Unequal gene sets: nrow(SCT) != nrow(newnames)"}
  obj@assays$SCT <- SCT
  return(obj)
}
###########################################



#### read spatial data
#### all half!!!

mous = readRDS("/cluster/share/hypo_group/Spatial/spatial_rds_20250328_with_anno/mouse.half.rds")
macs = readRDS("/cluster/share/hypo_group/Spatial/spatial_rds_20250328_with_anno/macaque.half.rds")
mars = readRDS("/cluster/share/hypo_group/Spatial/spatial_rds_20250328_with_anno/marmoset.half.rds")
uarea = unique(c(mous$mainAreaID, macs$mainAreaID, mars$mainAreaID))
##### rename
conversion_table = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/conv_tb.gprofiler.20230729.csv")
conversion_table = as_tibble(lapply(conversion_table, 2, FUN = gsub, pattern = '_', replacement = '-'))
##### add info
mous$specie = "mouse"
macs$specie = "macaque"
mars$specie = "marmoset"
mseul = list(mous, macs, mars)
names(mseul) = c("mouse", "macaque", "marmoset")


############## for single specie
ebdl = lapply(mseul, function(seu) {
  ##### remove gene with no ortholog expression
  cells = which(colSums(seu@assays$RNA@counts) != 0)
  if (length(cells) != 0) {
    seu = seu[, cells]
  }
  seu %<>%
    SCTransform(method="glmGamPoi") %>%
    RunPCA() %>%
    RunUMAP(dims=1:30) %>%
    FindNeighbors(dims=1:30) %>%
    FindClusters(resolution=1)

    print(table(seu$seurat_clusters))

    seu$cluster_label = Idents(seu)

    tb = seu[[]]
    tb$specie_maid = paste(tb$specie, tb$mainAreaID)
    tbl = tb |> rownames_to_column("cell") |> group_split(chip) 
    ebdl = future_lapply(tbl, get_nearby_embedding, k=100, cluster_labs="seurat_clusters", smooth=F, future.seed=T)
    acnames = Reduce(union, sapply(ebdl, colnames))
    ebdl = future_lapply(ebdl, function(x) {
      if (identical(sort(acnames), sort(colnames(x)))) {
        return(x)
      }
      mn = setdiff(acnames, colnames(x)) # miss name
      y = matrix(0, nrow = nrow(x), ncol = length(mn), 
        dimnames=list(rownames(x), mn)
      )
      x = cbind(x, y)
      x[, sort(colnames(x))]
    })
    ebd = Reduce(rbind, ebdl)
    ebd = ebd[rownames(tb), ]
    ebd_pca = prcomp(ebd)$x[,1:30]
  return(ebd_pca)
})

qsave(ebdl, "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/ebdl.split.2025.qs")


ebdl = qread("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/ebdl.split.2025.qs")
resos = 1:500
plan(multicore, workers=50)
res = lapply(ebdl, function(ebd) {
   # 2: 62 46 35 
  suppressMessages({
    kmclud = future_sapply(resos, function(reso) {
      km_mb = MiniBatchKmeans(ebd, clusters=reso, batch_size=40, num_init=200, max_iters=100, 
                            init_fraction=0.25, initializer='kmeans++', early_stop_iter=10,
                            verbose=F, seed=42)
      print(reso)
      kmclu = predict(object=km_mb, newdata=ebd)
    }, future.seed=T)
  })
  # kmclud
})

qsave(ebdl, "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/ebdl.split.res.2025.qs")


future_lapply(seq_along(res), function(i) {
  kmclul = res[[i]]
  colnames(kmclul) = paste0("clu_", resos)
  kmclul = as_tibble(kmclul)

  seu = mseul[[i]]
  cells = which(colSums(seu@assays$RNA@counts) != 0)
  if (length(cells) != 0) {
    seu = seu[, cells]
  }
  seu$cluster_label = Idents(seu)

  tb = seu[[]]
  kmclul$name = rownames(tb)
  kmclul$rx = tb$rx
  kmclul$ry = tb$ry
  kmclul$chip = tb$chip
  kmclul$mainAreaID = tb$mainAreaID
  kmclul$hid = tb$hid
  kmclul$specie = tb$specie
  write.csv(kmclul, file=gzfile(paste0("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/single_cell.c1_100.", tb$specie[1] , ".csv.gz")))
  seu[[paste0("clu_", resos[length(resos)])]] = as.factor(kmclul[[paste0("clu_", resos[length(resos)])]]) 
  mks = wilcoxauc(seu, paste0("clu_", resos[length(resos)]))
  write_csv(mks, paste0("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/single_cell.c1_100.", tb$specie[1] , ".mks.csv.gz"))
  mks_300 = top_markers(mks, 300)
  write_csv(mks_300, paste0("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/single_cell.c1_100.", tb$specie[1] , ".mks_300.csv.gz"))
})






##### for integrated data
plan(sequential) # prevent out of memory
mseul = list(mous, macs, mars)
names(mseul) = c("mouse", "macaque", "marmoset")

mseul[["mouse"]] = mseul[["mouse"]][rownames(mseul[["mouse"]]) %in% conversion_table$mouseGene,]
mseul[["macaque"]] = RenameGenesSeurat(mseul[["macaque"]], conversion_table, "macaqueGene", "mouseGene")
mseul[["marmoset"]] = RenameGenesSeurat(mseul[["marmoset"]], conversion_table, "marmosetGene", "mouseGene")

mseu_ori = merge(mseul[["mouse"]], c(mseul[["macaque"]], mseul[["marmoset"]]))
qsave(mseu_ori, "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/rds/mseu_ori.3s.qs")

mseul = list(mous, macs, mars)
names(mseul) = c("mouse", "macaque", "marmoset")

mseul[["mouse"]] = mseul[["mouse"]][rownames(mseul[["mouse"]]) %in% conversion_table$mouseGene,]
mseul[["macaque"]] = RenameGenesSeurat(mseul[["macaque"]], conversion_table, "macaqueGene", "mouseGene")
mseul[["marmoset"]] = RenameGenesSeurat(mseul[["marmoset"]], conversion_table, "marmosetGene", "mouseGene")

mseul = future_lapply(mseul, function(seu) {
  ##### remove gene with no ortholog expression?
  cells = which(colSums(seu@assays$RNA@counts) != 0)
  if (length(cells) != 0) {
    seu = seu[, cells]
  }
  seu %<>%
    SCTransform(method="glmGamPoi") %>%
    RunPCA() %>%
    RunUMAP(dims=1:30) %>%
    FindNeighbors(dims=1:30) %>%
    FindClusters(resolution=1)

    seu$cluster_label = Idents(seu)
    seu
})
qsave(mseul, "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/mseul_sct.qs")



options(future.globals.maxSize=1024^3*100)
pdims=30
inte_f = SelectIntegrationFeatures(mseul, nfeatures=1000)
iseu = PrepSCTIntegration(object.list=lapply(mseul, function(s) {s[inte_f,]}), anchor.features=inte_f) %>%
        FindIntegrationAnchors(normalization.method="SCT", anchor.features=inte_f, reduction = "rpca") %>%
        IntegrateData(normalization.method="SCT") %>%
        RunPCA(npcs=pdims, verbose=F) %>%
        FindNeighbors(dims=1:pdims, reduction="pca", nn.eps=0.0) %>%
        FindClusters(resolution=1) %>%
        RunUMAP(dims=1:pdims)

saveRDS(iseu, "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/rds/iseu.20250303.rpca.rds")


########### embedding processing
tb = iseu[[]]
tb$specie_maid = paste(tb$specie, tb$mainAreaID)
tbl = tb |> rownames_to_column("cell") |> group_split(chip) 
set.seed(42)
plan(multicore, workers=25)
ebdl = future_lapply(tbl, get_nearby_embedding, k=100, cluster_labs="seurat_clusters", smooth=F, future.seed=T)
acnames = Reduce(union, sapply(ebdl, colnames))
ebdl = future_lapply(ebdl, function(x) {
  if (identical(sort(acnames), sort(colnames(x)))) {
    return(x)
  }
  mn = setdiff(acnames, colnames(x)) # miss name
  y = matrix(0, nrow = nrow(x), ncol = length(mn), 
    dimnames=list(rownames(x), mn)
  )
  x = cbind(x, y)
  x[, sort(colnames(x))]
})
ebd = Reduce(rbind, ebdl)
ebd = ebd[rownames(tb), ]
ebd_pca = prcomp(ebd)$x[,1:30]
ebd = ebd_pca
qsave(ebd, "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/ebdl.inte.qs")


#### MBK computation
ebd = qread("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/ebdl.inte.qs")
set.seed(42)
clu_step = seq(1,100)#seq(51, 200, 20)
kmclul = future_sapply(clu_step, function(i) {
  km_mb = MiniBatchKmeans(ebd, clusters=i, batch_size=40, num_init=20, max_iters=100, 
                      init_fraction=0.25, initializer='kmeans++', early_stop_iter=10,
                      verbose=F, seed=42)
  predict(object=km_mb, newdata=ebd)
}, future.seed=42)
colnames(kmclul) = paste0("clu_", clu_step)
kmclul = as_tibble(kmclul)
kmclul$name = rownames(tb)
kmclul$rx = tb$rx
kmclul$ry = tb$ry
kmclul$chip = tb$chip
kmclul$mainAreaID = tb$mainAreaID
kmclul$hid = tb$hid
kmclul$specie = tb$specie
write.csv(kmclul, file=gzfile("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/single_cell.c1_100.HY.csv.gz"))

kmclul = read.csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/single_cell.c1_100.HY.csv.gz", row.names=1)
rownames(kmclul) = colnames(iseu)
iseuk = AddMetaData(iseu, kmclul)
iseukl = iseuk |> SplitObject("specie")
qsave(iseuk, "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/iseu.no_split.qs")
qsave(iseukl, "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/iseu.split.qs")