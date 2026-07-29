#! /usr/bin/env R

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
  library(svglite)
  library(tidyverse)
  library(reticulate)

  # specific lib
  library(Seurat)
  library(ComplexHeatmap)
})
###########################################
setwd("/public/group/weiwu_group/majianyun/data/Hypothalamus/raw.all.snRNA")
options(future.globals.maxSize=+Inf)
# theme set
theme_set(theme_light())
# plan(multicore, workers=16)



seu = readRDS("/cluster/share/hypo_group/snRNA_clean/all_cell/all_cell_sct_1122.rds")
seu = subset(seu, idents=c("GABA", "GLU"), invert=TRUE)
### anno4 is to be refined

anno_name = as.factor(seu$anno4)
levels(anno_name)
reps = c("a1"="TC_a1", "a2"="TC_a2", "b1"="TC_b1", "b2"="TC_b2",
         "as0" = "AS0", "as1"="AS1", "as2"= "AS2",
         "Cycling_opc"="OPC_cycling", 
         "epc0"="EPC0", "epc1"="EPC1",
         "mg0"="MG0", "mg1"="MG1", "mg2"="MG2", "mg3"="MG3", "mg4"="MG4", "mg5"="MG5",
         "NFO"="NFO0", "NFO2"="NFO2",
         "OD1"="OD0", "OD2"="OD1", "OD3"="OD2")

levels(anno_name) = ifelse(levels(anno_name) %in% names(reps), reps, levels(anno_name))
table(anno_name)
seu$annof = anno_name
# write_tsv(seu[[]] |> rownames_to_column("cellname"), "all_cell_sct_1122.tsv")

seuf = seu |> 
          RunUMAP(dims=1:50, reduction="integrated.rpca.full") |> 
          FindNeighbors(reduction="integrated.rpca.full") |>
          FindClusters(resolution=0.5) 

map_mat = as.matrix(table(seuf$anno3, seuf$seurat_clusters))
map_vec = setNames(rownames(map_mat)[apply(map_mat,2, which.max)], colnames(map_mat))
seuf$anno_new = map_vec[as.character(seuf$seurat_clusters)]
Idents(seuf) = "anno_new"

### plot
DimPlot(seuf, label=T) + NoLegend()
DimPlot(seuf, group.by="anno3", label=T) + NoLegend()
DimPlot(seuf, group.by="anno_new", label=T) + NoLegend()

### run SCT first?
pdims=30


### AS
as = subset(seuf, idents="AS")
asl = as |> SplitObject("specie") |> lapply(function(s) {SCTransform(s, method="glmGamPoi")})
inte_f = SelectIntegrationFeatures(asl, nfeatures=3000)
asri = PrepSCTIntegration(object.list=asl, anchor.features=inte_f) %>%
  FindIntegrationAnchors(normalization.method="SCT", anchor.features=inte_f) %>%
  IntegrateData(normalization.method="SCT") %>%
  RunPCA(npcs=pdims, verbose=F) %>%
  FindNeighbors(dims=1:pdims, reduction="pca", nn.eps=0.0) %>%
  FindClusters(resolution=0.1) %>%
  RunUMAP(dims=1:pdims)



DimPlot(asri, label=T) + NoLegend()
asm = FindAllMarkers(asri, only.pos=T)
asmr = asm |> group_by(cluster) |> filter(p_val_adj < 1e-3 & abs(avg_log2FC) > 2) |> ungroup()
table(asm$cluster)
table(asri$seurat_clusters, asri$specie)
as_df = asri[[]]
as_df %<>%
  mutate(clu=case_when(
    seurat_clusters == "0" ~ "AS1",
    seurat_clusters == "1" ~ "AS1",
    seurat_clusters == "4" ~ "AS1",
    seurat_clusters == "6" ~ "AS1",
    seurat_clusters == "2" ~ "AS2",
    seurat_clusters == "3" ~ "AS3",
    T ~ "problem"
  ))
asri$clu = as_df$clu
asri_wd = asri
asri = asri_wd[, asri_wd$clu != "problem"]
DimPlot(asri, group.by = "clu") 
DimPlot(asri, group.by = "specie") 
asri %<>% RunUMAP(dims = 1:20, seed.use = 42)
DimPlot(asri, label=T, group.by="clu") + NoLegend()

saveRDS(asri, "rds/AS.2024.rds")



### MG
mg = subset(seuf, idents="MG")
load("zh_MG.Rdata")
mgl = mg |> SplitObject("specie") |> lapply(function(s) {SCTransform(s, method="glmGamPoi")})
inte_f = SelectIntegrationFeatures(mgl, nfeatures=3000)
mgri = PrepSCTIntegration(object.list=mgl, anchor.features=inte_f) %>%
  FindIntegrationAnchors(normalization.method="SCT", anchor.features=inte_f) %>%
  IntegrateData(normalization.method="SCT") %>%
  RunPCA(npcs=pdims, verbose=F) %>%
  FindNeighbors(dims=1:pdims, reduction="pca", nn.eps=0.0) %>%
  FindClusters(resolution=0.8) %>%
  FindSubCluster(cluster = '8', graph.name = "integrated_snn", resolution = 2) %>%   ### partition MG5
  RunUMAP(dims=1:pdims)

mgri$clu = setNames(meta$cell, rownames(meta))[colnames(mgri)]
Idents(mgri) = mgri$clu
mgri = mgri[, !is.na(mgri$clu)]

saveRDS(mgri, "rds/MG.2024.rds")

### TC
tc = subset(seuf, idents="TC")
tcl = tc |> SplitObject("specie") |> lapply(function(s) {SCTransform(s, method="glmGamPoi")})
inte_f = SelectIntegrationFeatures(tcl, nfeatures=3000)
tcri = PrepSCTIntegration(object.list=tcl, anchor.features=inte_f) %>%
  FindIntegrationAnchors(normalization.method="SCT", anchor.features=inte_f) %>%
  IntegrateData(normalization.method="SCT") %>%
  RunPCA(npcs=pdims, verbose=F) %>%
  FindNeighbors(dims=1:pdims, reduction="pca", nn.eps=0.0) %>%
  FindClusters(resolution=0.1) %>%
  RunUMAP(dims=1:pdims)


DimPlot(tcri, label=T) + NoLegend()
table(tcri$seurat_clusters, tc$anno4)
tcm = FindAllMarkers(tcri, only.pos=T)
table(tcri$seurat_clusters, tc$specie)
tc_df = tcri[[]]
tc_df %<>%
  mutate(clu=case_when(
    seurat_clusters == "0" ~ "TCα",
    seurat_clusters == "1" ~ "TCβ",
    T ~ "problem"
  ))
tcri$clu = tc_df$clu
tcri = subset(tcri, subset = clu != "problem")
DimPlot(tcri, group.by = "clu")
Idents(tcri) = "clu"
##### a1 a2 b1 b2
FeaturePlot(tcri, c("Slc17a8", "Vcan", "Frzb", 
  "Scn7a", "Rax", "Ccdc153", "Col25a1",
  "Pdzph1", "Lyz2", "Mia",
  "Col23a1", "Slc16a2", "Rax"))

saveRDS(tcri, "rds/TC.2024.rds")
DimPlot(tcri, group.by = "specie")



### epc

epc = subset(seuf, idents="EpC")
epcl = epc |> SplitObject("specie") |> lapply(function(s) {SCTransform(s, method="glmGamPoi")})
inte_f = SelectIntegrationFeatures(epcl, nfeatures=3000)
epcri = PrepSCTIntegration(object.list=epcl, anchor.features=inte_f) %>%
  FindIntegrationAnchors(normalization.method="SCT", anchor.features=inte_f) %>%
  IntegrateData(normalization.method="SCT") %>%
  RunPCA(npcs=pdims, verbose=F) %>%
  FindNeighbors(dims=1:pdims, reduction="pca", nn.eps=0.0) %>%
  FindClusters(resolution=0.1) %>%
  RunUMAP(dims=1:pdims)
DimPlot(epcri, label=T) + NoLegend()
epcm = FindAllMarkers(epcri, only.pos=T)
table(epcri$seurat_clusters, epcri$specie)
epc_df = epcri[[]]
epc_df %<>%
  mutate(clu=case_when(
    seurat_clusters == "0" ~ "EpC1",
    seurat_clusters == "1" ~ "EpC2",
    T ~ "problem"
  ))
epcri$clu = epc_df$clu
Idents(epcri) = epcri$clu
FeaturePlot(epcri, c("B2m", "Dmd"))
saveRDS(epcri, "rds/EpC.2024.rds")


### OD
od = subset(seuf, idents="OD")
odl = od |> SplitObject("specie") |> lapply(function(s) {SCTransform(s, method="glmGamPoi")})
inte_f = SelectIntegrationFeatures(odl, nfeatures=3000)
odri = PrepSCTIntegration(object.list=odl, anchor.features=inte_f) %>%
  FindIntegrationAnchors(normalization.method="SCT", anchor.features=inte_f) %>%
  IntegrateData(normalization.method="SCT") %>%
  RunPCA(npcs=pdims, verbose=F) %>%
  FindNeighbors(dims=1:pdims, reduction="pca", nn.eps=0.0) %>%
  FindClusters(resolution=0.1) %>%
  RunUMAP(dims=1:pdims)
DimPlot(odri, label=T) + NoLegend()
odm = FindAllMarkers(odri, only.pos=T)
table(odri$seurat_clusters, odri$specie)
odri = subset(odri, idents=c("3"), invert=T)
od_df = odri[[]]
od_df %<>%
  mutate(clu=case_when(
    seurat_clusters == "0" ~ "OD1",
    seurat_clusters == "1" ~ "OD2",
    seurat_clusters == "2"~ "OD3",
    T ~ "problem"
  ))
odri$clu = od_df$clu
FeaturePlot(odri, "Spp1")
DimPlot(odri, group.by="specie")
saveRDS(odri, "rds/OD.2024.rds")

### OPC
opc = subset(seuf, idents="OPC")
opcl = opc |> SplitObject("specie") |> lapply(function(s) {SCTransform(s, method="glmGamPoi")})
inte_f = SelectIntegrationFeatures(opcl, nfeatures=3000)
opcri = PrepSCTIntegration(object.list=opcl, anchor.features=inte_f) %>%
  FindIntegrationAnchors(normalization.method="SCT", anchor.features=inte_f) %>%
  IntegrateData(normalization.method="SCT") %>%
  RunPCA(npcs=pdims, verbose=F) %>%
  FindNeighbors(dims=1:pdims, reduction="pca", nn.eps=0.0) %>%
  FindClusters(resolution=0.1) %>%
  RunUMAP(dims=1:pdims)

DimPlot(opcri, label=T) + NoLegend()
opcm = FindAllMarkers(opcri, only.pos=T)
table(opcri$seurat_clusters, opcri$specie)
opc_df = opcri[[]]
opc_df %<>%
  mutate(clu=case_when(
    seurat_clusters == "0" ~ "OPC",
    seurat_clusters == "1" ~ "DOPC",
    seurat_clusters == "2" ~ "OPC",
    seurat_clusters == "4" ~ "COPC",
    T ~ "problem"
  ))
opcri$clu = opc_df$clu

Idents(opcri) = "clu"
opcri = subset(opcri, idents="problem", invert=T)


FeaturePlot(opcri, c("Cst3", "Sox2"))

saveRDS(opcri, "rds/OPC.2024.rds")

### END
end = subset(seuf, idents="END")
endl = end |> SplitObject("specie") |> lapply(function(s) {SCTransform(s, method="glmGamPoi")})
inte_f = SelectIntegrationFeatures(endl, nfeatures=3000)
endri = PrepSCTIntegration(object.list=endl, anchor.features=inte_f) %>%
  FindIntegrationAnchors(normalization.method="SCT", anchor.features=inte_f) %>%
  IntegrateData(normalization.method="SCT") %>%
  RunPCA(npcs=pdims, verbose=F) %>%
  FindNeighbors(dims=1:pdims, reduction="pca", nn.eps=0.0) %>%
  FindClusters(resolution=0.1) %>%
  RunUMAP(dims=1:pdims)
DimPlot(endri, label=T) + NoLegend()
endm = FindAllMarkers(endri, only.pos=T)
table(endri$seurat_clusters, endri$specie)
end_df = endri[[]]
end_df %<>%
  mutate(clu=case_when(
    seurat_clusters == "0" ~ "END1", 
    seurat_clusters == "1" ~ "END2",
    seurat_clusters == "2" ~ "END3",
    seurat_clusters == "3" ~ "END4",
    T ~ "problem"
  ))
endri$clu = end_df$clu
Idents(endri) = "clu"
saveRDS(endri, "rds/END.2024.rds")

### VLMC
vlmc = subset(seuf, idents="VLMC")
vlmcl = vlmc |> SplitObject("specie") |> lapply(function(s) {SCTransform(s, method="glmGamPoi")})
inte_f = SelectIntegrationFeatures(vlmcl, nfeatures=3000)
DimPlot(vlmc)
vlmc$clu = "VLMC"
saveRDS(vlmc, "rds/VLMC.2024.rds")

### PTC
ptc = subset(seuf, idents="PTC")
ptcl = ptc |> SplitObject("specie") |> lapply(function(s) {SCTransform(s, method="glmGamPoi")})
inte_f = SelectIntegrationFeatures(ptcl, nfeatures=3000)
ptcri = PrepSCTIntegration(object.list=ptcl, anchor.features=inte_f) %>%
  FindIntegrationAnchors(normalization.method="SCT", anchor.features=inte_f) %>%
  IntegrateData(normalization.method="SCT") %>%
  RunPCA(npcs=pdims, verbose=F) %>%
  FindNeighbors(dims=1:pdims, reduction="pca", nn.eps=0.0) %>%
  FindClusters(resolution=0.1) %>%
  RunUMAP(dims=1:pdims)
DimPlot(ptcri)
ptcri$clu = "PTC"
saveRDS(ptcri, "rds/PTC.2024.rds")



### merge
seun = merge(asri, list(odri, opcri, vlmc, endri, ptcri, tcri, epcri, mgri))
seuf = seuf[,colnames(seun)]
seuf$clu = seun$clu
DimPlot(seuf, label=T, group.by="clu") + NoLegend()
Idents(seuf) = "clu"
saveRDS(seuf, "rds/non_neu.2024.rds")
seufm = seuf[, seuf$sex == "M"]
saveRDS(seufm, "rds/non_neu.2024.M.rds")



### extra 
OD = readRDS("rds/OD.2024.rds")
OPC = readRDS("rds/OPC.2024.rds")
oligo = merge(OD, OPC)
DefaultAssay(oligo) = "RNA"
oligo = DietSeurat(oligo, assay=c( "RNA"))

oligol = oligo |> SplitObject("specie") |> lapply(function(s) {SCTransform(s, method="glmGamPoi")})
inte_f = SelectIntegrationFeatures(oligol, nfeatures=3000)
oligori = PrepSCTIntegration(object.list=oligol, anchor.features=inte_f) %>%
  FindIntegrationAnchors(normalization.method="SCT", anchor.features=inte_f) %>%
  IntegrateData(normalization.method="SCT") %>%
  RunPCA(npcs=pdims, verbose=F) %>%
  FindNeighbors(dims=1:pdims, reduction="pca", nn.eps=0.0) %>%
  FindClusters(resolution=0.1) %>%
  RunUMAP(dims=1:pdims)
DimPlot(oligori, group.by="clu", label=T, ) + NoLegend()
saveRDS(oligori, "rds/Oligo.2024.rds")

