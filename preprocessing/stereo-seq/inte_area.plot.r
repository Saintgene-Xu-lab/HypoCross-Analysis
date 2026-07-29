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
  library(ggsci)
  library(future.apply)
  library(svglite)
  library(tidyverse)
  library(qs)
  
  # specific lib
  library(ComplexHeatmap)
  library(reshape2)
  library(ClusterR)
  # library(Ternary)
  library(uwot)
  library(aricode)
  library(ggrepel)
  library(ggridges)
  library(SpacoR)
  # library(aricode) # not compatible
  library(presto)
  library(fastcluster)
  library(circlize)
  library(readxl)
})
###########################################
setwd("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial")
theme_set(theme_classic())
plan(multicore, workers=3)
options(future.globals.maxSize=200000000000)


### source for rsf
source("/public/group/weiwu_group/majianyun/data/Hypothalamus/cross_inte_lib.r")
source("/public/group/weiwu_group/majianyun/data/Hypothalamus/draw_spatial_lib.r")


#### read spatial data
mous = readRDS("/gpfs/group/data_sci_group/weiwu_group/group_data/hypo_group/Spatial/spatial_rds_20250328_with_anno/mouse.half.rds")
macs = readRDS("/gpfs/group/data_sci_group/weiwu_group/group_data/hypo_group/Spatial/spatial_rds_20250328_with_anno/macaque.half.rds")
mars = readRDS("/gpfs/group/data_sci_group/weiwu_group/group_data/hypo_group/Spatial/spatial_rds_20250328_with_anno/marmoset.half.rds")

# load additional annotation
mous$clu5 = ifelse(mous$clu5 == "HY1", "HY", mous$clu5)
macs$clu5 = ifelse(macs$clu5 == "HY1", "HY", macs$clu5)
mars$clu5 = ifelse(mars$clu5 == "HY1", "HY", mars$clu5)

mous$clu6 = ifelse(mous$clu6 == "HY1", "HY", mous$clu6)
macs$clu6 = ifelse(macs$clu6 == "HY1", "HY", macs$clu6)
mars$clu6 = ifelse(mars$clu6 == "HY1", "HY", mars$clu6)

mouse_map =  read_csv("/public/group/weiwu_group/majianyun/share_ori/HD_data/Spatial_Hypothalamus/renewed_ST_data_250227/SpeciesBrainRegion-XSJ_location_v5.mouse.csv") |>
  mutate(specie="mouse") |>
  dplyr::distinct(Acronym, Assign_Simple)
macaque_map = read_csv("/public/group/weiwu_group/majianyun/share_ori/HD_data/Spatial_Hypothalamus/renewed_ST_data_250227/SpeciesBrainRegion-XSJ_location_v5.macaque.csv") |>
  mutate(specie="macaque") |>
  dplyr::distinct(Acronym, Assign_Simple)
marmoset_map =  read_csv("/public/group/weiwu_group/majianyun/share_ori/HD_data/Spatial_Hypothalamus/renewed_ST_data_250227/SpeciesBrainRegion-XSJ_location_v5.marmoset.csv") |>
  mutate(specie="marmoset") |>
  dplyr::distinct(Acronym, Assign_Simple)

all_map = unique(bind_rows(mouse_map, macaque_map, marmoset_map)$Assign_Simple)

combined_map = bind_rows(mouse_map, macaque_map, marmoset_map) #|>

all_map = unique(c(mous$clu6, macs$clu6, mars$clu6))



uarea = unique(combined_map$Assign_Simple)

##### rename
mous$mainAreaID_old = mous$mainAreaID
macs$mainAreaID_old = macs$mainAreaID
mars$mainAreaID_old = mars$mainAreaID

mous$hid = left_join(mous[[]], mouse_map, by = join_by(mainAreaID == Acronym))$Assign_Simple
macs$hid = left_join(macs[[]], macaque_map, by = join_by(mainAreaID == Acronym))$Assign_Simple
mars$hid = left_join(mars[[]], marmoset_map, by = join_by(mainAreaID == Acronym))$Assign_Simple

conversion_table = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/conv_tb.gprofiler.20230729.csv")
conversion_table = as_tibble(lapply(conversion_table, 2, FUN = gsub, pattern = '_', replacement = '-'))

mous$specie = "mouse"
macs$specie = "macaque"
mars$specie = "marmoset"

mseul = list(mous, macs, mars)
names(mseul) = c("mouse", "macaque", "marmoset")
mseul = lapply(mseul, NormalizeData)







mseul_sct = future_lapply(mseul, function(seu) {
  ##### remove gene with no ortholog gene expression?
  cells = which(colSums(seu@assays$RNA@counts) != 0)
  if (length(cells) != 0) {
    seu = seu[, cells]
  }
  seu %>%
    SCTransform(method="glmGamPoi") %>%
    RunPCA() %>%
    RunUMAP(dims=1:30) %>%
    FindNeighbors(dims=1:30) %>%
    FindClusters(resolution=1)
})
qsave(mseul_sct, "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/mseul.sct.qs")


mseul_sct = qread("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/mseul.sct.qs")

################ read integrated data
iseu = readRDS("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/rds/iseu.20250303.rpca.rds")
tb = iseu[[]]
tb$specie_maid = paste(tb$specie, tb$mainAreaID)
################ read cluster label

kmclul = read.csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/single_cell.c1_100.HY.csv.gz", row.names=1)
rownames(kmclul) = colnames(iseu)



################################### banksy result load
iseuk = qread("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/iseu.no_split.qs")
iseukl = qread("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/iseu.split.qs")

iseuk$clu5 = ifelse(iseuk$clu5 == "HY1", "HY", iseuk$clu5)
iseukl[[1]]$clu5 = ifelse(iseukl[[1]]$clu5 == "HY1", "HY", iseukl[[1]]$clu5)
iseukl[[2]]$clu5 = ifelse(iseukl[[2]]$clu5 == "HY1", "HY", iseukl[[2]]$clu5)
iseukl[[3]]$clu5 = ifelse(iseukl[[3]]$clu5 == "HY1", "HY", iseukl[[3]]$clu5)

iseuk$clu6 = ifelse(iseuk$clu6 == "HY1", "HY", iseuk$clu6)
iseukl[[1]]$clu6 = ifelse(iseukl[[1]]$clu6 == "HY1", "HY", iseukl[[1]]$clu6)
iseukl[[2]]$clu6 = ifelse(iseukl[[2]]$clu6 == "HY1", "HY", iseukl[[2]]$clu6)
iseukl[[3]]$clu6 = ifelse(iseukl[[3]]$clu6 == "HY1", "HY", iseukl[[3]]$clu6)

res = list()
res$mouse = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/mouse.k30.csv.gz") |> column_to_rownames("rownames")
res$mouse["specie"] = "mouse"
res$mouse = res$mouse[colnames(mous),]
res$mouse$mainAreaID = mous[["mainAreaID"]][rownames(res$mouse), ]
res$mouse$hid = mous[["hid"]][rownames(res$mouse), ]
res$macaque = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/macaque.k60.csv.gz") |> column_to_rownames("rownames")
res$macaque$specie = "macaque"
res$macaque = res$macaque[colnames(macs),]
res$macaque$mainAreaID = macs[["mainAreaID"]][rownames(res$macaque), ]
res$macaque$hid = macs[["hid"]][rownames(res$macaque), ]
res$marmoset = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/marmoset.k60.csv.gz") |> column_to_rownames("rownames")
res$marmoset$specie = "marmoset"
res$marmoset = res$marmoset[colnames(mars),]
res$marmoset$mainAreaID = mars[["mainAreaID"]][rownames(res$marmoset), ]
res$marmoset$hid = mars[["hid"]][rownames(res$marmoset), ]

color_map = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/banksy_colorv2.csv")
color_map = setNames(color_map$color, color_map$nucle)
specie_color = c("mouse"="#00b0ebff", "macaque"="#e2424aff", "marmoset"="#f57f17ff")
#### define used chip
chips = c("T789", "T016_T964", "T318", 
          "T309", "T786", "T010_T968",
          "T008_T970", "T780", "T305")
######################### load over




#################### PCA correlation
iseuk_rna_pca = iseuk
DefaultAssay(iseuk_rna_pca) = "RNA"
iseuk_rna_pca = iseuk_rna_pca |>
  NormalizeData() |>
  FindVariableFeatures(nfeatures=3000) |>
  ScaleData() |>
  RunPCA(npcs=30)
ipca = iseuk_rna_pca[["pca"]]@cell.embeddings
ipca = CreateSeuratObject(t(ipca), meta.data=iseuk[[]])
ipca$clu6_specie = str_c(ipca$clu6, "_", ipca$specie)
ipca = AverageExpression(ipca, group.by="clu6_specie", slot="count")$RNA
ipcor = cor(ipca)
ipcor = ipcor[rownames(ipcor) != "MarMG_marmoset",rownames(ipcor) != "MarMG_marmoset"]
mean_v = c()
mean_mou_mac = c()
mean_mou_mar = c()
mean_mac_mar = c()
for (i in 1:(ncol(ipcor)/3)) {
  i = (i - 1)*3 + 1:3
  ipcor[i, setdiff(1:nrow(ipcor),i)] = NA
  ipcor[setdiff(1:nrow(ipcor),i), i] = NA
  scor = ipcor[i,i]
  mean_v = c(mean_v,  mean(c(scor[2,1], scor[3,1], scor[3,2])))
  mean_mou_mac = c(mean_mou_mac, scor[3,1])
  mean_mou_mar = c(mean_mou_mar, scor[3,2])
  mean_mac_mar = c(mean_mac_mar, scor[2,1])
}
phm = Heatmap(ipcor, cluster_columns=F, cluster_rows=F)
draw(phm)
dev.off()
resp_cor = tibble(correlation=c(mean_mou_mac, mean_mou_mar, mean_mac_mar),
                  correlation_type = rep(c("mouse-macaque correlation", 
                                            "mouse-marmoset correlation", 
                                            "macaque-marmoset correlation"), 
                                                each=length(mean_mou_mac)),
                  clu6=rep(unique(str_split_i(rownames(ipcor), '_', 1)), 3)) |>
              arrange(clu6, correlation_type)
mean_cor_tb = tibble(mean_correlation=mean_v, clu6=unique(str_split_i(rownames(ipcor), '_', 1))) |>
  arrange(desc(mean_correlation)) |>
  mutate(clu6=factor(clu6, levels=clu6))
resp_cor$clu6 = factor(resp_cor$clu6, levels=levels(mean_cor_tb$clu6))
ggplot() +
    geom_col(data=mean_cor_tb, aes(clu6, mean_correlation, fill=clu6)) +
    geom_point(aes(clu6, correlation, shape=correlation_type), data=resp_cor) +
    scale_fill_manual(values=color_map) +R
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/etp.bar.rank.svg", height=5, width=7, fix_text_size=F)

############### marker finding
iseuk$clu_pa = paste(iseuk$clu6, iseuk$specie)
iseuk_spl = mseul |> lapply(NormalizeData)
for (i in seq_along(iseuk_spl)) {
  iseuk_spl[[i]]$clu_pa = paste(iseuk_spl[[i]]$clu6, iseuk_spl[[i]]$specie)
}


################## define marker set
conversion_table = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/conv_tb.gprofiler.20230729.csv")
conversion_table = as_tibble(lapply(conversion_table, 2, FUN = gsub, pattern = '_', replacement = '-'))
macaque_map = setNames(conversion_table$mouseGene, conversion_table$macaqueGene)
marmoset_map = setNames(conversion_table$mouseGene, conversion_table$marmosetGene)
macaque_map_rev = setNames(conversion_table$macaqueGene, conversion_table$mouseGene)
marmoset_map_rev = setNames(conversion_table$marmosetGene, conversion_table$mouseGene)


mouse_clu6_mks = read_csv(Sys.glob("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/mouse.anno.250722.clu6one_to_other*.csv.gz")) |>
  mutate(feature_ortholog = gene, specie="mouse", cluster=ifelse(cluster == "HY1", "HY", cluster))
macaque_clu6_mks = read_csv(Sys.glob("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/macaque.anno.250722.clu6one_to_other*.csv.gz")) |>
  mutate(feature_ortholog = macaque_map[gene], specie="macaque", cluster=ifelse(cluster == "HY1", "HY", cluster))
marmoset_clu6_mks = read_csv(Sys.glob("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/Banksy/marmoset.anno.250722.clu6one_to_other*.csv.gz")) |>
  mutate(feature_ortholog = marmoset_map[gene], specie="marmoset", cluster=ifelse(cluster == "HY1", "HY", cluster))
del_ori = list(mouse_clu6_mks, macaque_clu6_mks, marmoset_clu6_mks) 
qsave(del_ori, "rds/test_all.log.qs")






################ load all data
del_ori = qread("rds/test_all.log.qs")
############ load data and normalize it
named_area = qread("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/rds/mseu_ori.3s.qs")
named_area$clu5 = ifelse(named_area$clu5 == "HY1", "HY", named_area$clu5)
named_area$clu6 = ifelse(named_area$clu6 == "HY1", "HY", named_area$clu6)
named_area$specie_area = str_c(named_area$clu6, ' - ', named_area$specie)
named_area = NormalizeData(named_area)


######################### NHP specific automatic cutofff
del = bind_rows(del_ori) |> # consensus markers
  group_by(cluster) |>
  mutate(cluster_nsp = length(unique(specie))) |>
  ungroup() |>
  filter(!is.na(feature_ortholog) & feature_ortholog %in% conversion_table$mouseGene) |>
  filter(p_val_adj<0.05 & avg_log2FC > 2 & pct.1 > 0.1)
del |> filter(feature_ortholog == "Wif1") 
del |> filter(feature_ortholog == "Cd74") 
NHP_feature_tb = del |> 
  group_by(feature_ortholog) |>
  filter(!any(specie == "mouse") & cluster != "MarMG") |>
  ungroup() |>
  arrange(cluster, feature_ortholog) |>
  mutate(specie_area = str_c(cluster, ' - ', specie),
          ortho_gene_specie_area = str_c(feature_ortholog, specie_area))
write_csv(bind_rows(del_ori)  |>
  filter(p_val_adj<0.05 & avg_log2FC > 2 & pct.1 > 0.1) |> 
  filter(feature_ortholog %in% NHP_feature_tb$feature_ortholog), "NHP_mks.20260429.csv")
NHP_feature_tb_marmg_sup = del |> 
  filter(cluster == "MarMG" & feature_ortholog %in% NHP_feature_tb$feature_ortholog) |>
  arrange(cluster, feature_ortholog) |>
  mutate(specie_area = str_c(cluster, ' - ', specie),
          ortho_gene_specie_area = str_c(feature_ortholog, specie_area))

NHP_feature = unique(NHP_feature_tb$feature_ortholog)
##### NHP
avg_exp_na = AverageExpression(named_area, assays="RNA", slot="data", features=NHP_feature, 
                    group.by="specie_area")
avg_exp_na = avg_exp_na$RNA
NHP_feature_tb_ori = rbind(NHP_feature_tb, NHP_feature_tb_marmg_sup)
NHP_feature_tb = NHP_feature_tb |> distinct(gene, cluster, .keep_all=T)
idx_mat = rbind(matrix(NHP_feature_tb$feature_ortholog, nrow=1),
                matrix(str_c(rep(NHP_feature_tb$cluster, each=3), " - ", 
                  rep(c("mouse", "marmoset", "macaque"), length(NHP_feature_tb$feature_ortholog))), nrow=3))
avg_exp_na2 = apply(idx_mat, 2, function(v) {
  print(avg_exp_na[v[1], v[2:4]])
  avg_exp_na[v[1], v[2:4]]
})
colnames(avg_exp_na2) = NHP_feature_tb$feature_ortholog
rownames(avg_exp_na2) = c("mouse", "marmoset", "macaque")
avg_exp_na = scale(t(avg_exp_na))
avg_exp_na2 = scale(avg_exp_na2)




################ NHP version 1 (plot all)
avg_exp_na_macaque = avg_exp_na
colnames(avg_exp_na_macaque) = macaque_map_rev[colnames(avg_exp_na)]
NHP_feature_sig_region = str_c(NHP_feature_tb_ori$feature_ortholog, NHP_feature_tb_ori$specie_area)
svglite("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/NHP.area.marker.svg", 
        width=ncol(avg_exp_na)/5+2, height=nrow(avg_exp_na)/5+2, fix_text_size=F)
p = Heatmap(avg_exp_na_macaque, cluster_rows=F, 
              heatmap_legend_param = list(
                title="scaled expression",
                direction = "horizontal"),
                cell_fun= function(j, i, x, y, width, height, fill) {
                  if(str_c(colnames(avg_exp_na)[j], rownames(avg_exp_na)[i]) %in% NHP_feature_sig_region) { #
                    grid.text(
                      sprintf("*")#, (t(avg_exp_na))[i, j]
                      , x, y, gp = gpar(fontsize = 15))
                }}
)
draw(p, heatmap_legend_side = "bottom")
dev.off()




###################################################
####### optional step for manual selection ########
###################################################
NHP_filter_tb = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/nhp-specific251204.csv") 
NHP_filter_tb = NHP_filter_tb |>
  uncount(3) |>
  mutate(specie=rep(c("macaque", "marmoset", "mouse"), length(unique(NHP_filter_tb$feature_ortholog)))) |>
  mutate(specie_area = str_c(cluster, ' - ', specie),
          ortho_gene_specie_area = str_c(feature_ortholog, specie_area))
NHP_feature_tb = NHP_feature_tb |>
  mutate(specie_area = str_c(cluster, ' - ', specie),
          ortho_gene_specie_area = str_c(feature_ortholog, specie_area)) |>
  filter(ortho_gene_specie_area %in% NHP_filter_tb$ortho_gene_specie_area) |>
  arrange(cluster, feature_ortholog) 


NHP_feature = unique(NHP_feature_tb$feature_ortholog)
##### NHP
avg_exp_na = AverageExpression(named_area, assays="RNA", slot="data", features=NHP_feature, 
                    group.by="specie_area")
avg_exp_na = avg_exp_na$RNA
NHP_feature_tb_ori = NHP_feature_tb
NHP_feature_tb = NHP_feature_tb |> distinct(gene, cluster, .keep_all=T)
idx_mat = rbind(matrix(NHP_feature_tb$feature_ortholog, nrow=1),
                matrix(str_c(rep(NHP_feature_tb$cluster, each=3), " - ", 
                  rep(c("mouse", "marmoset", "macaque"), length(NHP_feature_tb$feature_ortholog))), nrow=3))
avg_exp_na2 = apply(idx_mat, 2, function(v) {
  print(avg_exp_na[v[1], v[2:4]])
  avg_exp_na[v[1], v[2:4]]
})
colnames(avg_exp_na2) = NHP_feature_tb$feature_ortholog
rownames(avg_exp_na2) = c("mouse", "marmoset", "macaque")
avg_exp_na = scale(t(avg_exp_na))
avg_exp_na2 = scale(avg_exp_na2)


################ NHP version 2 (plot avg_exp_na2)
crange = c(-1, 0, 1)
crg = crange[c(1,3)]
col_fun = colorRamp2(crange, c("skyblue", "white", "red"))
# NHP_feature_tb$cluster = factor(NHP_feature_tb$cluster, levels=unique(NHP_feature_tb$cluster))
nhp_ca = columnAnnotation(
  BKS_region = anno_block(
    gp = gpar(fill = setNames(unname(color_map[unique(NHP_feature_tb$cluster)]), 1:length(unique(NHP_feature_tb$cluster)))),
    labels=factor(unique(NHP_feature_tb$cluster), levels=unique(NHP_feature_tb$cluster)),
    labels_gp = gpar(fontsize = 7)
  )
)

avg_exp_na2_macaque = avg_exp_na2
colnames(avg_exp_na2_macaque) = macaque_map_rev[unname(colnames(avg_exp_na2))]
colnames(avg_exp_na2_macaque) = ifelse(colnames(avg_exp_na2_macaque) == "ENSMFAG00000003837", "ECELL*", colnames(avg_exp_na2_macaque))
NHP_feature_sig_region = str_c(NHP_feature_tb_ori$feature_ortholog, NHP_feature_tb_ori$specie_area)
svglite("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/NHP.area.marker2.svg", 
        width=ncol(avg_exp_na2_macaque)/4+4, height=nrow(avg_exp_na2_macaque)/4+2, fix_text_size=F)
p = Heatmap(avg_exp_na2_macaque, cluster_rows=F, cluster_columns=F, 
              top_annotation=nhp_ca, col=col_fun,
              column_split=factor(NHP_feature_tb$cluster, levels=unique(NHP_feature_tb$cluster)),
              heatmap_legend_param = list(
                title="scaled expression",
                legend_direction = "horizontal", 
                title_position="topcenter",
                at = crg,
                labels = c("Low", "High"),
                legend_width = unit(5, "cm"),
                direction = "horizontal"),
                cell_fun= function(j, i, x, y, width, height, fill) {
                  if(str_c(colnames(avg_exp_na2)[j], NHP_feature_tb$cluster[j],' - ', rownames(avg_exp_na2)[i]) %in% NHP_feature_sig_region) {
                    grid.text(
                      sprintf("*"), x, y, gp = gpar(fontsize = 15))
                }},
                column_title = NULL
)
draw(p, heatmap_legend_side = "bottom")
dev.off()






####### Wif1
library(ggnewscale)
library(ggspatial)
mseul[["mouse"]]$WIF1 = mseul[["mouse"]]@assays$RNA@data["Wif1",]
mseul[["macaque"]]$WIF1 = mseul[["macaque"]]@assays$RNA@data["WIF1",]
mseul[["marmoset"]]$WIF1 = mseul[["marmoset"]]@assays$RNA@data["WIF1",]

scc = c("#a7a7a7", "#1da7b6", "#17b6b1", "#01dccd", "#01ffee")
limits=range(c(mseul[["mouse"]]$WIF1, mseul[["macaque"]]$WIF1, mseul[["marmoset"]]$WIF1))
p1 = ggplot(mseul[["mouse"]][[]] |> filter(chip == chips[3]) |> arrange(WIF1), aes(x=rx, y=ry, col=WIF1)) +
        geom_point_rast(size=1) +
        coord_fixed() +
        scale_y_reverse() +
        scale_color_gradientn(colors=scc, limits=limits) +
        ggtitle("mouse") +
        facet_zoom(x = clu6 %in% c("ARH", "ME"),
                y = clu6  %in% c("ARH", "ME"),
                horizontal = F, zoom.size=1) +
        guides(color="none")  +
        annotation_scale(line_col="black", text_col="black", style="ticks") +
        theme(
          panel.border = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank()
        ) 
p2 = ggplot(mseul[["macaque"]][[]] |> filter(chip == chips[1]) |> arrange(WIF1), aes(x=rx, y=ry, col=WIF1)) +
        geom_point_rast(size=1) +
        coord_fixed() +
        scale_y_reverse() +
        scale_color_gradientn(colors=scc, limits=limits) +
        ggtitle("macaque") +
        facet_zoom(x = clu6 %in% c("ARH", "ME"),
                y = clu6  %in% c("ARH", "ME"),
                horizontal = F, zoom.size=1) +
        guides(color="none")  +
        annotation_scale(line_col="black", text_col="black", style="ticks") +
        theme(
          panel.border = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank()
        ) 
p3 = ggplot(mseul[["marmoset"]][[]] |> filter(chip == chips[2]) |> arrange(WIF1), aes(x=rx, y=ry, col=WIF1)) +
        geom_point_rast(size=1) +
        coord_fixed() +
        scale_y_reverse() +
        scale_color_gradientn(colors=scc, limits=limits) +
        ggtitle("marmoset") +
        facet_zoom(x = clu6 %in% c("ARH", "ME"),
                y = clu6  %in% c("ARH", "ME"),
                shrink=T,
                horizontal = F, zoom.size=1) +
        annotation_scale(line_col="black", text_col="black", style="ticks") +
        theme(
          panel.border = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank()
        ) 
p1 + p2 + p3
ggsave("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/Wif1.exp.png", 
        width=7*2, height=7*1.5, limitsize=F)
ggsave("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/Wif1.exp.svg", 
        width=7*2, height=7*1.5, fix_text_size=F)


############################ draw concensus marker

library(circlize)
crange = c(-3, 0, 3)
crg = crange[c(1,3)]
col_fun = colorRamp2(crange, c("skyblue", "white", "red"))

del = bind_rows(del_ori) |>
  filter(!is.na(feature_ortholog) & feature_ortholog %in% conversion_table$mouseGene) |>
  filter(p_val_adj<0.05 & avg_log2FC > 0.5)
conserved_feature_3 = del |> 
  group_by(feature_ortholog, cluster) |>
  filter(n() == 3) |>
  mutate(consensus_mean=mean(avg_log2FC), order=order(consensus_mean)) |>
  ungroup() |>
  arrange(cluster)
conserved_feature  = conserved_feature_3 |>
  filter(specie == "mouse") |>
  mutate(
    gene_mouse=gene, 
    gene_macaque=macaque_map_rev[gene],
    gene_marmoset=marmoset_map_rev[gene]) |>
  filter(cluster != "MarMG") |>
  arrange(cluster, gene) |>
  dplyr::distinct(gene, .keep_all=T) 

write_csv(conserved_feature_3, "conserved_feature.csv")
write_csv(bind_rows(del_ori) |> left_join(conserved_feature_3[, c("cluster", "feature_ortholog")] |> distinct() |> mutate(conserved=T)) |> filter(!is.na(conserved)) |> arrange(feature_ortholog, cluster, specie), "conserved_feature_mks.20260429.csv")

conserved_feature_sub = read_xlsx("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/conserved_feature.260123.xlsx") |>
  dplyr::rename(feature_macaque=feature_ortholog) |>
  mutate(feature_macaque=toupper(feature_macaque)) |>
  filter(cluster != "POA")

conserved_feature = conserved_feature[conserved_feature$gene_macaque %in% conserved_feature_sub$feature_macaque,]




br_order = c(as.character(unique(conserved_feature$cluster)), setdiff(unique(iseuk$clu6), c("MarMG", as.character(unique(conserved_feature$cluster)))))
rowsp = factor(conserved_feature$cluster, levels=unique(conserved_feature$cluster))
clu6_mean_exp = mseul |> sapply(AverageExpression, group.by="clu6")
mouse_mat = clu6_mean_exp[["mouse.RNA"]][conserved_feature$gene_mouse, br_order]
macaque_mat = clu6_mean_exp[["macaque.RNA"]][conserved_feature$gene_macaque, br_order]
marmoset_mat = clu6_mean_exp[["marmoset.RNA"]][conserved_feature$gene_marmoset, br_order]


ha = rowAnnotation(
  BKS_region = anno_block(
    gp = gpar(fill = setNames(unname(color_map[unique(conserved_feature$cluster)]), 1:length(unique(conserved_feature$cluster)))),
    labels=factor(unique(conserved_feature$cluster), levels=unique(conserved_feature$cluster)),
    width = unit(20, "mm")
  )
)

svglite("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/conserved_markers.mat.svg", 
          width=10, height=20*nrow(conserved_feature)/100, fix_text_size=F)
pmou = Heatmap(t(scale(t(mouse_mat))), col=col_fun, 
            cluster_rows=F, cluster_columns=F,
            row_split=rowsp, row_title=NULL, 
            name = "Z-scored expression per specie",
            heatmap_legend_param = list(
              title_gp=gpar(fontsize=17.5),
              legend_direction = "horizontal", 
              title_position="topcenter",
              at = crg,
              labels = c("Low", "High"),
              legend_width = unit(20, "cm")),
            column_title="mouse")
pmac = Heatmap(t(scale(t(macaque_mat))), col=col_fun, 
            cluster_rows=F, cluster_columns=F,
            row_split=rowsp, row_title=NULL, 
            name = "Z-scored expression per specie",
            heatmap_legend_param = list(
              title_gp=gpar(fontsize=17.5),
              legend_direction = "horizontal", 
              title_position="topcenter",
              at = crg,
              labels = c("Low", "High"),
              legend_width = unit(20, "cm")),
            column_title="macaque",
            right_annotation=ha)
pmar = Heatmap(t(scale(t(marmoset_mat))), col=col_fun, 
            cluster_rows=F, cluster_columns=F,
            row_split=rowsp, row_title=NULL, 
            name = "Z-scored expression per specie",
            heatmap_legend_param = list(
              title_gp=gpar(fontsize=17.5),
              legend_direction = "horizontal", 
              title_position="topcenter",
              at = crg,
              labels = c("Low", "High"),
              legend_width = unit(20, "cm")),
            column_title="marmoset"
            )
draw(pmou + pmar + pmac, merge_legend = TRUE, heatmap_legend_side="bottom")
dev.off()



################# mouse specific features
mouse_feature = bind_rows(del_ori) |>
  mutate(feature_ortholog=case_when(
    specie == "mouse" ~ gene,
    specie == "macaque" ~ unname(macaque_map[gene]),
    specie == "marmoset" ~ unname(marmoset_map[gene]),
  )) |>
  filter(!is.na(feature_ortholog)) |>
  filter(p_val_adj<0.05 & avg_log2FC > 2) |> 
  group_by(feature_ortholog) |>
  filter(!any(specie %in% c("macaque", "marmoset"))) |>
  ungroup()



################################################
#########     overlap & enrich    ##############
################################################

library(eulerr)
csmks_not_merge_pre = bind_rows(del_ori) |> # consensus markers
  group_by(cluster) |>
  mutate(cluster_nsp = length(unique(specie))) |>
  ungroup() |>
  filter(!is.na(feature_ortholog)) |>
  filter(p_val_adj<0.05 & avg_log2FC>0.5)
for (br in unique(csmks_not_merge_pre$cluster)) {
  br_mks = csmks_not_merge_pre |> filter(cluster == br)
  sfeature = split(br_mks$feature_ortholog, br_mks$specie)
  mat = matrix(F, nrow=length(unique(br_mks$feature_ortholog)), ncol=length(sfeature), dimnames=list(unique(br_mks$feature_ortholog), names(sfeature)))
  mat[Pair(br_mks$feature_ortholog, br_mks$specie)] = T
  if (length(sfeature) > 1) {
    svglite(str_c("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/venn_region_mks.", br, ".svg"), fix_text_size=F, width=3, height=3)
    p=plot(euler(mat), quantities = TRUE, fills = list(fill = specie_color[colnames(mat)], col=specie_color[colnames(mat)], alpha = 0.5), main=br)
    plot(p)
    dev.off()
  }
}

################# enrich of brain area specie comparison (include specie-specific marker)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationDbi)
csmks_not_merge_pre_all_gene = bind_rows(del_ori) |> # consensus markers
  group_by(cluster) |>
  mutate(cluster_nsp = length(unique(specie))) |>
  ungroup() |>
  filter(p_val_adj<0.05 & avg_log2FC>0.5)

conv_tb = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/scvi/region_comparison/conversion_table.ma.20241016.dedup.csv")
conv_tb = left_join(conv_tb, AnnotationDbi::select(org.Mm.eg.db,
       keys = conv_tb$mouseID, 
       columns = c("ENSEMBL", "SYMBOL"),
       keytype = "ENSEMBL", multiVals="first") |> dplyr::rename(mouseID=ENSEMBL))
mouse_map = setNames(conv_tb$humanGene, conv_tb$SYMBOL)
macaque_map = setNames(conv_tb$humanGene, conv_tb$macaqueGene)
marmoset_map = setNames(conv_tb$humanGene, conv_tb$marmosetGene)




csmks_not_merge_pre_all_gene$human_gene = case_when(
  csmks_not_merge_pre_all_gene$specie == "mouse" ~ mouse_map[csmks_not_merge_pre_all_gene$gene],
  csmks_not_merge_pre_all_gene$specie == "macaque" ~ macaque_map[csmks_not_merge_pre_all_gene$gene],
  csmks_not_merge_pre_all_gene$specie == "marmoset" ~ marmoset_map[csmks_not_merge_pre_all_gene$gene]
)
csmks_not_merge_pre_all_gene = csmks_not_merge_pre_all_gene |> filter(!is.na(human_gene))
future_lapply(split(csmks_not_merge_pre_all_gene, csmks_not_merge_pre_all_gene$cluster), function(tbs) {
  tbs = left_join(tbs, AnnotationDbi::select(org.Hs.eg.db,
       keys = unique(tbs$human_gene), 
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL") |> dplyr::rename(human_gene=SYMBOL))
  tbs = tbs |>
    group_by(human_gene) |>
    summarize(n=n(), specie=case_when(
      n == 3 ~ "Conserved",
      n == 2 & all(c("marmoset", "macaque") %in% specie) ~ "NHP",
      n == 1 & specie == "macaque" ~ "Macaque",
      n == 1 & specie == "marmoset" ~ "Marmoset",
      T ~ "Other"
    ), cluster=cluster[1], ENTREZID=ENTREZID[1]) |>
    ungroup() |>
    arrange(n) |>
    filter(specie != "Other")


  br = tbs$cluster[1]
  print(br)
  res <- compareCluster(ENTREZID~specie, data=tbs, fun='enrichGO', OrgDb='org.Hs.eg.db')
  if(!is.null(res)) {
    svglite(str_c("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/comp_cluster.", br, ".go.svg"), fix_text_size=F, width=9, height=length(res@compareClusterResult$Description)/3+2)
    p=dotplot(res, title=str_c(br, " GO"), showCategory=30)
    plot(p)
    dev.off()
    write_csv(as_tibble(res@compareClusterResult), str_c("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/comp_cluster.", br, ".go.csv"))
  }
  
  res <- compareCluster(ENTREZID~specie, data=tbs, fun='enrichGO', OrgDb='org.Hs.eg.db', ont="BP")
  if(!is.null(res)) {
    svglite(str_c("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/comp_cluster.", br, ".go_BP.svg"), fix_text_size=F, width=9, height=length(res@compareClusterResult$Description)/3+2)
    p=dotplot(res, title=str_c(br, " GO-BP"), showCategory=30)
    plot(p)
    dev.off()
    write_csv(as_tibble(res@compareClusterResult), str_c("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/comp_cluster.", br, ".go_BP.csv"))
  }

  res <- compareCluster(ENTREZID~specie, data=tbs, fun='enrichGO', OrgDb='org.Hs.eg.db', ont="CC")
  if(!is.null(res)) {
    svglite(str_c("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/comp_cluster.", br, ".go_CC.svg"), fix_text_size=F, width=9, height=length(res@compareClusterResult$Description)/3+2)
    p=dotplot(res, title=str_c(br, " GO-CC"), showCategory=30)
    plot(p)
    dev.off()
    write_csv(as_tibble(res@compareClusterResult), str_c("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/comp_cluster.", br, ".go_CC.csv"))
  }


  res <- compareCluster(ENTREZID~specie, data=tbs, fun='enrichGO', OrgDb='org.Hs.eg.db', ont="MF")
  if(!is.null(res)) {
    svglite(str_c("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/comp_cluster.", br, ".go_MF.svg"), fix_text_size=F, width=9, height=length(res@compareClusterResult$Description)/3+2)
    p=dotplot(res, title=str_c(br, " GO-MF"), showCategory=30)
    plot(p)
    dev.off()
    write_csv(as_tibble(res@compareClusterResult), str_c("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/comp_cluster.", br, ".go_MF.csv"))
  }
  res <- compareCluster(ENTREZID~specie, data=tbs, fun='enrichKEGG', organism="hsa")
  if(!is.null(res)) {
    svglite(str_c("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/comp_cluster.", br, ".kegg.svg"), fix_text_size=F, width=9, height=length(res@compareClusterResult$Description)/3+2)
    p=dotplot(res, title=str_c(br, " KEGG"), showCategory=30)
    plot(p)
    dev.off()
    write_csv(as_tibble(res@compareClusterResult), str_c("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/comp_cluster.", br, ".KEGG.csv"))
  }
})
#################### draw filtered
bmgo_filter = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/filtered8_selected_comp_cluster_BM_GO.csv")
bmgo = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/comp_cluster.BM.go_MF.csv") |>
  filter(Description %in% bmgo_filter$Description) |>
  arrange(desc(Cluster), desc(p.adjust))
bmgo$Description = factor(bmgo$Description, levels=unique(bmgo$Description))
ggplot(bmgo, aes(Cluster, Description, color=Count, size=-log10(p.adjust))) +
  geom_point() +
  scale_color_gradientn(colors=plasma(100)) +
  xlab("Specie") +
  ylab("Term") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("BM-GO:MF")
ggsave("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/enrich.BM.go.svg", 
        height=length(unique(bmgo$Description))/3, width=5.75, fix_text_size=F)


pego_filter = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/filter8_selected_comp_cluster.Pe.go_CC.csv")
pego = read_csv("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/comp_cluster.Pe.go_CC.csv") |>
  filter(Description %in% pego_filter$Description) |>
  arrange(desc(Cluster), desc(p.adjust))
pego$Description = factor(pego$Description, levels=unique(pego$Description))
ggplot(pego, aes(Cluster, Description, color=Count, size=-log10(p.adjust))) +
  geom_point() +
  scale_color_gradientn(colors=plasma(100)) +
  xlab("Specie") +
  ylab("Term") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Pe-GO:CC")
ggsave("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/enrich.Pe.go.svg", 
        height=length(unique(pego$Description))/3+0.3, width=6, fix_text_size=F)










###### cluster level (2025)
get_lab = function(mseu, lab_path) {
  mseu$cell = colnames(mseu) # str_extract(colnames(mseu), "^\\d+")
  nt_pred = lapply(Sys.glob(paste0(lab_path, "/*T*.csv")), read.csv)
  nt_pred = Reduce(rbind, nt_pred)
  as.factor(inner_join(mseu[[]], nt_pred, by="cell")$celltype_pred)
}

mous$celltype_pred = get_lab(mous, "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/label_transfer_evaluation/label/spaitalID/mouse_non8x_56xOD_2xTC_0.4xEN1_0.2xIN1_IN13mk_85c_csv")
macs$celltype_pred = get_lab(macs, "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/label_transfer_evaluation/label/spaitalID/macaque_non20x_40xOD_5xTC_0.5xEN1_0.2xIN1_IN13mk_85c_csv")
mars$celltype_pred = get_lab(mars, "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/label_transfer_evaluation/label/spaitalID/marmoset_non8x_10xOD_0.4xEN1_0.2xIN1_IN13mk_85c_csv")

ctm = bind_rows(list(mous[[]], macs[[]], mars[[]]))

ctm = as.data.frame.matrix(table(str_c(ctm$clu6, rep(' ', length(ctm$clu6)), ctm$specie), ctm$celltype_pred))
ctm_ori = ctm
ctm[, str_ends(rownames(ctm), "mouse")] = ctm[, str_ends(rownames(ctm), "mouse")]/rowSums(ctm[, str_ends(rownames(ctm), "mouse")])
ctm[, str_ends(rownames(ctm), "macaque")] = ctm[, str_ends(rownames(ctm), "macaque")]/rowSums(ctm[, str_ends(rownames(ctm), "macaque")])
ctm[, str_ends(rownames(ctm), "marmoset")] = ctm[, str_ends(rownames(ctm), "marmoset")]/rowSums(ctm[, str_ends(rownames(ctm), "marmoset")])

ctm = ctm / rowSums(ctm)
ctm = ctm[, colSums(ctm_ori) > 50]

rownames(ctm) = str_c("st-", rownames(ctm))
colnames(ctm) = ifelse(str_starts(colnames(ctm), "EN|IN|CHAT"), toupper(colnames(ctm)), colnames(ctm))


bks_name = str_split_i(rownames(ctm), ' ', 1)
specie_name = str_split_i(rownames(ctm), ' ', 2)
color_map2 = color_map
names(color_map2) = str_c("st-", names(color_map))
ra <- HeatmapAnnotation(
  df = tibble(Specie=specie_name, Region_name=bks_name), 
  col = list(Specie=specie_color, Region_name=color_map2),
  which = "row",
  annotation_legend_param = list(nrow = 2)
)



svglite("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/area.celltype_aboundance.svg", 
        width=15, height=14, fix_text_size=F)
p = Heatmap(ctm, cluster_rows=F, row_split=bks_name, rect_gp = gpar(col = "black", lwd = 1),
              right_annotation = ra, show_row_names = F,
              heatmap_legend_param = list(
                title="Scaled Aboundance",
                direction = "horizontal")
)
draw(p, heatmap_legend_side = "bottom")
dev.off()



############ draw cell type
ctps = bind_rows(list(mous[[]], macs[[]], mars[[]]))
lapply(c("EpC1", "AS2", "OD1"), function(ct) {
  ctps$celltype = ctps$celltype == ct
  pls = lapply(ctps |> 
                mutate(rx=ifelse(specie == "mouse", -rx, rx),
                        ry=-ry) |> group_split(specie), 
  function(df) {
    point_size = case_when(df$specie[1] == "macaque" ~ 0.3, 
                          df$specie[1] == "mouse" ~  0.7,
                          T ~ 0.6)
    ggplot(df, aes(x=rx, y=ry, col=celltype)) +
      geom_point_rast(data = df[!df$celltype,], size=point_size, stroke = 0, color="#EFEEED") +
      geom_point_rast(data = df[df$celltype,], size=point_size, stroke = 0, color="#d63b52") +
      facet_wrap(~chip) +
      guides(colour = guide_colorbar(barwidth = 5, barheight = 25)) +
      NoLegend() +
      coord_fixed() +
      ggtitle(ct) +
      ggspatial::annotation_scale(line_col="black", text_col="black", style="ticks") +
      theme_void() +
      theme(legend.title = element_text(size=30), legend.text = element_text(size=30))
  })
  Reduce('+', pls) + plot_layout(ncol=3, guides = "collect")
  prefix_fn = str_c("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/celltype.plot.", ct)
  ggsave(str_c(prefix_fn, ".svg"), width=45, height=15, fix_text_size=F)
  ggsave(str_c(prefix_fn, ".pdf"), width=45, height=15)
})










############## draw entropy on integration cluster

##### define entropy
etp = function(v) {
  v = (v)/sum(v)
  # print(v)
  -sum(v*log(v))
}

#################### discard one cluster
iseuk$clu_1 = NULL
####################
etp_res = apply(iseuk[[]][, str_c("clu_", 2:50)], 2, function(x) {
  clu_spe = as.data.frame.matrix(table(x, iseuk$specie)) + 1
  clu_spe_norm = apply(clu_spe, 1, function(x) x/colSums(clu_spe))
  e_res = apply(clu_spe_norm, 2, etp)
  e_res[x]
})

################ calculate entropy of mouse vs monky & macaque vs marmoset
mou_mon_score = apply(iseuk[[]][, str_c("clu_", 2:50)], 2, function(x) {
  clu_spe = as.data.frame.matrix(table(x, iseuk$specie)) + 1
  clu_spe_norm = apply(clu_spe, 1, function(x) x/colSums(clu_spe))
  
  monkey=apply(clu_spe_norm[which(rownames(clu_spe_norm) != "mouse"),], 2, mean)
  clu_spe_norm_r = tibble(monkey=monkey,
    mouse=clu_spe_norm[rownames(clu_spe_norm) == "mouse",])
  e_res = apply(clu_spe_norm_r, 1, etp)
  e_res[x]
})
mac_mar_score = apply(iseuk[[]][, str_c("clu_", 2:50)], 2, function(x) {
  clu_spe = as.data.frame.matrix(table(x, iseuk$specie)) + 1
  clu_spe_norm = apply(clu_spe, 1, function(x) x/colSums(clu_spe))
  e_res = apply(t(clu_spe_norm[rownames(clu_spe_norm) != "mouse",]) / colSums(clu_spe_norm[rownames(clu_spe_norm) != "mouse",]), 1, etp)
  e_res[x]
})
iseuk$cell_mou_mon_entropy_mean = rowMeans(mou_mon_score) 
iseuk$cell_mou_mon_entropy_mean_scaled = scale(iseuk$cell_mou_mon_entropy_mean)
iseuk$cell_mac_mar_entropy_mean = rowMeans(mac_mar_score) 
iseuk$cell_mac_mar_entropy_mean_scaled = scale(iseuk$cell_mac_mar_entropy_mean)


xdown = quantile(iseuk$cell_mou_mon_entropy_mean, 0.05)[[1]]
xup = quantile(iseuk$cell_mou_mon_entropy_mean, 0.95)[[1]]

pls = lapply(iseuk[[]] |> 
              mutate(rx=ifelse(specie == "mouse", -rx, rx),
                      ry=-ry) |> group_split(specie), function(df) {
  point_size = case_when(df$specie[1] == "macaque" ~ 0.3, 
                        df$specie[1] == "mouse" ~  0.7,
                        T ~ 0.6)
  ggplot(df, aes(x=rx, y=ry, col=cell_mou_mon_entropy_mean)) +
    geom_point_rast(size=point_size, stroke = 0) +
    facet_wrap(~chip) +
    guides(colour = guide_colorbar(barwidth = 5, barheight = 25)) +
    NoLegend() +
    coord_fixed() +
    scale_color_gradientn(colors=rev(viridis(1e3)), limits = c(xdown, xup), oob = scales::squish) +
    ggspatial::annotation_scale(line_col="black", text_col="black", style="ticks") +
    theme_void() +
    theme(legend.title = element_text(size=30), legend.text = element_text(size=30))
})
Reduce('+', pls) + plot_layout(ncol=3, guides = "collect")
ggsave("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/cell_mou_mon_entropy_mean.svg", width=45, height=15, fix_text_size=F)
ggsave("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/cell_mou_mon_entropy_mean.pdf", width=45, height=15)


xdown = quantile(iseuk$cell_mac_mar_entropy_mean, 0.05)[[1]]
xup = quantile(iseuk$cell_mac_mar_entropy_mean, 0.95)[[1]]


pls = lapply(iseuk[[]] |> 
              mutate(rx=ifelse(specie == "mouse", -rx, rx),
                      ry=-ry) |> group_split(specie), function(df) {
  point_size = case_when(df$specie[1] == "macaque" ~ 0.3, 
                        df$specie[1] == "mouse" ~  0.7,
                        T ~ 0.6)
  ggplot(df, aes(x=rx, y=ry, col=cell_mac_mar_entropy_mean)) +
    geom_point_rast(size=point_size, stroke = 0) +
    facet_wrap(~chip) +
    guides(colour = guide_colorbar(barwidth = 5, barheight = 25)) +
    NoLegend() +
    coord_fixed() +
    scale_color_gradientn(colors=rev(viridis(1e3)), limits = c(xdown, xup), oob = scales::squish) +
    ggspatial::annotation_scale(line_col="black", text_col="black", style="ticks") +
    theme_void() +
    theme(legend.title = element_text(size=30), legend.text = element_text(size=30))
})
Reduce('+', pls) + plot_layout(ncol=3, guides = "collect")
ggsave("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/cell_mac_mar_entropy_mean.svg", width=45, height=15, fix_text_size=F)
ggsave("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/cell_mac_mar_entropy_mean.pdf", width=45, height=15)
##### draw case (chips)
transfer_result_plot_chip_numeric(iseuk[[]], chips, "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/chip.mou_mon.etp.png", 
                        coln="cell_mou_mon_entropy_mean", many_color=viridis(1e3), limits=range(iseuk$cell_mou_mon_entropy_mean))
transfer_result_plot_chip_numeric(iseuk[[]], chips, "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/chip.mac_mar.etp.png", 
                        coln="cell_mac_mar_entropy_mean", many_color=parula(1e3), limits=range(iseuk$cell_mac_mar_entropy_mean))

######### xy plot
offset_region = iseuk[[]] |>
  mutate(rx = ifelse(specie == "mouse", -rx, rx)) |>
  group_by(specie, chip) |>
  mutate(ry = ry + max(iseuk$ry) - max(ry), rx = rx + max(iseuk$rx) - max(rx)) |>
  ungroup() |>
  group_by(specie) |>
  mutate(minx=quantile(rx, 0.01), miny=quantile(ry, 0.01), 
          maxx=quantile(rx, 0.99), maxy=quantile(ry, 0.99),
          ofxo = (rx-minx)/(maxx-minx), ofx = ifelse(ofxo > 0, 1 - ofxo, 0), 
          ofyo = (ry-miny)/(maxy-miny), ofy = ifelse(ofyo > 0, 1 - ofyo, 0),
          specie_count=n()) |>
  ungroup() |>
  group_by(specie, clu6) |>
  summarize(ofx=mean(ofx), ofy=mean(ofy), count=n(), cluster_ratio=sum(count)/(specie_count[1])) |>
  ungroup() |>
  group_by(clu6) |>
  summarize(ofx=mean(ofx), ofy=mean(ofy), cell_mean_ratio=mean(cluster_ratio)) |>
  ungroup() |> 
  mutate(ofx = ofx-min(ofx), ofy = ofy-min(ofy)) |>
  mutate(ofx=ofx/max(ofx), ofy=ofy/max(ofy))
etp_v2_ba = iseuk[[]] |>
  mutate(ratio_mean=apply(mac_mar_score, 1, mean) / apply(mou_mon_score, 1, mean)) |> #ratio_mean=apply(mac_mar_score/mou_mon_score, 1, median)) |>
  group_by(clu6) |>
  summarize(`mouse-NHP entropy`= mean(cell_mou_mon_entropy_mean), 
            `macaque-marmoset entropy`= mean(cell_mac_mar_entropy_mean),
            entropy_ratio = log(median(ratio_mean)),
            NMER = entropy_ratio, # NHP - mouse entropy ratio
            count = n()
            ) |>
  ungroup()
etp_v2_ba = left_join(etp_v2_ba |> filter(clu6 != "HY"), offset_region) |>
  arrange(desc(NMER)) |>
  mutate(clu6=factor(clu6, levels=clu6)) |>
  dplyr::rename(X_offset=ofx, Y_offset=ofy)
ggplot(etp_v2_ba, aes(X_offset, Y_offset,label=clu6)) +
  geom_point(aes(X_offset, Y_offset, size=cell_mean_ratio*1.01), color="black") +
  geom_point(aes(X_offset, Y_offset, color=NMER, size=cell_mean_ratio)) +
  scale_color_gradientn(colors=brewer.rdpu(100), name="NMER\n(NHP-mouse entropy ratio)") +
  geom_text_repel(size=3, seed=42, max.time=10) +
  labs(size = "Region ratio in HY")
ggsave("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/etp.v2.ba.bar.svg", height=4, width=6, fix_text_size=F)
ggsave("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/etp.v2.ba.bar.pdf", height=4, width=6)
ggplot(etp_v2_ba, aes(x=`mouse-NHP entropy`, y=`macaque-marmoset entropy`, 
                        label=clu6, col=clu6)) +  
  geom_point(size=4) +
  geom_text_repel(size=5, force_pull=30, min.segment.length=0, max.overlaps=Inf) +
  scale_color_manual(values=color_map) +
  NoLegend()
ggsave("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/etp.v2.ba.svg", height=4, width=4, fix_text_size=F)

################ calculate per-specie spatial complexity
iseu_diversity = bind_rows(lapply(iseu |> SplitObject("chip"), function(isc) {
  k=50
  isc_pca = isc@reductions$pca@cell.embeddings
  isc_coord = isc[[c("rx", "ry")]]
  nn = nn2(data=as.matrix(isc_coord), k=k, treetype="kd", searchtype="radius", radius=3000, eps=1)$nn.idx
  isc$cell_diversity = apply(nn, 1, function(v) {sum(apply(isc_pca[v,], 2, var))})
  isc[[]]
}))
iseu_diversity = iseu_diversity |>
  group_by(specie) |>
  mutate(cell_diversity_scaled = scale(cell_diversity)) |>
  ungroup()

xdown = quantile(iseu_diversity$cell_diversity_scaled, 0.05)[[1]]
xup = quantile(iseu_diversity$cell_diversity_scaled, 0.95)[[1]]

pls = lapply(iseu_diversity |> 
              mutate(rx=ifelse(specie == "mouse", -rx, rx),
                      ry=-ry) |>
              group_split(specie), function(df) {
  point_size = case_when(df$specie[1] == "macaque" ~ 0.3, 
                        df$specie[1] == "mouse" ~  0.7,
                        T ~ 0.6)
  ggplot(df, aes(x=rx, y=ry, col=cell_diversity_scaled)) +
    geom_point_rast(size=point_size, stroke = 0) +
    facet_wrap(~chip) +
    guides(colour = guide_colorbar(barwidth = 5, barheight = 25)) +
    NoLegend() +
    coord_fixed() +
    scale_color_gradientn(colors=plasma(1e3), limits = c(xdown, xup), oob = scales::squish) +
    ggspatial::annotation_scale(line_col="black", text_col="black", style="ticks") +
    theme_void() +
    theme(legend.title = element_text(size=30), legend.text = element_text(size=30)) #, limits=range(iseu_diversity$cell_diversity))
})
Reduce('+', pls) + plot_layout(ncol=3, guides = "collect")
ggsave("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/cell_diversity.svg", width=45, height=15, fix_text_size=F)
ggsave("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/cell_diversity.pdf", width=45, height=15)

ggplot(iseu_diversity, aes(x=clu6, y=cell_diversity, fill=specie)) + 
  geom_violin(width=1) +
  scale_fill_manual(values=specie_color)  + 
  theme(axis.text.x=element_blank(), # Hides the tick labels (text)
    axis.ticks.x = element_blank()  # Hides the tick marks (small lines)
  ) +
  facet_wrap(~ clu6, ncol=5, scales="free_x")
ggsave("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/cell_diversity.violin.svg", width=10, height=10, fix_text_size=F)


############### sliding window draw of cell diversity across x axis
bin_num=100
swindow_agg_across_x = function(x, v, bin_num) {
  cutoff = quantile(x, probs=seq(1, bin_num)/bin_num)
  grp = sapply(x, function(x) {vc = x > cutoff; rev(names(cutoff)[which(vc)])[1]})
  v |> split(grp) |> sapply(median)
}
idv = iseu_diversity |>
  group_by(specie, chip) |>
  reframe(dv=swindow_agg_across_x(rx, cell_diversity_scaled, bin_num=bin_num), order=seq(1, bin_num-1)/bin_num) |>
  ungroup()
idv = idv |>
  group_by(specie, chip) |>
  mutate(order=ifelse(specie == "mouse", order, rev(order)))

ggplot(idv, aes(order, dv, color=specie)) +
  geom_smooth(method=loess) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_line(arrow = arrow(length = unit(0.4, "cm"), type="closed"))
  ) + 
  scale_color_manual(values=specie_color) +
  xlab("offside position (from third ventricle to ZI)") +
  ylab("Scaled cell diversity across species")
ggsave("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/parcellation/figures/cell_diversity.curve.svg", width=5.5, height=3.5, fix_text_size=F)





####### non-merge result
csmks_not_merge = bind_rows(del_ori) |> # consensus markers
  group_by(cluster) |>
  mutate(cluster_nsp = length(unique(specie))) |>
  ungroup() |>
  filter(!is.na(feature_ortholog)) |>
  filter(p_val_adj<0.05 & avg_log2FC>0.5) |>
  arrange(cluster,feature_ortholog) |>
  group_by(cluster, feature_ortholog, cluster_nsp) |>
  summarize(n=n(), p_val_adj=mean(p_val_adj), avg_log2FC=mean(avg_log2FC)) |>
  ungroup()

write_csv(csmks_not_merge, "mks.region.no_specie_filter.csv")

for (gene in unique(csmks_not_merge$feature_ortholog)) {
  draw_3specie_moug(mseul, gene_name=gene, slot="data", dir_name="gene_plot/marker_clu6")
}




##################### new result
pe_mks = csmks_not_merge_pre |> 
  group_by(cluster, gene_mouse) |>
  filter(all(cluster == "Pe") & all(specie != "mouse") & n() == 2) 
write_csv(tibble(gene=unique(pe_mks$gene)), "pe_mks.csv")

zi_mks = csmks_not_merge_pre |> 
  group_by(cluster, gene_mouse) |>
  filter(all(cluster == "ZI") & all(specie != "mouse") & n() == 2) 
write_csv(tibble(gene=unique(zi_mks$gene)), "ZI_mks.csv")

mbo_mks = csmks_not_merge_pre |> 
  group_by(cluster, gene_mouse) |>
  filter(all(cluster == "MBO") & all(specie != "mouse") & n() == 2) 
write_csv(tibble(gene=unique(mbo_mks$gene)), "MBO_mks.csv")

iseuk$specie2 = ifelse(iseuk$specie == "mouse", "mouse", "NHP")
ZI = iseuk[, iseuk$clu6 == "ZI"]
Idents(ZI) = "specie2"
LHA = iseuk[, iseuk$clu6 == "LHA"]
Pe = iseuk[, iseuk$clu6 == "Pe"]
zi_mks = FindMarkers(ZI, ident.2="mouse", ident.1="NHP", group.by="specie2", only.pos=F) |>
  filter(p_val_adj < 0.05 & avg_log2FC > 1)
lha_mks = FindMarkers(LHA, ident.2="mouse", ident.1="NHP", group.by="specie2", only.pos=F) |>
  filter(p_val_adj < 0.05 & avg_log2FC > 1)
pe_mks = FindMarkers(Pe, ident.2="mouse", ident.1="NHP", group.by="specie2", only.pos=F) |>
  filter(p_val_adj < 0.05 & avg_log2FC > 1)
NHP_mgs = intersect(rownames(zi_mks), rownames(lha_mks))
NHP_mgs_u = union(rownames(zi_mks), rownames(lha_mks))

write_csv(tibble(gene=NHP_mgs), "LHA_ZI_mks.csv")
write_csv(tibble(gene=unique(rownames(zi_mks))), "ZI_mks.csv")
