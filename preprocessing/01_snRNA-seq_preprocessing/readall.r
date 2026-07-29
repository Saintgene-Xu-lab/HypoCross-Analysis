#! /usr/bin/env R

# read all snRNA data

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
  library(ggspatial)

  # specific lib
})
###########################################
source("/public/group/weiwu_group/majianyun/data/Hypothalamus/draw_spatial_lib.r")
setwd("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial")
options(future.globals.maxSize=+Inf)
# theme set
theme_set(theme_light())
plan(multicore, workers=8)

# read snRNA
mouse = readRDS("/public/group/weiwu_group/majianyun/data/Hypothalamus/raw.mouse.snRNA/20221128_czh/rds/mouse_czh_own_doublet_rm.rf1.rds")
macaque = readRDS("/public/group/weiwu_group/majianyun/data/Hypothalamus/raw.macaque.snRNA/20221214_all/rds/macaque_all_doublet_rm.rf1.rds")
marmoset = readRDS("/public/group/weiwu_group/majianyun/data/Hypothalamus/raw.marmoset.snRNA/20221201_all/rds/marmoset_all_doublet_rm.rf1.rds")

# read all data (202505)
aseu = readRDS("/cluster/share/hypo_group/snRNA_20250306_clena_rm_female/all_cell_male_677320_20250305.rds")
aseu$annotation = aseu[[]]$annotation
mg = aseu[, str_starts(aseu[[]]$annotation, "MG")]
mg = mg |> NormalizeData()
aseu$MG = str_starts(aseu[[]]$annotation, "MG")
aseu$HTR2C = aseu@assays$RNA@layers$data[which(rownames(aseu) == "Htr2c"),]
aseu$ms = str_c(aseu$MG, aseu$specie)
aseu = aseu |> NormalizeData()

hms = split(aseu$HTR2C, aseu$ms)
wilcox.test(hms$FALSEmacaque, hms$TRUEmacaque, alternative="less")
wilcox.test(hms$FALSEmouse, hms$TRUEmouse, alternative="less")
wilcox.test(hms$FALSEmarmoset, hms$TRUEmarmoset, alternative="less")
VlnPlot(aseu, "Htr2c", group.by="specie", split.by="MG", pt.size=0)

VlnPlot(mg, "Htr2c", group.by="specie", pt.size=0)
ggsave("~/t.jpg", height=3)


aseu = readRDS("/cluster/share/hypo_group/snRNA_20250306_clena_rm_female/all_cell_male_677320_20250305.rds")
as = aseu[, str_starts(aseu[[]]$annotation, "AS")]
as = as |> NormalizeData()
VlnPlot(as, "Hspb8", group.by="annotation", split.by="specie", pt.size=0) +
  scale_fill_manual(values=c("mouse"="#00B0EB", "macaque"="#E2424A", "marmoset"="#F57F17", "human"="#2fd897"))
ggsave("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/label_transfer_evaluation/HSPB8.snRNA.svg", width=10, height=3, fix_text_size=F)

# read spatial 
mous = readRDS("/cluster/share/hypo_group/Spatial/spatial_rds_20250328_with_anno/mouse.half.rds")
macs = readRDS("/cluster/share/hypo_group/Spatial/spatial_rds_20250328_with_anno/macaque.half.rds")
mars = readRDS("/cluster/share/hypo_group/Spatial/spatial_rds_20250328_with_anno/marmoset.half.rds") 
Idents(mous) = "mainAreaID"
Idents(macs) = "mainAreaID"
Idents(mars) = "mainAreaID"
mous_mks = FindAllMarkers(mous, only.pos=T)
macs_mks = FindAllMarkers(macs, only.pos=T)
mars_mks = FindAllMarkers(mars, only.pos=T)
write_csv(mous_mks, "mouse.half.mainAreaID.mks.csv")
write_csv(macs_mks, "macaque.half.mainAreaID.mks.csv")
write_csv(mars_mks, "marmoset.half.mainAreaID.mks.csv")
### draw HSPB8
####### AS
draw_gene_expression_density(mous[, mous$chip == "T322"], "Hspb8", 
      "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/label_transfer_evaluation/HSPB8", 
      gname_replace="Hspb8",  specie_name="mou", slot="counts", extension="svg", clu_prefix="AS")
draw_gene_expression_density(macs[, macs$chip == "T794"], "HSPB8", #T800
      as.character(str_glue("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/label_transfer_evaluation/HSPB8")), 
      gname_replace="HSPB8",  specie_name="mac", slot="counts", extension="svg", clu_prefix="AS")
draw_gene_expression_density(mars[, mars$chip == "T012_T623"], "HSPB8", 
      as.character(str_glue("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/label_transfer_evaluation/HSPB8")), 
      gname_replace="HSPB8",  specie_name="mar", slot="counts", extension="svg", clu_prefix="AS")
sort(AverageExpression(macs, feature="HSPB8",group.by="celltype_pred")$SCT)


### draw HTR2C
draw_gene_expression_density(mous[, mous$chip == "T322"], "Htr2c", 
      "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/label_transfer_evaluation/HTR2C", 
      gname_replace="HTR2C",  specie_name="mou", slot="data", extension="svg", clu_prefix="MG")
draw_gene_expression_density(macs[, macs$chip == "T800"], "HTR2C", 
      as.character(str_glue("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/label_transfer_evaluation/HTR2C")), 
      gname_replace="HTR2C",  specie_name="mac", slot="data", extension="svg", clu_prefix="MG")
draw_gene_expression_density(mars[, mars$chip == "T012_T623"], "HTR2C", 
      as.character(str_glue("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/label_transfer_evaluation/HTR2C")), 
      gname_replace="HTR2C",  specie_name="mar", slot="data", extension="svg", clu_prefix="MG")



### draw all
draw_gene_expression_density(mous[, mous$chip == "T309"], "Hspb8", 
      "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/label_transfer_evaluation/HSPB8", 
      gname_replace="Hspb8",  specie_name="mou", slot="data", extension="svg")
draw_gene_expression_density(macs[, macs$chip == "T786"], "HSPB8", 
      as.character(str_glue("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/label_transfer_evaluation/HSPB8")), 
      gname_replace="HSPB8",  specie_name="mac", slot="data", extension="svg")
draw_gene_expression_density(mars[, mars$chip == "T011_T624"], "HSPB8", 
      as.character(str_glue("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial/label_transfer_evaluation/HSPB8")), 
      gname_replace="HSPB8",  specie_name="mar", slot="data", extension="svg")





for (i in c("P2RY12")) {
  draw_gene_expression(mous, str_to_title(i), "~/mouse", extension="png")
  draw_gene_expression(macs, toupper(i), "~/macaque", extension="png")
  draw_gene_expression(mars, toupper(i), "~/marmoset", extension="png")
}





all_df = bind_rows(mouse[[]], macaque[[]], marmoset[[]])
rdf = all_df |>
    mutate(cutoff=nFeature_RNA > 1500) |>
    group_by(batch) |>
    summarise(cut_ratio=sum(cutoff)/n()) |>
    ungroup() |>
    arrange(cut_ratio)
print(rdf, n=100)
ln = rdf |>
  filter(cut_ratio<0.3)
ln = ln$batch
ggplot(all_df %>% filter(batch %in% ln), aes(x=batch, y=nFeature_RNA)) +
  geom_violin()
ggsave("~/batch.vln.png", width=14)
macaque$lb = macaque$batch %in% ln
FeaturePlot(macaque, "lb", raster=F)
ggsave("~/lb.png")
ggplot(macaque[[]] |> filter(lb), aes(x=seurat_clusters)) + geom_bar()
ggsave("~/lb.bar.png", height=3)


# plot spatial nFeature/nCount

p1=ggplot(mous[[]], aes(y=nFeature_RNA, x=mainAreaID, fill=after_stat(n))) +
  geom_violin() +
  coord_flip()
p2=ggplot(macs[[]], aes(y=nFeature_RNA, x=mainAreaID, fill=after_stat(n))) +
  geom_violin() +
  coord_flip() 
p3=ggplot(mars[[]], aes(y=nFeature_RNA, x=mainAreaID, fill=after_stat(n))) +
  geom_violin() +
  coord_flip()

p1 | p2 | p3
ggsave("~/j.jpg", height=10, width=15)
