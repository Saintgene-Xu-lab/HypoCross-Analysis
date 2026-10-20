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
setwd("/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial")
options(future.globals.maxSize=+Inf)
# theme set
theme_set(theme_light())
plan(multicore, workers=8)

# read snRNA
mouse = readRDS("/public/group/weiwu_group/majianyun/data/Hypothalamus/raw.mouse.snRNA/20221128_czh/rds/mouse_czh_own_doublet_rm.rf1.rds")
macaque = readRDS("/public/group/weiwu_group/majianyun/data/Hypothalamus/raw.macaque.snRNA/20221214_all/rds/macaque_all_doublet_rm.rf1.rds")
marmoset = readRDS("/public/group/weiwu_group/majianyun/data/Hypothalamus/raw.marmoset.snRNA/20221201_all/rds/marmoset_all_doublet_rm.rf1.rds")

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
