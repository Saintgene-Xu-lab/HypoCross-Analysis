if(T){
  library(SeuratWrappers)
  library(Seurat)
  library(BPCells)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(Seurat)
  library(SeuratData)
  library(Azimuth)
  library(ggplot2)
  library(patchwork)
  library(COSG)
  library(RColorBrewer)
  # set this option when analyzing large datasets
  options(future.globals.maxSize = 3e+09)
  options(Seurat.object.assay.version = "v5")
}

path <- ""

hy_neuron <- readRDS("")

meta <- readRDS("")
DefaultAssay(hy_neuron) <- "RNA"
meta <- meta[colnames(hy_neuron),]
hy_neuron
dim(meta)
table(colnames(hy_neuron) == rownames(meta))
hy_neuron$celltype <- hy_neuron$anno6
set.seed(12580)
DefaultAssay(hy_neuron) <- "RNA"

Idents(hy_neuron) <- hy_neuron$anno6
col <- c(colorRampPalette(brewer.pal(12,'Paired'))(32),colorRampPalette(brewer.pal(8,'Set1'))(30))
names(col) <- levels(hy_neuron$anno6)

color3  <- c("#e20612","#F57F17","#00b0eb")
names(color3) <- c("macaque","marmoset","mouse")

p1 <- DimPlot(hy_neuron, reduction = "umap",alpha = 0.1,group.by = "anno6",cols = col)+theme_void()+NoLegend()
ggsave(p1,filename = paste0(path,"/umap_anno6.png"),width = 5,height =5)
ggsave(p1,filename = paste0(path,"/umap_anno6.pdf"),width = 5,height =5)


p2 <- DimPlot(hy_neuron, reduction = "umap",alpha = 0.1,group.by = "anno6",cols = col)
ggsave(p2,filename = paste0(path,"/umap_anno6_2.png"),width = 8,height =5)
ggsave(p2,filename = paste0(path,"/test1_anno6_2.pdf"),width = 8,height =5)

p3 <- DimPlot(hy_neuron, reduction = "umap",alpha = 0.1,label = T,group.by = "specie",cols = color3)
ggsave(p3,filename = paste0(path,"/umap_g_specie.png"),width = 6,height =5)
ggsave(p3,filename = paste0(path,"/test_g_specie.pdf"),width = 6,height =5)

p4 <- DimPlot(hy_neuron, reduction = "umap",alpha = 0.1,label = T,split.by = "specie",cols = color3)+NoLegend()
ggsave(p4,filename = paste0(path,"/umap_s_specie.png"),width = 15,height =5)
ggsave(p4,filename = paste0(path,"/test_s_specie.pdf"),width = 15,height =5)


rm(list=ls())









