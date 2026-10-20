if(T){
  library(ggplot2)
  library(ggpubr)
  library(Seurat)
  library(reshape2)
  library(cowplot)
  library(qs)
  library(stringr)
}

path <- ""

####### mar
mar1 <- qread("")
mar2 <- readRDS("")
mar1 <- NormalizeData(mar1)
mar2 <- NormalizeData(mar2)

add_umap <- function(st){
  umap <- data.frame(uamp1 = st$rx,umap2 = -st$ry)
  rownames(umap) <- colnames(st)
  DefaultAssay(st) <- "RNA"
  st[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap),key = "umap")
  return(st)
}
mar2 <- add_umap(st = mar2)
DimPlot(mar2,label = T)
mar2 <- NormalizeData(mar2)
FeaturePlot(mar2,features = c("AGRP"))
Idents(mar2) <- mar2$chip

####
mar1_ave <- data.frame(AverageExpression(mar1))
mar2_ave <- data.frame(AverageExpression(mar2))
common_genes <- intersect(rownames(mar1_ave),rownames(mar2_ave))

ave2 <- cbind(mar1_ave[common_genes ,],mar2_ave[common_genes,])
library(ggpubr)
library(ggpmisc)
colnames(ave2) <- c("chip1","chip2")
p <- ggplot(ave2,aes(x = log1p(chip1),y = log1p(chip2)))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(label.y = 5.4)+
  theme_bw()+
  xlim(c(0,4))+
  ylim(c(0,4))
p
ggsave(p,filename = paste0(path,"marmoset_sample_cor_All_slice.pdf"),width = 4,height = 4)
rm(mar1)
rm(mar2)


################################################################################################################
################################### macaque -------------------------------------------------------------------
################################################################################################################
mac <- qread("")
macaque <- readRDS("")
add_umap <- function(st){
  umap <- data.frame(uamp1 = st$rx,umap2 = -st$ry)
  rownames(umap) <- colnames(st)
  DefaultAssay(st) <- "RNA"
  st[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap),key = "umap")
  return(st)
}

macaque <- add_umap(macaque)
DimPlot(macaque,label = T)
macaque <- NormalizeData(macaque)
mac <- NormalizeData(mac)

Idents(mac) <- mac$orig.ident
mac_ave <- data.frame(AverageExpression(mac))

Idents(macaque) <- macaque$orig.ident
mac2_ave <- data.frame(AverageExpression(macaque))
common_genes <- intersect(rownames(mac),rownames(macaque))
ave <- cbind(mac_ave[common_genes,],mac2_ave[common_genes,])
colnames(ave) <- c("chip1","chip2")
p <- ggplot(ave,aes(x = log1p(chip1),y = log1p(chip2)))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(label.y = 5.4)+
  theme_bw()

p
ggsave(p,filename = paste0(path,"macaque_sample_cor_ALL.pdf"),width = 4,height = 4)

################################################################################################################
################################### mouse -------------------------------------------------------------------
################################################################################################################
mou2 <- qread("")
mou_st2 <- readRDS("")

add_umap <- function(st){
  umap <- data.frame(uamp1 = st$rx,umap2 = -st$ry)
  rownames(umap) <- colnames(st)
  DefaultAssay(st) <- "RNA"
  st[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap),key = "umap")
  return(st)
}
mou_st2 <- add_umap(mou_st2)
mou_st2 <- NormalizeData(mou_st2)

mou_st2$orig.ident
Idents(mou_st2) <- mou_st2$orig.ident
mou_st2_ave <- data.frame(AverageExpression(mou_st2))

mou2$orig.ident
Idents(mou2) <- mou2$orig.ident
mou2 <- NormalizeData(mou2)
mou2_ave <- data.frame(AverageExpression(mou2))
common_genes <- intersect(rownames(mou_st2),rownames(mou2))
ave <- cbind(mou2_ave[common_genes,],mou_st2_ave[common_genes,])
colnames(ave) <- c("chip1","chip2")
p <- ggplot(ave,aes(x = log1p(chip1),y = log1p(chip2)))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(label.y = 5.4)+
  theme_bw()

p
ggsave(p,filename = paste0(path,"mou_sample_cor_All.pdf"),width = 4,height = 4)





