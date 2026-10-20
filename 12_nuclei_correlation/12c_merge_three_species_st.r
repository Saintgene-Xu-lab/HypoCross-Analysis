if(T){
  library(Seurat)
  library(qs)
  library(viridis)
  # 查看画板颜色
  library(scales)
  library(tidyverse)
  library(ggplot2)
  library(reshape2)
}



mac <- qread("mac.qs")
set.seed(12580)
mac <- subset(mac,cells = sample(colnames(mac),size = 34000,replace = F) )

mar <- qread("mar.qs")
set.seed(12580)
mar <- subset(mar,cells = sample(colnames(mar),size = 27000,replace = F) )

mou <- qread("mou2.qs")
set.seed(12580)
mou <- subset(mou,cells = sample(colnames(mou),size = 31000,replace = F) )

genes <- intersect(str_to_upper(rownames(mou)),rownames(mac))
genes <- intersect(genes,rownames(mar))
length(genes)

mac <- CreateSeuratObject(counts = mac@assays$RNA@counts[genes,],meta.data = mac@meta.data)
mar <- CreateSeuratObject(counts = mar@assays$RNA@counts[genes,],meta.data = mar@meta.data)
mar$chip <- mar$chip2

mat <- mou@assays$RNA@counts
mat <- mat[str_to_title(genes),]
table(rownames(mat) == str_to_title(genes))
rownames(mat) <- genes
mou <- CreateSeuratObject(counts =mat,meta.data = mou@meta.data)

mac$spe <- "mac"
mar$spe <- "mar"
mou$spe <- "mou"

int <- merge(mou,y = list(mar,mac))
int <- JoinLayers(int)
qsave(int,file = "int_three_spe.qs")
