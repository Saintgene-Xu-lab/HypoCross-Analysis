
if(T){
  library(Seurat)
  library(qs)
  library(viridis)
  library(ComplexHeatmap)
  library(ggplot2)
  library(RColorBrewer)
  library(COSG)
}

mar <- qread("mar.qs")
mac <- qread("mac.qs")
mou <- qread("mou2.qs")

################ v4  macaque -----------------------------------------------------------
meta <- mac@meta.data
sample_Cell_for_chip <- function(dat){
  chip <- names(table(dat$chip))
  cells <- list()
  for (i in 1:length(chip)) {
    set.seed(12580)
    cells[[chip[i]]] <- sample(rownames(dat[dat$chip %in% chip[i],]),size = 1000,replace = F)
    
  }
  return(cells)
  
  
}
cells <- sample_Cell_for_chip(dat = meta)
mac <- subset(mac,cells = unlist(cells))

get_hvgs <- function(dat,hvgs){
  
  chip <- names(table(dat$chip))
  cells <- list()
  for (i in 1:length(chip)) {
    set.seed(12580)
    dat2 <- subset(dat,idents = chip[i])
    dat2 <- NormalizeData(dat2)
    dat2 <- FindVariableFeatures(dat2,nfeatures = hvgs)
    cells[[chip[i]]] <- VariableFeatures(dat2)
    
  }
  return(cells)
  
  
  
}

a <- get_hvgs(dat = mac,hvgs = 100)
a <- unique(unlist(a))
get_empeading2 <- function(seu,gene){
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu,features = gene)
  seu <- RunPCA(seu)
  return(seu)
}
mac <- get_empeading2(seu = mac,gene = a)

pca_embeddings <- Embeddings(mac, "pca")

metadata <- mac@meta.data
metadata$sample_id <- metadata$chip
sample_ids <- names(table(metadata$chip))
sample_centroids <- t(sapply(sample_ids, function(s) {
  colMeans(pca_embeddings[rownames(metadata[metadata$sample_id == s,]), ])
}))

# 计算样本间的相关性（基于PCA）
cor_matrix_pca <- cor(t(sample_centroids), method = "pearson")
#cor_matrix_pca[cor_matrix_pca < 0 ] <- 0
pheatmap::pheatmap(cor_matrix_pca,cluster_rows = F,cluster_cols = F, main = "Macaque Sample Correlation (PCA Space)",border_color = "white",color = colorRampPalette(viridis_pal(option = "B")(10))(100))
qsave(mac,file = "mac_v4_cor.qs")
saveRDS(cor_matrix_pca,file = "mac_v4_cor_plot.rds")

################### mouse 

meta <- mou@meta.data
sample_Cell_for_chip <- function(dat){
  chip <- names(table(dat$chip))
  cells <- list()
  for (i in 1:length(chip)) {
    set.seed(12580)
    cells[[chip[i]]] <- sample(rownames(dat[dat$chip %in% chip[i],]),size = 700,replace = F)
    
  }
  return(cells)
  
  
}
cells <- sample_Cell_for_chip(dat = meta)
mou <- subset(mou,cells = unlist(cells))

get_hvgs <- function(dat,hvgs){
  
  chip <- names(table(dat$chip))
  cells <- list()
  for (i in 1:length(chip)) {
    set.seed(12580)
    dat2 <- subset(dat,idents = chip[i])
    dat2 <- NormalizeData(dat2)
    dat2 <- FindVariableFeatures(dat2,nfeatures = hvgs)
    cells[[chip[i]]] <- VariableFeatures(dat2)
    
  }
  return(cells)
  
  
  
}

a <- get_hvgs(dat = mou,hvgs = 100)
a <- unique(unlist(a))
get_empeading2 <- function(seu,gene){
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu,features = gene)
  seu <- RunPCA(seu)
  return(seu)
}
mou <- get_empeading2(seu = mou,gene = a)

pca_embeddings <- Embeddings(mou, "pca")

metadata <- mou@meta.data
metadata$sample_id <- metadata$chip
sample_ids <- names(table(metadata$chip))
sample_centroids <- t(sapply(sample_ids, function(s) {
  colMeans(pca_embeddings[rownames(metadata[metadata$sample_id == s,]), ])
}))

# 计算样本间的相关性（基于PCA）
cor_matrix_pca <- cor(t(sample_centroids), method = "pearson")
#cor_matrix_pca[cor_matrix_pca < 0 ] <- 0
pheatmap::pheatmap(cor_matrix_pca,cluster_rows = F,cluster_cols = F, main = "mouse Sample Correlation (PCA Space)",border_color = "white",color = colorRampPalette(viridis_pal(option = "B")(10))(100))

qsave(mou,file = "mou_v4_cor.qs")
saveRDS(cor_matrix_pca,file = "mou_v4_cor_plot.rds")
