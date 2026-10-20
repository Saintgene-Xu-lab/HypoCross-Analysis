
if(T){
  library(Seurat)
  library(qs)
  library(viridis)
  library(ComplexHeatmap)
  library(ggplot2)
  library(RColorBrewer)
  library(COSG)
  library(stringr)
}

mar <- qread()
mar$batch <- "1"
mar@meta.data[str_detect(mar$chip,pattern = "T6"),]$batch <- "2"

mar$chip <- mar$chip2
Idents(mar) <- mar$chip

meta <- mar@meta.data
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
mar <- subset(mar,cells = unlist(cells))

mar2 <- SCTransform(mar, vars.to.regress = "batch", verbose = FALSE) 
mar2 <- RunPCA(mar2, verbose = FALSE)
pca_embeddings <- Embeddings(mar2, "pca")

metadata <- mar@meta.data
metadata$sample_id <- metadata$chip
sample_ids <- names(table(metadata$chip))
sample_centroids <- t(sapply(sample_ids, function(s) {
  colMeans(pca_embeddings[rownames(metadata[metadata$sample_id == s,]), ])
}))

# 计算样本间的相关性（基于PCA）
cor_matrix_pca <- cor(t(sample_centroids), method = "pearson")
#cor_matrix_pca[cor_matrix_pca < 0 ] <- 0
pheatmap::pheatmap(cor_matrix_pca,cluster_rows = F,cluster_cols = F, main = "Maraque Sample Correlation (PCA Space)",border_color = "white",color = colorRampPalette(viridis_pal(option = "B")(10))(100))
saveRDS(cor_matrix_pca,file = "mar_v5_cor_plot_2.rds")

