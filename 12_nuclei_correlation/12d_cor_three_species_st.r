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

int <- qread("int_three_spe.qs")

int$group <- paste0(int$spe,"_",int$chip)
Idents(int) <- int$group
ave <- data.frame(AverageExpression(int))
hvgs <- VariableFeatures(int)
a <- scale(ave[hvgs,])
cor_dat <- cor(a)
pheatmap::pheatmap(cor_dat,cluster_rows = F,
                   cluster_cols = F, 
                   main = "ample Correlation (PCA Space)",
                   border_color = "white",
                   color = colorRampPalette(viridis_pal(option = "B")(10))(100))
