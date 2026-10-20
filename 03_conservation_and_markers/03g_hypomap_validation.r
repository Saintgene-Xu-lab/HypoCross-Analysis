if(T){
  library(Seurat)
  library(tidyverse)
  library(reshape2)
  library(stringr)
  library(clustree)
  library(ape)
  library(ggplot2)
  library(tidyverse)
  library(reshape2)
  library(qs)
  library(ggtree)
  library(ape)
  library(ggplot2)
  gc()
  library(COSG)
}

### function 
dot_plot_to_heat <- function(seu,gene){
  p <- DotPlot(seu,features = str_to_upper(gene))
  a <- p$data
  a <- a[,c("features.plot","id", "avg.exp.scaled")] %>% 
    spread(key = id,value = avg.exp.scaled)
  a <- data.frame(a)
  rownames(a) <- a$features.plot
  a <- a[,-1] 
  return(a)
}

## ranking 
ranking_mat <- function(data){
  library(zoo)
  ## 基因排序
  window <- 1
  step <- 1
  dat_ts <- data.frame(data)
  dat2_ts <- dat_ts[order(apply(t(rollapply(t(dat_ts), width=window, by=step, FUN=mean)), 1, which.max)), ]
  return(dat2_ts)
}

### 
single_feature <- readRDS("deg_in_single_specie.rds")
single_feature_v2 <- lapply(single_feature, function(x){
  x <- x[x$avg_log2FC > 0.5 & x$p_val_adj < 0.05,]
})


cells <- names(table(single_feature_v2$macaque$cluster))
conservate_cells <- list()
for (i in 1:length(cells)) {
  conservate_cells[[cells[i]]] <- intersect(intersect(single_feature_v2[[1]][single_feature_v2[[1]]$cluster %in% cells[i],]$gene,
                                                      single_feature_v2[[2]][single_feature_v2[[2]]$cluster %in% cells[i],]$gene),single_feature_v2[[3]][single_feature_v2[[3]]$cluster %in% cells[i],]$gene)
}

unique(as.character(unlist(conservate_cells)))


hum_hypomap <- qread("")
hum_neu <- qread("")
hypomap <- readRDS("")

Idents(hypomap) <- hypomap$C2_named
hypomap <- subset(hypomap,idents = c("C2-1: Neurons"))
Idents(hypomap) <- as.character(hypomap$C66_named)

common_genes <- intersect(rownames(hypomap),unique(as.character(unlist(conservate_cells))))
common_genes <- intersect(str_to_upper(common_genes),rownames(hum_hypomap))
common_genes <- intersect(str_to_upper(common_genes),rownames(hum_neu))


Idents(hum_hypomap) <- hum_hypomap$cell


neu_hum_heat <- dot_plot_to_heat(seu = hum_hypomap,gene = common_genes)
neu_hum_heat <- ranking_mat(data = neu_hum_heat)
pheatmap::pheatmap(neu_hum_heat,cluster_rows = F,border_color = "white",
                   cluster_cols = F,show_rownames = F,show_colnames = T,
                   color = colorRampPalette(c("#FAFAFA","#F3E5F5","#AB47BC"))(100))

neu_hum_heat2 <- dot_plot_to_heat(seu = hum_neu,gene = common_genes)
neu_hum_heat2 <- ranking_mat(data = neu_hum_heat2)
pheatmap::pheatmap(neu_hum_heat2,cluster_rows = F,border_color = "white",
                   cluster_cols = F,show_rownames = F,show_colnames = T,
                   color = colorRampPalette(c("#FAFAFA","#F3E5F5","#AB47BC"))(100))

neu_mou_heat <- dot_plot_to_heat(seu = hypomap,gene = common_genes)
neu_mou_heat <- ranking_mat(data = neu_mou_heat)
pheatmap::pheatmap(neu_mou_heat,cluster_rows = F,border_color = "white",
                   cluster_cols = F,show_rownames = F,show_colnames = T,
                   color = colorRampPalette(c("#FAFAFA","#F3E5F5","#AB47BC"))(100))

path <- ""

pheatmap::pheatmap(neu_hum_heat,cluster_rows = F,border_color = "white",filename = paste0(path,"human_hypomap_common.pdf"),
                   cluster_cols = F,show_rownames = F,show_colnames = T,
                   color = colorRampPalette(c("#FAFAFA","#F3E5F5","#AB47BC"))(100))
pheatmap::pheatmap(neu_hum_heat2,cluster_rows = F,border_color = "white",filename = paste0(path,"human_data_common.pdf"),
                   cluster_cols = F,show_rownames = F,show_colnames = T,
                   color = colorRampPalette(c("#FAFAFA","#F3E5F5","#AB47BC"))(100))
pheatmap::pheatmap(neu_mou_heat,cluster_rows = F,border_color = "white",filename = paste0(path,"mouse_hypomap_common.pdf"),
                   cluster_cols = F,show_rownames = F,show_colnames = T,
                   color = colorRampPalette(c("#FAFAFA","#F3E5F5","#AB47BC"))(100))

save(neu_hum_heat,neu_hum_heat2,neu_mou_heat,file = paste0(path,"plot_common_heatmap.Rdata"))
saveRDS(common_genes,file = paste0(path,"common.rds"))
common_genes <- data.frame(common_genes)
colnames(common_genes) <- "conservative features"
writexl::write_xlsx(common_genes,path = paste0(path,"conservative_features.xlsx"))





