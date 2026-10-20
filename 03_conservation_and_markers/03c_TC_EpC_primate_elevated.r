library(Seurat) 
library(dplyr)

dir('~/')
nn <- readRDS(paste(dir,"non_neu_male_69142_20250305.rds",sep=''))

tc <- subset(nn ,subset=major_class== 'TC')
species_list <- unique(tc$specie)  
library(dplyr)
library(ggplot2)
library(princurve)

results <- lapply(species_list, function(sp) {
  
  tc_sp <- subset(tc, subset = specie == sp)
  
  tc_sp <- FindVariableFeatures(tc_sp, nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE,npcs=26) %>%
    FindNeighbors(dims = 1:10) %>%
    RunUMAP(dims = 1:10)
  
  pcembed <- Embeddings(tc_sp, "pca")[, 1:2]
  y <- principal_curve(pcembed)
  
  df <- data.frame(y$s[order(y$lambda), ], species = sp)
  colnames(df)[1:2] <- c("x", "y")
  
  points <- data.frame(
    id = tc_sp$annotation,
    PC_1 = pcembed[, 1],
    PC_2 = pcembed[, 2],
    species = sp
  )
  
  list(
    plot_data = list(curve = df, points = points),
    object = tc_sp,
    principal_curve = y
  )
})
names(results) <- species_list
all_curve <- do.call(rbind, lapply(results, function(x) x$plot_data$curve))
all_points <- do.call(rbind, lapply(results, function(x) x$plot_data$points))
combined_plot <- ggplot(all_curve, aes(x, y)) +
  geom_point(
    data = all_points,
    aes(PC_1, PC_2, color = id),
    alpha = 0.5, size = 1.5
  ) +
  geom_line(
    aes(group = species),
    arrow = arrow(type = "closed", length = unit(0.2, "cm")),
    color = "orange", linewidth = 1
  ) +
  facet_wrap(~species, scales = "free") +
  theme_pubr(legend = "none") +
  scale_color_jco() +
  labs(x = "PC1", y = "PC2") +
  coord_flip()
setwd('../TC/')
ggsave("combined_species_tc_3and3.png", combined_plot, width = 3, height = 3)


#### ***************************************



cts <- c("TCα", "TCβ", "EpC2", "EpC1")
st_sc_li_1 <- list()

for (ct in cts){
  conserved_markers_TCa <- readRDS(paste('./conserved_markers_' , ct,'.RDS',sep=''))
  
  conserved_markers_TCa <- conserved_markers_TCa %>%
    tibble::rownames_to_column("gene")  
  

  higher_in_mm <- conserved_markers_TCa %>%
    filter(
      macaque_p_val_adj < 0.05 &   
        marmoset_p_val_adj < 0.05 &  
        macaque_avg_log2FC > mouse_avg_log2FC,  
      marmoset_avg_log2FC > mouse_avg_log2FC   
    )
  
  higher_in_mm <- higher_in_mm %>% 
    mutate(
      macaque_diff = macaque_avg_log2FC - mouse_avg_log2FC,
      marmoset_diff = marmoset_avg_log2FC - mouse_avg_log2FC,
      total_diff = macaque_diff + marmoset_diff
    ) %>% 
    arrange(desc(total_diff)) 
  
  head(higher_in_mm)
    top_50_genes <- higher_in_mm %>%
    slice_max(total_diff, n = 25) 
  
  heatmap_data <- top_50_genes[, c("mouse_avg_log2FC", "macaque_avg_log2FC", "marmoset_avg_log2FC")]
  rownames(heatmap_data) <- top_50_genes$gene
  names(heatmap_data) <- c('mouse', 'macaque','marmoset')
  p <- pheatmap(
    heatmap_data,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,  
    fontsize_row = 10,
    fontsize_col =15,
    main = paste0("Log2FC of Selected Genes ",ct)
  )
  ggsave(paste("heatmap_",ct,"_.pdf",sep=''), p, width = 3, height = 5)
  
  top_50_genes$gene
  st_sc_li_1[[ct]] <- st_marker[st_marker$gene %in% toupper(top_50_genes$gene),]
}


library(RColorBrewer)
library(scales)st_sc_tca <- st_sc_li_1$TCα[st_sc_li_1$TCα$celltype_pred == 'TCα',]
st_sc_tcb <- st_sc_li_1$TCβ[st_sc_li_1$TCβ$celltype_pred == 'TCβ',]

st_sc_EpC1 <- st_sc_li_1$EpC1[st_sc_li_1$EpC1$celltype_pred == 'EpC1',]
st_sc_EpC2 <- st_sc_li_1$EpC2[st_sc_li_1$EpC2$celltype_pred == 'EpC2',]
dput(st_sc_EpC2$gene)
########

capitalize_first_letter <- function(x) {
  sapply(x, function(s) {
    if (nchar(s) == 0) return(s)  
    paste0(
      toupper(substr(s, 1, 1)),  
      tolower(substr(s, 2, nchar(s)))  
    )
  }, USE.NAMES = FALSE) 
}


my_vector <-st_sc_EpC2$gene


result <- capitalize_first_letter(my_vector)


dput(result)




