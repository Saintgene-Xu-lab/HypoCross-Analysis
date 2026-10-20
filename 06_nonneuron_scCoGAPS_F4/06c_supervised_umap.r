library(spatstat)
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(cowplot)
library(RColorBrewer)
library(CoGAPS)
library(projectR)
library(pheatmap)
library(viridis)
library(ggrepel)
library(umap)
library(ggrepel)

dir <- '~/'

macaque_aucMat <- readRDS('./macaque_aucMat_CoGAPS.rds')
mouse_aucMat <- readRDS('./mouse_aucMat_CoGAPS.rds')
marmoset_aucMat <- readRDS('./marmoset_aucMat_CoGAPS.rds')


rownames(macaque_aucMat) <- paste0("Macaque_", rownames(macaque_aucMat))
rownames(mouse_aucMat) <- paste0("Mouse_", rownames(mouse_aucMat))
rownames(marmoset_aucMat) <- paste0("Marmoset_", rownames(marmoset_aucMat))
macaque_df <- data.frame(macaque_aucMat, Species = "Macaque")
mouse_df <- data.frame(mouse_aucMat, Species = "Mouse")
marmoset_df <- data.frame(marmoset_aucMat, Species = "Marmoset")

combined_mat <- rbind(macaque_df, mouse_df, marmoset_df)
species_labels <- combined_mat$Species  
combined_mat <- combined_mat[, -ncol(combined_mat)]  


############################################################################################################
library(tidyverse)

combined_mat <- bind_rows(
  macaque_aucMat %>% as.data.frame() %>% mutate(Species = "Macaque"),
  mouse_aucMat %>% as.data.frame() %>% mutate(Species = "Mouse"),
  marmoset_aucMat %>% as.data.frame() %>% mutate(Species = "Marmoset")
)

features <- combined_mat %>% select(-Species)
labels <- factor(combined_mat$Species)  
library(uwot)
labels <- factor(species_labels)  
set.seed(321) 

supervised_umap <- umap(
  features,
  y = labels,  
  n_neighbors = 15,
  min_dist = 0.8,
  metric = "euclidean"
)


umap_df <- data.frame(
  UMAP1 = supervised_umap[,1],
  UMAP2 = supervised_umap[,2],
  Species = labels
)

# mac  #E11020 mou #88CAE3   mar  #F2802C 

p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Species)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Macaque" = "#E11020", "Mouse" = "#88CAE3", "Marmoset" = "#F2802C")) +
  labs(title = "Supervised UMAP Projection") +

  theme_minimal(base_size = 10) + 
  theme(
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    axis.text.x = element_text(
      angle = 0, 
      hjust = 1,  
      vjust = 1, 
      size = 12,
      color = "black"
    ),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.5),  # 

    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    legend.position = "right",  #
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(5, 5, 5, 5), "mm"),  #
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) 



library(RColorBrewer)


library(viridis)
colors <- c("#0073C2FF", "#CD534CFF", "#7AA6DCFF", "#003C67FF", "#8F7700FF", "#3B3B3BFF", "#A73030FF", 
   "#20854EFF",  "#984EA3", "#FFFF33",  "#FF9900FF", "#FFCC00FF", "#6699FFFF",
    "#CC33FFFF", "#99991EFF", "#999999", "#66C2A5",
"#F05C3BFF", "#46732EFF", "#71D0F5FF", "#370335FF", "#075149FF", 
 "#00B5E2FF", "#FC8D62", "#8DA0CB", 
"#E78AC3", "#A6D854",  "#E5C494", "#A6CEE3", "#1F78B4",  "#33A02C", "#FB9A99", "#FDBF6F", 
"#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928",  "#D95F02", 
 "#E7298A",  "#E6AB02", "#A6761D", "#666666", "#6F99ADFF", "#FFDC91FF",  
 "#FF00CCFF", "#CC0000FF", "#8C564BFF",  
 "#BCBD22FF", "#17BECFFF", "#D2AF81FF", "#00FF00FF"
)



umap_df$Pattern <- gsub("^_", "", umap_df$Pattern)  #
p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Pattern, shape = Species), size = 2) +  
  geom_mark_ellipse(
    aes(fill = Species), 
    alpha = 0.1, linetype = "dashed", 
    color = NA, show.legend = FALSE  
  ) +
  scale_color_manual(values = colors) +  
  scale_shape_manual(values = c("Macaque" = 16, "Mouse" = 17, "Marmoset" = 15)) +
  scale_fill_manual(values = c("Macaque" = "#E11020", "Mouse" = "#88CAE3", "Marmoset" = "#F2802C")) +
  labs(title = "Supervised UMAP Projection",size=6) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid = element_blank(),
    axis.text = element_text(size = 8),            
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 5)
  )



p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Pattern, shape = Species), size = 2) +  
   geom_mark_hull(  
    aes(fill = Species),  # 
    alpha = 0.1,                        
    concavity = 3,                       
    radius = unit(3, "mm"),              
    show.legend = FALSE                  
  ) +
  scale_color_manual(values = colors) +  
  scale_shape_manual(values = c("Macaque" = 16, "Mouse" = 17, "Marmoset" = 15)) +
  scale_fill_manual(values = c("Macaque" = "#E11020", "Mouse" = "#88CAE3", "Marmoset" = "#F2802C")) +
  labs(title = "Supervised UMAP Projection",size=6) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid = element_blank(),
    axis.text = element_text(size = 8),            
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 5)
  )



pdf('./xu_umap_20260112/8.4_UMAP_w6h5_hull.pdf', h=5, w=6)
print(p)
dev.off()


