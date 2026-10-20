
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

dir <- './'
nn <- readRDS(paste(dir,"non_neu_male_69142_20250305.rds",sep=''))
mac <- subset(nn,subset=specie=='macaque')
mar <- subset(nn,subset=specie=='marmoset')
mou <- subset(nn,subset=specie=='mouse')

mac_count_count <- as.matrix(mac@assays$RNA@counts)
macaque_count_norm <- log1p(mac_count_count)
macaque.meta <- mac@meta.data

mar_count_count <- as.matrix(mar@assays$RNA@counts)
marmoset_count_norm <- log1p(mar_count_count)
marmoset.meta <- mar@meta.data

mou_count_count <- as.matrix(mou@assays$RNA@counts)
mouse_count_norm <- log1p(mou_count_count)
mouse.meta <- mou@meta.data

rm(nn,mac,mou,mar)
gc()

## Run CoGAPS ref Mouse
result <- CoGAPS(as.matrix(mouse_count_norm), nPatterns=30, nIterations=500, sparseOptimization=TRUE, seed=830, nThreads=15, distributed="single-cell")
saveRDS(result,'./mouse_CoGAPS_result.rds')
myCols <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)

sort_cells <- c()
for(i in sort(unique(mouse.meta$annotation))){ sort_cells <- c(sort_cells, rownames(mouse.meta[mouse.meta$annotation == i,])) }
projectR(data=as.matrix(mouse_count_norm), loadings=result@featureLoadings, dataNames = rownames(mouse_count_norm)) -> projected2mouse
mouse_aucMat <- aucMat(mouse.meta[sort_cells, 'annotation'], projected2mouse[,sort_cells]) 
pdf('1_mouse.pdf')
p <- pheatmap(t(mouse_aucMat), color = myCols, border_color = NA, scale = "column", cluster_cols = T, cluster_rows = T, angle_col = 90, clustering_method = "ward.D2")
dev.off()
sort_cols <- p$tree_col$labels[p$tree_col$order]
pdf('2_mouse.pdf')
pheatmap(t(mouse_aucMat), color = myCols, border_color = NA, scale = "column", cluster_cols = T, cluster_rows = F, angle_col = 90, clustering_method = "ward.D2")
dev.off()

sort_cells <- c()
for(i in sort(unique(macaque.meta$annotation))){ sort_cells <- c(sort_cells, rownames(macaque.meta[macaque.meta$annotation == i,])) }
projectR(data=as.matrix(macaque_count_norm), loadings=result@featureLoadings, dataNames = rownames(macaque_count_norm)) -> projected2macaque
macaque_aucMat <- aucMat(macaque.meta[sort_cells, 'annotation'], projected2macaque[,sort_cells]) 
pdf('3_macaque.pdf')
pheatmap(t(macaque_aucMat)[,sort_cols], color = myCols, border_color = NA, scale = "column", cluster_cols = F, cluster_rows = F, angle_col = 90, clustering_method = "ward.D2")
dev.off()
sort_cells <- c()
for(i in sort(unique(marmoset.meta$annotation))){ sort_cells <- c(sort_cells, rownames(marmoset.meta[marmoset.meta$annotation == i,])) }
projectR(data=as.matrix(marmoset_count_norm), loadings=result@featureLoadings, dataNames = rownames(marmoset_count_norm)) -> projected2marmoset
marmoset_aucMat <- aucMat(marmoset.meta[sort_cells, 'annotation'], projected2marmoset[,sort_cells]) 
pdf('4_marmoset.pdf')
pheatmap(t(marmoset_aucMat)[,sort_cols], color = myCols, border_color = NA, scale = "column", cluster_cols = F, cluster_rows = F, angle_col = 90, clustering_method = "ward.D2")
dev.off()
saveRDS(macaque_aucMat,'./macaque_aucMat_CoGAPS.rds')
saveRDS(mouse_aucMat,'./mouse_aucMat_CoGAPS.rds')
saveRDS(marmoset_aucMat,'./marmoset_aucMat_CoGAPS.rds')



