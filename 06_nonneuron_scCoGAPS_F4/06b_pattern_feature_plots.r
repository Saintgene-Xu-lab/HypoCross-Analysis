
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


dir_cogapsresult <- './scCoGAPS_allGlial_Version2_byMou/'
result <- readRDS(paste(dir_cogapsresult,'mouse_CoGAPS_result.rds',sep=''))

sort_cells <- c()
for(i in sort(unique(mouse.meta$annotation))){ sort_cells <- c(sort_cells, rownames(mouse.meta[mouse.meta$annotation == i,])) }
projectR(data=as.matrix(mouse_count_norm), loadings=result@featureLoadings, dataNames = rownames(mouse_count_norm)) -> projected2mouse
mouse_aucMat <- aucMat(mouse.meta[sort_cells, 'annotation'], projected2mouse[,sort_cells]) 
saveRDS(projected2mouse,'projected2mouse.RDS')

sort_cells <- c()
for(i in sort(unique(macaque.meta$annotation))){ sort_cells <- c(sort_cells, rownames(macaque.meta[macaque.meta$annotation == i,])) }
projectR(data=as.matrix(macaque_count_norm), loadings=result@featureLoadings, dataNames = rownames(macaque_count_norm)) -> projected2macaque
macaque_aucMat <- aucMat(macaque.meta[sort_cells, 'annotation'], projected2macaque[,sort_cells]) 
saveRDS(projected2macaque,'projected2macaque.RDS')

sort_cells <- c()
for(i in sort(unique(marmoset.meta$annotation))){ sort_cells <- c(sort_cells, rownames(marmoset.meta[marmoset.meta$annotation == i,])) }
projectR(data=as.matrix(marmoset_count_norm), loadings=result@featureLoadings, dataNames = rownames(marmoset_count_norm)) -> projected2marmoset
marmoset_aucMat <- aucMat(marmoset.meta[sort_cells, 'annotation'], projected2marmoset[,sort_cells]) 

saveRDS(projected2marmoset,'projected2marmoset.RDS')


# ========================================================================================

dir_cogapsresult <- './scCoGAPS_allGlial_Version2_byMou/'
cogapsresult <- readRDS(paste(dir_cogapsresult,'mouse_CoGAPS_result.rds',sep=''))
mac_result <- readRDS('projected2macaque.RDS')
mar_result <- readRDS('projected2marmoset.RDS')

dir <- './'
nn <- readRDS(paste(dir,"non_neu_male_69142_20250305.rds",sep=''))
speci <- 'macaque'
af <- subset(nn,subset=specie==speci)

# mou_af <- subset(nn,subset=specie=='mouse')
# patterns_in_order <-t(cogapsresult@sampleFactors[colnames(mou_af),])
af[["CoGAPS"]] <- CreateAssayObject(counts = mar_result)
DefaultAssay(af) <- "CoGAPS"
pattern_names <- rownames(af@assays$CoGAPS)
color_palette <- viridis(n=10)
pdf(paste0(speci,'FeaturePlot_pattern_names_1.pdf'),height=40,width=25)
FeaturePlot(af, pattern_names, cols=color_palette, reduction = "umap") & NoLegend()
dev.off()





create_clean_feature_plot <- function(seurat_obj, feature,
                                      pt_size = 0.1, col_palette = "viridis") {
  
  # 基础绘图
  p <- FeaturePlot(
    object = seurat_obj,
    features = feature,
    pt.size = pt_size,           
    order = TRUE,                
    combine = TRUE,             
    cols = col_palette           
  )
  
  p <- p + theme_classic() + 
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      
      legend.position = "right",
      
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      
      plot.margin = unit(c(-0.5, -0.5, -0.5, -0.5), "cm"),
      aspect.ratio = 1 
    )
  
  p <- p + coord_fixed(
    xlim = c(min(seurat_obj@reductions$umap@cell.embeddings[,1]), 
             max(seurat_obj@reductions$umap@cell.embeddings[,1]) * 1.05),
    ylim = c(min(seurat_obj@reductions$umap@cell.embeddings[,2]),
             max(seurat_obj@reductions$umap@cell.embeddings[,2]))
  )
  
  return(p)
}



for (i in pattern_names){
  print(i)

  clean_plot <- create_clean_feature_plot(
    seurat_obj = af,
    feature = i,      
    pt_size = 0.05,              
    col_palette = color_palette
  )

  ggsave(paste("./pattern_comp_FeaturePlot_clean" ,"_",speci, "/fp_", i ,"_",speci, "_20260109.pdf",sep=""), clean_plot,
         width = 3, height = 3, dpi = 300, units = "in",
         bg = "transparent") 
}




