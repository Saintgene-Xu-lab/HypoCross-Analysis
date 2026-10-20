if(T){
  library(Seurat)
  library(qs)
  library(viridis)
  library(ComplexHeatmap)
  library(ggplot2)
  library(RColorBrewer)
  library(COSG)
  library(cowplot)
  library(stringr)
}

mar <- qread("")
mac <- qread("")
mou <- qread("")

mar <- NormalizeData(mar)
mou <- NormalizeData(mou)
mac <- NormalizeData(mac)

mar$chip <- as.character(mar$chip2)

Idents(mou) <- mou$chip
Idents(mar) <- mar$chip
Idents(mac) <- mac$chip


plot_gene_location <- function(mac,mar,mou,gene,color){
  ### macaque
  p <- VlnPlot(object = mac,features = gene)
  colnames(p$data)[1] <- "gene"
  cells <- rownames(p$data[p$data$gene > 0.5,])
  p1 <- DimPlot(object = mac,
                cells.highlight = cells,reduction = "umap2",
                sizes.highlight = 0.5,raster = F,cols = "#EEEEEE",cols.highlight = color)+theme_nothing()+
    ggtitle(paste0("Macaque",":",gene,":"))
  
  ### macaque
  p <- VlnPlot(object = mar,features = gene)
  colnames(p$data)[1] <- "gene"
  cells <- rownames(p$data[p$data$gene > 0.5,])
  p2 <- DimPlot(object = mar,
                cells.highlight = cells,reduction = "umap2",
                sizes.highlight = 0.5,raster = F,cols = "#EEEEEE",cols.highlight = color)+theme_nothing()+
    ggtitle(paste0("Marmoset",":",gene,":"))
  
  ### mou
  p <- VlnPlot(object = mou,features = str_to_title(gene))
  colnames(p$data)[1] <- "gene"
  cells <- rownames(p$data[p$data$gene > 0.5,])
  p3 <- DimPlot(object = mou,
                cells.highlight = cells,reduction = "umap2",
                sizes.highlight = 0.5,raster = F,cols = "#EEEEEE",cols.highlight = color)+theme_nothing()+
    ggtitle(paste0("Mouse",":",str_to_title(gene),":"))
  
  return(p1+p2+p3)
  gc()
  
}

### macaque_specific
path <- ""
p <- plot_gene_location(mac = mac,mar = mar,mou = mou,gene = "TERB1",color = "#fd624d")
ggsave(p,file = paste0(path,"mac_TERB1.pdf"),width = 27,height = 8)

### primate
p <- plot_gene_location(mac = mac,mar = mar,mou = mou,gene = "BBX",color = "#97AAE8")
p
ggsave(p,file = paste0(path,"pri_BBX.pdf"),width = 27,height = 8)
ggsave(p,file = paste0(path,"pri_BBX.png"),width = 27,height = 8)

### mar 

p <- plot_gene_location(mac = mac,mar = mar,mou = mou,gene = "PEG10",color = "#FF9523")
p
ggsave(p,file = paste0(path,"mar_PEG10.pdf"),width = 27,height = 8)
ggsave(p,file = paste0(path,"mar_PEG10.png"),width = 27,height = 8)

p <- plot_gene_location(mac = mac,mar = mar,mou = mou,gene = str_to_upper("Rbfox1"),color = "#FF9523")
p
ggsave(p,file = paste0(path,"mar_INPP5F.pdf"),width = 27,height = 8)
ggsave(p,file = paste0(path,"mar_INPP5F.png"),width = 27,height = 8)

p <- plot_gene_location(mac = mac,mar = mar,mou = mou,gene = str_to_upper("ACRBP"),color = "#FF9523")

ggsave(p,file = paste0(path,"mar_ACRBP.pdf"),width = 27,height = 8)
ggsave(p,file = paste0(path,"mar_ACRBP.png"),width = 27,height = 8)

### mouse
p <- plot_gene_location(mac = mac,mar = mar,mou = mou,gene = "CD5L",color = "#FF9523")
p
ggsave(p,file = paste0(path,"mou_CD5L.pdf"),width = 27,height = 8)
ggsave(p,file = paste0(path,"mou_CD5L.png"),width = 27,height = 8)
