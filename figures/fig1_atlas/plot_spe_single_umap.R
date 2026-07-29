if(T){
  library(Seurat)
  library(ggplot2)
  library(qs)
  library(reshape2)
  library(cowplot)
  library(scales)
  library(RColorBrewer)
  display.brewer.all()
  display.brewer.all(type="seq")
  display.brewer.all(type="div")
}


color3  <- c("#e20612","#F57F17","#00b0eb")
names(color3) <-  c("macaque", "marmoset", "mouse")

cell_col <- c("#13a983","#44c1f0","#3f60aa","#e8490f","#f18800","#ffa631","#9ec417","#f47983","#fcbad3","#9b5de5","#cca4e3","#7bcfa6")
names(cell_col) <-  c("GABA","GLU","CHAT","AS","OPC","OD","TC","EpC","MG","END","VLMC","PTC")

mar <- qread("")
meta <- qread("")
meta2 <- meta[colnames(mar),]
mar$donors <- meta2$donors
Idents(mar) <- mar$anno

path <- ""

color4 <- c("#3cb346","#00abf0",
            "#2e409a","#942d8d")

p <- DimPlot(mar,label = T,cols = color4,
             reduction = "umap2",
             group.by = "donors",
             raster = T,alpha = 0.3)
p
ggsave(p,filename = paste0(path,"mar_color2_legend.pdf"),width = 5,height = 5)


p <- DimPlot(mar,label = F,cols = color4,
             reduction = "umap2",
             group.by = "donors",
             raster = F,alpha = 0.1)+theme_nothing()
p
ggsave(p,filename = paste0(path,"mar_color2.pdf"),width = 7,height = 7)

donors <- names(table(meta2$donors))
for (i in 1:length(donors)){
  p <- DimPlot(object = mar,reduction = "umap2",
               cells = rownames(meta2[meta2$donors %in% donors[i],]),cols = cell_col)+theme_nothing()
  ggsave(p,filename = paste0(path,"mar_donor",donors[i],".pdf"),width = 4,height = 4)
}
################ -----------------------------------------------------------------------------------
mac <- qread("")
meta <- qread("")
meta2 <- meta[colnames(mac),]
mac$donors <- meta2$donors
Idents(mac) <- mac$anno

path <- ""
color4 <- c("#3cb346","#00abf0",
            "#2e409a","#942d8d")

p <- DimPlot(mac,label = T,cols = color4,
             reduction = "umap2",
             group.by = "donors",
             raster = T,alpha = 0.3)
p
ggsave(p,filename = paste0(path,"mac_color2_legend.pdf"),width = 5,height = 5)


p <- DimPlot(mac,label = F,cols = color4,
             reduction = "umap2",
             group.by = "donors",
             raster = F,alpha = 0.1)+theme_nothing()
p
ggsave(p,filename = paste0(path,"mac_color2.pdf"),width = 7,height = 7)

donors <- names(table(meta2$donors))
for (i in 1:length(donors)){
  p <- DimPlot(object = mac,reduction = "umap2",
               cells = rownames(meta2[meta2$donors %in% donors[i],]),cols = cell_col)+theme_nothing()
  ggsave(p,filename = paste0(path,"mac_donor",donors[i],".pdf"),width = 4,height = 4)
}

####################### ----------------------------------------------------------------------------
mou <- qread("")
mou$area <- mou$sample_area

p <- DimPlot(mou,label = T,
             reduction = "umap2",
             group.by = "batch",
             raster = T,alpha = 0.3)
p

p <- DimPlot(mou,label = T,
             reduction = "umap2",
             group.by = "sample_area",
             raster = T,alpha = 0.3)
p

color9 <- colorRampPalette(brewer.pal(12,'Paired'))(9)


p <- DimPlot(mou,label = T,
             reduction = "umap2",
             group.by = "sample_area",cols = color9,
             raster = F,alpha = 0.3)+theme_nothing()
p


p <- DimPlot(mou,label = T,cols = color9,
             reduction = "umap2",
             group.by = "sample_area",
             raster = T,alpha = 0.3)
p
ggsave(p,filename = paste0(path,"mou_color2_legend.pdf"),width = 5,height = 5)


p <- DimPlot(mou,label = F,cols = color9,
             reduction = "umap2",
             group.by = "sample_area",
             raster = F,alpha = 0.3)+theme_nothing()
p
ggsave(p,filename = paste0(path,"mou_color2.pdf"),width = 7,height = 7)
