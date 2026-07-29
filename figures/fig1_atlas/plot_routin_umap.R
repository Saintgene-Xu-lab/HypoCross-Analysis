
path <- ""
if(T){
  library(Seurat)
  library(ggplot2)
  library(reshape2)
}

dat <- qread(file = "")

mou <- qread("")
mou <- subset(mou,cells = intersect(rownames(dat),colnames(mou)))
qsave(mou,file = "")

mar <- qread("")
mac <- qread("")


############## -----------
all_cell$anno2 <- factor(all_cell$anno2,levels = c("GABA","GLU","CHAT","AS","OPC","OD","TC","EpC","MG","END","VLMC","PTC"))
Idents(all_cell) <- all_cell$anno2




if(T){
  library(Seurat)
  library(ggplot2)
  library(reshape2)
}
library(scales)
library(cowplot)
path <- ""
DimPlot(mou,label = T,reduction = "umap2",cols = color,raster = F,alpha = 0.1)+NoLegend()

### color set
color <- c("#13a983","#44c1f0","#3f60aa","#e8490f","#f18800","#ffa631","#9ec417","#f47983","#fcbad3","#9b5de5","#cca4e3","#7bcfa6")
names(color) <- c("GABA","GLU","CHAT","AS","OPC","OD","TC","EpC","MG","END","VLMC","PTC")

### data 
mou <- qread("")
mar <- qread("")
mac <- qread("")

p1 <- DimPlot(mou,label = T,reduction = "umap2",cols = color,raster = F,alpha = 0.2)+theme_void()+NoLegend()
p1
ggsave(p1,filename = paste0(path,"umap_mou.png"),width = 4,height =4)
ggsave(p1,filename = paste0(path,"umap_mou.pdf"),width = 4,height =4)

p2 <- DimPlot(mou,label = F,reduction = "umap2",cols = color,raster = F,alpha = 0.2)+theme_void()+NoLegend()
p2
ggsave(p2,filename = paste0(path,"umap_mou2.png"),width = 4,height =4)
ggsave(p2,filename = paste0(path,"umap_mou2.pdf"),width = 4,height =4)

p3 <- DimPlot(mou,label = F,reduction = "umap2",cols = color,raster = F,alpha = 0.1)+theme_void()+NoLegend()
p3
ggsave(p3,filename = paste0(path,"umap_mou3.png"),width = 4,height =4)
ggsave(p3,filename = paste0(path,"umap_mou3.pdf"),width = 4,height =4)

####### ----- mar

p1 <- DimPlot(mar,label = T,reduction = "umap2",cols = color,raster = F,alpha = 0.2)+theme_void()+NoLegend()
p1
ggsave(p1,filename = paste0(path,"umap_mar.png"),width = 4,height =4)
ggsave(p1,filename = paste0(path,"umap_mar.pdf"),width = 4,height =4)

p2 <- DimPlot(mar,label = F,reduction = "umap2",cols = color,raster = F,alpha = 0.2)+theme_void()+NoLegend()
p2
ggsave(p2,filename = paste0(path,"umap_mar2.png"),width = 4,height =4)
ggsave(p2,filename = paste0(path,"umap_mar2.pdf"),width = 4,height =4)

p3 <- DimPlot(mar,label = F,reduction = "umap2",cols = color,raster = F,alpha = 0.1)+theme_void()+NoLegend()
p3
ggsave(p3,filename = paste0(path,"umap_mar3.png"),width = 4,height =4)
ggsave(p3,filename = paste0(path,"umap_mar3.pdf"),width = 4,height =4)

### mac 
p1 <- DimPlot(mac,label = T,reduction = "umap2",cols = color,raster = F,alpha = 0.2)+theme_void()+NoLegend()
p1
ggsave(p1,filename = paste0(path,"umap_mac.png"),width = 4,height =4)
ggsave(p1,filename = paste0(path,"umap_mac.pdf"),width = 4,height =4)

p2 <- DimPlot(mac,label = F,reduction = "umap2",cols = color,raster = F,alpha = 0.2)+theme_void()+NoLegend()
p2
ggsave(p2,filename = paste0(path,"umap_mac2.png"),width = 4,height =4)
ggsave(p2,filename = paste0(path,"umap_mac2.pdf"),width = 4,height =4)

p3 <- DimPlot(mac,label = F,reduction = "umap2",cols = color,raster = F,alpha = 0.1)+theme_void()+NoLegend()
p3
ggsave(p3,filename = paste0(path,"umap_mac3.png"),width = 4,height =4)
ggsave(p3,filename = paste0(path,"umap_mac3.pdf"),width = 4,height =4)

rm(mar)
rm(mac)
rm(mou)
### all cell 

dat <- qread("")
Idents(dat) <- dat$anno
p1 <- DimPlot(dat,label = T,reduction = "umap",cols = color,raster = F,alpha = 0.2)+theme_void()+NoLegend()

ggsave(p1,filename = paste0(path,"umap_all_cell.png"),width = 4,height =4)
ggsave(p1,filename = paste0(path,"umap_all_cell.pdf"),width = 4,height =4)

p2 <- DimPlot(dat,label = F,reduction = "umap",cols = color,raster = F,alpha = 0.2)+theme_void()+NoLegend()
ggsave(p2,filename = paste0(path,"umap_all_cell2.png"),width = 4,height =4)
ggsave(p2,filename = paste0(path,"umap_all_cell2.pdf"),width = 4,height =4)

p3 <- DimPlot(dat,label = F,reduction = "umap",cols = color,raster = F,alpha = 0.1)+theme_void()+NoLegend()
ggsave(p3,filename = paste0(path,"umap_all_cell3.png"),width = 4,height =4)
ggsave(p3,filename = paste0(path,"umap_all_cell3.pdf"),width = 4,height =4)

color3  <- c("#e20612","#F57F17","#00b0eb")
names(color3) <- c("macaque","marmoset","mouse")

Idents(dat) <- dat$specie
p4 <- DimPlot(dat,reduction = "umap",cols = color3,raster = F,alpha = 0.1)+theme_void()+NoLegend()
p4
ggsave(p4,filename = paste0(path,"umap_all_cell4_spe.png"),width = 4,height =4)
ggsave(p4,filename = paste0(path,"umap_all_cell4_spe.pdf"),width = 4,height =4)












