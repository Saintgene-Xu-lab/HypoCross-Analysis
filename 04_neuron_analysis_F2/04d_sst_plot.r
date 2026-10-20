

library(Seurat)
library(qs)
library(ggplot2)

dat <- qread("")
Idents(dat) <- dat$anno7
DefaultAssay(dat) <- "RNA"

p <- FeaturePlot(dat,features = "Sst",raster = F,pt.size = 0.1)+
  scale_color_gradientn(colours = colorRampPalette(c("#EEEEEE","#EF6C00"))(100))
p
dpath <- " "
ggsave(p,filename = paste0(path,"Sst.pdf"),width = 10,height = 10)
ggsave(p,filename = paste0(path,"Sst.png"),width = 10,height = 10)

#### 
p <- FeaturePlot(dat,features = "Sst",raster = F,alpha = 0.1,pt.size = 0.1,order = T,
                 cols = colorRampPalette(c("#EEEEEE","#FFF3E0","#FFCC80","#EF6C00"))(100))+theme_nothing()
path <- " "
ggsave(p,filename = paste0(path,"Sst2.pdf"),width = 10,height = 10)
ggsave(p,filename = paste0(path,"Sst2.png"),width = 10,height = 10)

p <- FeaturePlot(dat,features = "Sst",raster = F,alpha = 0.1,pt.size = 0.1,order = T,
                 cols = colorRampPalette(c("#EEEEEE","#FFF3E0","#FFCC80","#EF6C00"))(100))+theme_nothing()
path <- " "
ggsave(p,filename = paste0(path,"Sst3.pdf"),width = 6,height = 6)
ggsave(p,filename = paste0(path,"Sst3.png"),width = 6,height = 6)

p <- FeaturePlot(dat,features = "Sst",raster = F,alpha = 0.1,pt.size = 0.1,order = T)+
  scale_color_gradientn(colours = colorRampPalette(c("#EEEEEE","#EF6C00"))(100))+theme_nothing()
path <- " "
ggsave(p,filename = paste0(path,"Sst4.pdf"),width = 6,height = 6)
ggsave(p,filename = paste0(path,"Sst4.png"),width = 6,height = 6)

p <- FeaturePlot(dat,features = "Sst",
                 cols = colorRampPalette(c("#EEEEEE","#FFF3E0","#FFCC80","#EF6C00"))(100))
path <- " "
ggsave(p,filename = paste0(path,"Sst_legend.pdf"),width = 3.5,height = 3.5)


