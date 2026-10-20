if(T){
  library(Seurat)
  library(ggplot2)
  library(tidyverse)
  library(reshape2)
  library(qs)
  library(RColorBrewer)
  library(cowplot)
  library(viridis)
  library(COSG)
  display.brewer.all()
  library(ggdensity)
	library(ggpubr)
}

load("use_color8_20241113_update.Rdata")

mar <- qread("")
Idents(mar) <- mar$chip2
mar_t008 <- subset(mar,idents = "T008")
Idents(mar_t008) <- mar_t008$Assign_Simple2
names(use_color_update2)[22] <- "SCN"
p2 <- DimPlot(mar_t008,cols = use_color_update2,label = T)+theme_nothing()
p2
ggsave(p2,filename = paste0(path,"umap.pdf"),width = 3,height = 3)
Idents(mar_t008) <- mar_t008$BANKSY_clu5
p1 <- DimPlot(mar_t008,cells.highlight = WhichCells(mar_t008,idents = "MarMG"),cols = "#EEEEEE",sizes.highlight = 0.5,cols.highlight = "#7f5e9d")+theme_nothing()
path <- ""
ggsave(p1,filename = paste0(path,"MarMG_DOMAIN2.pdf"),width = 3,height = 3)

MarMG <- mar_t008@meta.data
MarMG <- MarMG[MarMG$BANKSY_clu5 %in% "MarMG",]
dat <- data.frame(table(MarMG$laber_transfer_1109))
dat <- dat[order(dat$Freq,decreasing = T),]

dat$Var1
dat$Var2 <- paste0("C",1:length(dat$Var1))
dat$Var2 <- factor(dat$Var2,levels = rev(dat$Var2))

p <- ggplot(dat,aes(x = Freq,y = Var2,fill = Freq,color = Freq))+
  geom_bar(stat = "identity",size = 0.1)+
  scale_color_gradientn(colours = colorRampPalette(brewer.pal(9,'RdPu')[2:9])(50))+
  scale_fill_gradientn(colours = colorRampPalette(brewer.pal(9,'RdPu')[2:9])(50))+
  theme_bw()
p
path <- ""
ggsave(p,filename = paste0(path,"mg_number_inMG_domain.pdf"),width = 4,height = 15)

