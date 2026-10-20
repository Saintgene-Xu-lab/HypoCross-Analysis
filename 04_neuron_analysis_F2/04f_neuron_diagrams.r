
if(T){
  library(Seurat)
  library(qs)
  library(viridis)
  library(ComplexHeatmap)
  library(ggplot2)
  library(RColorBrewer)
}

neu <- qread("neu_counts_rm_female_20241029.qs")
DimPlot(neu,label = T)
cells <- levels(neu$anno6)

color3  <- c("#e20612","#F57F17","#00b0eb")
names(color3) <- c("macaque","marmoset","mouse")

###############1 Cell abundance 
meta <- qread("./meta_rm_female_20241030.qs")
meta_neu <- meta[meta$anno %in% c("NEU"),]
color3  <- c("#e20612","#F57F17","#00b0eb")
names(color3) <- c("macaque","marmoset","mouse")

col <- c(colorRampPalette(brewer.pal(12,'Paired'))(32),colorRampPalette(brewer.pal(8,'Set1'))(30))
names(col) <- c(paste0("EN",1:30),paste0("IN",1:31),"CHAT")

#meta <- all_cell@meta.data
meta2 <- meta_neu[!(meta_neu$sample_area %in% c("Arc","ArcC","LH","MM","POA","PVN","SCN","VMH","ZI")),]

set.seed(12580)
mou_cell <- rownames(meta2[meta2$specie %in% "mouse",])

##
CR_species <- function(meta,cell,species){
  Cellratio <- prop.table(table(meta[,cell],meta[,species]), margin = 2)
  Cellratio <- as.data.frame(Cellratio)
  Cellratio$Var1 <- factor(Cellratio$Var1,levels = c(paste0("EN",1:30),paste0("IN",1:31),"CHAT"))
  p <- ggplot(Cellratio) + 
    geom_bar(aes(x =Var1, y= Freq, fill = Var2),stat = "identity",width = 0.7,linewidth = 0.5,colour = 'white',position = 'fill')+ # position = 'fill' 
    scale_fill_manual(values =  c("#e20612","#F57F17","#00b0eb"))+
    labs(x='Sample',y = 'Ratio')+
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90,hjust = 1))
  return(p)
} # 
p1 <- CR_species(meta = meta2,cell ="anno_neu" ,species = "specie") # 3 12
p1

CR_species3 <- function(meta,cell,species){
  Cellratio <- prop.table(table(meta[,cell],meta[,species]), margin = 2)
  Cellratio <- as.data.frame(Cellratio)
  Cellratio$Var1 <- factor(Cellratio$Var1,levels = c(paste0("EN",1:30),paste0("IN",1:31),"CHAT"))
  p <- ggplot(Cellratio) + 
    geom_bar(aes(x =Var1, y= Freq, fill = Var2),stat = "identity",width = 0.7,linewidth = 0.5,colour = 'white',position = 'stack')+ # position = 'fill' 
    scale_fill_manual(values =  c("#e20612","#F57F17","#00b0eb"))+
    labs(x='Sample',y = 'Ratio')+
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90,hjust = 1))
  return(p)
}
p3 <- CR_species3(meta = meta2,cell ="anno_neu" ,species = "specie")

p <- p1/p3
ggsave(p,filename = paste0(out_path,"cell_ratio.pdf"),width = 10,height = 4)

library(ggpubr)
a <- data.frame(table(meta_neu$anno_neu))
a$Var1 <- factor(a$Var1,levels = c(paste0("EN",1:30),paste0("IN",1:31),"CHAT"))
a$lgFreq <- log1p(a$Freq)


p_dot  <- ggplot(a, aes(x = Var1, y = lgFreq)) +
  geom_segment( aes(x = Var1, xend = Var1, y = 0, yend = lgFreq),color = "grey40")+ 
  geom_point(aes(color = Var1),size = 3) + 
  theme_bw()+
  scale_color_manual(values = col)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),legend.position = "none")
p_dot

p_dot2  <- ggplot(a, aes(x = Var1, y = Freq)) +
  geom_segment( aes(x = Var1, xend = Var1, y = 0, yend = Freq),color = "grey40")+ 
  geom_point(aes(color = Var1),size = 3) + 
  theme_bw()+
  scale_color_manual(values = col)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),legend.position = "none")
p_dot2

p <- p1/p_dot2# 4 10
p
ggsave(p,filename = paste0(out_path,"cell_ratio3.pdf"),width = 20,height = 4)

p_bar  <- ggplot(a, aes(x = Var1, y = Freq)) +
  geom_segment( aes(x = Var1, xend = Var1, y = 0, yend = Freq),color = "grey40")+ 
  geom_bar(aes(fill = Var1),stat = "identity") + 
  theme_bw()+
  scale_fill_manual(values = col)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),legend.position = "none")
p_bar 

p_bar2  <- ggplot(a, aes(x = Var1, y = lgFreq)) +
  geom_segment( aes(x = Var1, xend = Var1, y = 0, yend = lgFreq),color = "grey40")+ #
  geom_bar(aes(fill = Var1),stat = "identity") + 
  theme_bw()+
  scale_fill_manual(values = col)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),legend.position = "none")
p_bar2 

################ one cluster vs other species
neu_sample <- PrepSCTFindMarkers(neu_sample)

mouse_vs_other <- list()
marmoset_vs_other <- list()
macaque_vs_other <- list()

cells <- c(paste0("EN",1:30),paste0("IN",1:31),"CHAT")
specie <- c("mouse","macaque","marmoset")

neu_sample$anno_neu <- as.character(neu_sample$anno6)

for (i in 1:length(cells)){
  neu_sample$anno7 <- neu_sample$specie
  neu_sample@meta.data[neu_sample$specie %in% c("macaque","marmoset"),]$anno7 <- "other"
  Idents(neu_sample) <- paste0(neu_sample$anno_neu,"_",neu_sample$anno7)
  ### get deg 
  mouse_vs_other[[cells[i]]] <- FindMarkers(object = neu_sample,
                                            ident.1 = paste0(cells[i],"_mouse"),
                                            ident.2 = paste0(cells[i],"_other"),logfc.threshold = 0.25,min.pct = 0.25,only.pos = F)
}
for (i in 1:length(cells)){
  neu_sample$anno7 <- neu_sample$specie
  neu_sample@meta.data[neu_sample$specie %in% c("macaque","mouse"),]$anno7 <- "other"
  Idents(neu_sample) <- paste0(neu_sample$anno_neu,"_",neu_sample$anno7)
  ### get deg 
  marmoset_vs_other[[cells[i]]] <- FindMarkers(object = neu_sample,
                                               ident.1 = paste0(cells[i],"_marmoset"),
                                               ident.2 = paste0(cells[i],"_other"),logfc.threshold = 0.25,min.pct = 0.25,only.pos = F)
}
for (i in 1:length(cells)){
  neu_sample$anno7 <- neu_sample$specie
  neu_sample@meta.data[neu_sample$specie %in% c("marmoset","mouse"),]$anno7 <- "other"
  Idents(neu_sample) <- paste0(neu_sample$anno_neu,"_",neu_sample$anno7)
  ### get deg 
  macaque_vs_other[[cells[i]]] <- FindMarkers(object = neu_sample,
                                              ident.1 = paste0(cells[i],"_macaque"),
                                              ident.2 = paste0(cells[i],"_other"),logfc.threshold = 0.25,min.pct = 0.25,only.pos = F)
}

test <- lapply(mouse_vs_other, function(x){
  x <- x[x$avg_log2FC > 0.5 & x$p_val_adj < 0.05,]
})
unlist(lapply(test, function(x){x <- nrow(x)}))

test2 <- lapply(marmoset_vs_other, function(x){
  x <- x[x$avg_log2FC > 0.5 & x$p_val_adj < 0.05,]
})
test3 <- lapply(macaque_vs_other, function(x){
  x <- x[x$avg_log2FC > 0.5 & x$p_val_adj < 0.05,]
})
unlist(lapply(test, function(x){x <- nrow(x)}))

library(reshape2)
deg_num <- data.frame(mouse= unlist(lapply(test, function(x){x <- nrow(x)})),
                      marmoset =unlist(lapply(test2, function(x){x <- nrow(x)})),
                      macaque = unlist(lapply(test3, function(x){x <- nrow(x)})))
deg_num$cell <- rownames(deg_num)
deg_num2 <- reshape2::melt(deg_num)
deg_num2 <- deg_num2[deg_num2$cell %in% cells[1:62],]

color3  <- c("#e20612","#F57F17","#00b0eb")
names(color3) <- c("macaque","marmoset","mouse")
deg_num2$cell <- factor(deg_num2$cell,levels = cells[1:62])
deg_num2$radio <- deg_num2$value / length(rownames(neu_sample2))

deg_num2 %>% ggplot(aes(x =cell,y =value,color = variable))+
  geom_point(size = 3)+
  scale_color_manual(values = color3)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))-> p3 # 3 30

deg_num2 %>% ggplot(aes(x =cell,y =radio,color = variable))+
  geom_point(size = 5)+
  scale_color_manual(values = color3)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1)) # 3 30

deg_num2 %>% ggplot(aes(x =cell,y =  log1p(value),color = variable))+
  geom_point(size = 5)+
  scale_color_manual(values = color3)+
  theme_bw()+
  scale_y_continuous(limits = c(0,10))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1)) # 3 30

p <-  p1/p_dot2/p3# 4 10
p
ggsave(p,filename = paste0(out_path,"cell_ratio4.pdf"),width = 20,height = 6)

