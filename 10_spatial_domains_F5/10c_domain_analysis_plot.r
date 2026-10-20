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
}

###################################### Aim ##############################################
### 1 get umap.rpca of spatial domain 
### 2 get ratio of spatial domain 
### 3 get location of spatial domain 


###################################### loading data #####################################
int_bin100 <- qread("int.qs")
mou <- qread("mou3.qs")
mar <- qread("mar3.qs")
mac <- qread("mac3.qs")

if(T){
  color3  <- c("#e20612","#F57F17","#00b0eb")
  names(color3) <-  c("macaque", "marmoset", "mouse")
 
  colorlist <- c("#a4cde1","#67a4cc","#277fb8","#549da3","#96cb8f","#8bc96d","#4dae47","#5c9e43","#b79973","#f38989","#ec5051","#e32427","#ef6a45","#f9b769","#f9a341","#f48521","#ee8e46","#d4a6a8","#af93c4","#8660a8","#815e99","#c6b598","#f6f28f","#d4a55b","#b05a28")
  colorlist <- colorRampPalette(colorlist)(26)
  names(colorlist) <- paste0("D",1:22)
}

path <- ""
############################################### 1 plot umap.rpca #############################################################
if(T){
  Idents(int_bin100) <- int_bin100$domain_motif
  p <- DimPlot(int_bin100,label = F,reduction = "umap.rpca",raster = F,cols = colorlist,alpha = 0.3)+theme_nothing()
  p
  ggsave(p,filename = paste0(path,"bin100_spatial_domain.pdf"),width = 4,height = 4)
  
  p <- DimPlot(int_bin100,label = F,reduction = "umap.rpca",raster = F,cols = colorlist,alpha = 0.1)+theme_nothing()
  p
  ggsave(p,filename = paste0(path,"bin100_spatial_domain2.pdf"),width = 4,height = 4)
  
  
  p <- DimPlot(int_bin100,label = T,reduction = "umap.rpca",raster = T,cols = colorlist,alpha = 0.3)+theme_nothing()
  p
  ggsave(p,filename = paste0(path,"bin100_spatial_domain_label.pdf"),width = 4,height = 4)
  
  p <- DimPlot(int_bin100,label = T,reduction = "umap.rpca",raster = T,cols = colorlist,alpha = 0.3)
  p
  ggsave(p,filename = paste0(path,"bin100_spatial_domain_label2.pdf"),width = 4,height = 4)
  
  
  p <- DimPlot(int_bin100,label = F,group.by = "specie",reduction = "umap.rpca",raster = F,cols =color3,alpha = 0.5)+theme_nothing()
  p
  ggsave(p,filename = paste0(path,"bin100_spatial_domain_split_spe.pdf"),width = 4,height = 4)
  
  meta <- int_bin100@meta.data
  Idents(int_bin100) <- int_bin100$specie
  p <- DimPlot(int_bin100,label = F,split.by = "specie",reduction = "umap.rpca",raster = F,cols =color3,alpha = 0.5)+theme_nothing()
  p
  ggsave(p,filename = paste0(path,"bin100_spatial_domain_split_spe2.pdf"),width = 12.5,height = 4)
}

############################################### 2 plot ratio #############################################################
meta <- int_bin100@meta.data

CR_species <- function(meta,cell,species){
  Cellratio <- prop.table(table(meta[,cell],meta[,species]), margin = 2)
  Cellratio <- as.data.frame(Cellratio)
  Cellratio$Var1 <- factor(Cellratio$Var1,levels = paste0("D",1:22))
  p <- ggplot(Cellratio) + 
    geom_bar(aes(x =Var1, y= Freq, fill = Var2),stat = "identity",width = 0.7,linewidth = 0.5,colour = 'white',position = 'fill')+ # position = 'fill' 
    scale_fill_manual(values =  c("#e20612","#F57F17","#00b0eb"))+
    labs(x='Sample',y = 'Ratio')+
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90,hjust = 1))
  return(p)
} 

p <- CR_species(meta = meta,cell ="domain_motif" ,species = "specie")
p
ggsave(p,filename = paste0(path,"domian_ratio.pdf"),width = 10,height = 4)

############################################### 3 domain location #########################################################

id <- levels(int_bin100$domain_motif)
plot_domain_location <- function(mou_dat,mar_dat,mac_dat,domain,id,color){
  
  library(stringr)
  p_mou <- DimPlot(mou_dat,reduction = "umap2",cells.highlight = str_remove(WhichCells(domain,idents = id),"mou_"),
                   raster = F,sizes.highlight = 0.1,cols ="#EEEEEE",cols.highlight = color)+NoLegend()
  p_mar <- DimPlot(mar_dat,reduction = "umap2",cells.highlight = str_remove(WhichCells(domain,idents = id),"mar_"),
                   raster = F,sizes.highlight = 0.1,cols ="#EEEEEE",cols.highlight = color)+NoLegend()
  p_mac <- DimPlot(mac_dat,reduction = "umap2",cells.highlight = str_remove(WhichCells(domain,idents = id),"mac_"),
                   raster = F,sizes.highlight = 0.1,cols ="#EEEEEE",cols.highlight = color)+NoLegend()
  library(patchwork)
  p <- p_mou + p_mar + p_mac
  return(p)
  
}


Idents(int_bin100) <- int_bin100$domain_motif
path <- ""
for (i in 1:length(id)) {
  print(paste0("Doing",id[i]))
  p <- plot_domain_location(mou_dat = mou,mar_dat = mar,mac_dat = mac,id = id[i],domain = int_bin100,color = colorlist[i])
  ggsave(p,filename = paste0(path,id[i],".png"),width = 45,height = 15,limitsize = F)
  ggsave(p,filename = paste0(path,id[i],".pdf"),width = 45,height = 15,limitsize = F)
  print(paste0("Done",id[i]))
}


############################################# 4 plot domain sperate chip ################################3
####### mouse 

nucle_color <- readRDS("use_color.rds")
plot_location_domain3 <- function(dat,domain,spe,col,domian_id,path,nucle_color){
  
  Idents(domain) <- domain$domain_motif
  cells <- WhichCells(domain,idents = domian_id)
  chip <- names(table(dat$chip))
  Idents(dat) <- dat$chip
  
  for (i in 1:length(chip)) {
    ncell <- intersect(str_remove(cells,paste0(spe,"_")),WhichCells(dat,idents = chip[i]))
    
    if (length(ncell) > 10) {
      p <- DimPlot(dat,cells.highlight = ncell,cells = WhichCells(dat,idents = chip[i]),pt.size = 2,sizes.highlight = 2,
                   cols = "#eaeaea",raster = F,cols.highlight = col[domian_id])+theme_nothing()
      ggsave(p,filename = paste0(path,domian_id,"_",chip[i],"_",spe,".pdf"),width = 3.3,height = 3)
      ggsave(p,filename = paste0(path,domian_id,"_",chip[i],"_",spe,".png"),width = 3.3,height = 3)
      
      p <- DimPlot(dat,cells = WhichCells(dat,idents = chip[i]),pt.size = 2,label = T,group.by = "Assign_Simple2",
                   cols = nucle_color,raster = F)+theme_nothing()
      ggsave(p,filename = paste0(path,chip[i],"_nucle_",spe,".pdf"),width = 4.5,height = 4)
    }else{
      print(paste0("Skip : ",chip[i],"_",spe))
    }
    
    
  }
  
}
path <- "/mou/"
id <- names(table(int_bin100$domain_motif))
for (i in 1:length(id)) {
  print(paste0("Doing",id[i]))
  plot_location_domain3(dat = mou,domain = int_bin100,spe = "mou",col = colorlist,domian_id = id[i],path = path,nucle_color = nucle_color)
  print(paste0("Done",id[i]))
}

###### marmoset
nucle_color <- readRDS("use_color.rds")
plot_location_domain3 <- function(dat,domain,spe,col,domian_id,path,nucle_color){
  
  Idents(domain) <- domain$domain_motif
  cells <- WhichCells(domain,idents = domian_id)
  chip <- names(table(dat$chip))
  Idents(dat) <- dat$chip
  for (i in 1:length(chip)) {
    ncell <- intersect(str_remove(cells,paste0(spe,"_")),WhichCells(dat,idents = chip[i]))
    
    if (length(ncell) > 10) {
      p <- DimPlot(dat,cells.highlight = ncell,cells = WhichCells(dat,idents = chip[i]),pt.size = 2,sizes.highlight = 2,
                   cols = "#eaeaea",raster = F,cols.highlight = col[domian_id])+theme_nothing()
      ggsave(p,filename = paste0(path,domian_id,"_",chip[i],"_",spe,".pdf"),width = 3,height = 3)
      ggsave(p,filename = paste0(path,domian_id,"_",chip[i],"_",spe,".png"),width = 3,height = 3)
      
      p <- DimPlot(dat,cells = WhichCells(dat,idents = chip[i]),pt.size = 2,label = T,group.by = "Assign_Simple2",
                   cols = nucle_color,raster = F)+theme_nothing()
      ggsave(p,filename = paste0(path,chip[i],"_nucle_",spe,".pdf"),width = 4.5,height = 4)
    }else{
      print(paste0("Skip : ",chip[i],"_",spe))
    }
    
    
  }
  
}
path <- "/mar/"
id <- names(table(int_bin100$domain_motif))
for (i in 1:length(id)) {
  print(paste0("Doing",id[i]))
  plot_location_domain3(dat = mar,domain = int_bin100,spe = "mar",col = colorlist,domian_id = id[i],path = path,nucle_color = nucle_color)
  print(paste0("Done",id[i]))
}

###### macaque

plot_location_domain3 <- function(dat,domain,spe,col,domian_id,path,nucle_color){
  
  Idents(domain) <- domain$domain_motif
  cells <- WhichCells(domain,idents = domian_id)
  chip <- names(table(dat$chip))
  Idents(dat) <- dat$chip
  
  for (i in 1:length(chip)) {
    ncell <- intersect(str_remove(cells,paste0(spe,"_")),WhichCells(dat,idents = chip[i]))
    
    if (length(ncell) > 10) {
      p <- DimPlot(dat,cells.highlight = ncell,cells = WhichCells(dat,idents = chip[i]),pt.size = 2,sizes.highlight = 2,
                   cols = "#eaeaea",raster = F,cols.highlight = col[domian_id])+theme_nothing()
      ggsave(p,filename = paste0(path,domian_id,"_",chip[i],"_",spe,".pdf"),width = 3,height = 3)
      ggsave(p,filename = paste0(path,domian_id,"_",chip[i],"_",spe,".png"),width = 3,height = 3)
      
      p <- DimPlot(dat,cells = WhichCells(dat,idents = chip[i]),pt.size = 2,label = T,group.by = "Assign_Simple2",
                   cols = nucle_color,raster = F)+theme_nothing()
      ggsave(p,filename = paste0(path,chip[i],"_nucle_",spe,".pdf"),width = 4.5,height = 4)
    }else{
      print(paste0("Skip : ",chip[i],"_",spe))
    }
    
    
  }
  
}

path <- "/mac/"
id <- names(table(int_bin100$domain_motif))
for (i in 1:length(id)) {
  print(paste0("Doing",id[i]))
  plot_location_domain3(dat = mac,domain = int_bin100,spe = "mac",col = colorlist,domian_id = id[i],path = path,nucle_color = nucle_color)
  print(paste0("Done",id[i]))
}


path <- "1_routin_plot/"
colorlist2 <- c(colorlist,"#EEEEEE")
names(colorlist2)[27] <- "HY"
Idents(mou) <- mou$domian
p <- DimPlot(mou,reduction = "umap2",cols = colorlist2,raster = F,pt.size = 0.5)
ggsave(p,filename = paste0(path,"mou_domain_in_st.png"),width = 20,height = 15)
Idents(mar) <- mar$domian
p <- DimPlot(mar,reduction = "umap2",cols = colorlist2,raster = F,pt.size = 0.5,
             cells = rownames(mar@meta.data[mar$domian %in% paste0("D",1:27),]))
ggsave(p,filename = paste0(path,"mar_domain_in_st.png"),width = 20,height = 15)
Idents(mac) <- mac$domian
p <- DimPlot(mac,reduction = "umap2",cols = colorlist2,raster = F,pt.size = 0.5,
             cells = rownames(mac@meta.data[mac$domian %in% paste0("D",1:27),]))
ggsave(p,filename = paste0(path,"mac_domain_in_st.png"),width = 20,height = 15)


################################################# 4 plot_domian compare with nucle annotation  
mou <- qread("")
mar <- qread("")
mac <- qread("")

int_bin100 <- qread("int_20251102.qs")

mou_meta <- mou@meta.data[intersect(str_remove(colnames(int_bin100)[str_detect(colnames(int_bin100),"mou")],"mou_"),colnames(mou)),]
mac_meta <- mac@meta.data[intersect(str_remove(colnames(int_bin100)[str_detect(colnames(int_bin100),"mac")],"mac_"),colnames(mac)),]
mar_meta <- mar@meta.data[intersect(str_remove(colnames(int_bin100)[str_detect(colnames(int_bin100),"mar")],"mar_"),colnames(mar)),]

get_nucle_and_domain_compare <- function(dat,spe){
  #meta <- dat@meta.data
  #meta <- meta[!str_detect(meta$domian,"HY"),]
  
  a <- t(prop.table(x = table(as.character(dat$domian_modif2),dat$Assign_Simple2),margin = 1))
  ranking_mat <- function(data){
    library(zoo)
    window <- 1
    step <- 1
    dat_ts <- data.frame(data)
    dat2_ts <- dat_ts[order(apply(t(rollapply(t(dat_ts), width=window, by=step, FUN=mean)), 1, which.max)), ]
    return(dat2_ts)
  }
  a2 <- melt(a)
  a2 <- spread(a2,key = Var2,value = value)
  head(a2)
  rownames(a2) <- a2$Var1
  a2 <- a2[,-1]
  a2 <- ranking_mat(a2)
  return(a2)
}

mou_com <- get_nucle_and_domain_compare(dat = mou_meta,spe = "mouse")
mar_com <- get_nucle_and_domain_compare(dat = mar_meta,spe = "marmoset")
mac_com <- get_nucle_and_domain_compare(dat = mac_meta,spe = "macaque")

nucle <- unique(c(rownames(mou_com),rownames(mar_com),rownames(mac_com)))
nucle <- c("POA","AHN","SO","LHA","fiber","ZI","PeVH","SCN","PVH","ARH","VMH","DMH","TU","CVO","PH","PM","MBO","SUM","HY")
nucle2 <- nucle[-c(5)]

pheatmap::pheatmap(mou_com[nucle2,],cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("#FAFAFA",rev(viridis_pal(option = "G")(10))[1:5]))(100),breaks = seq(from = 0, to = 0.8, length.out = 100),
                   border_color = "white",angle_col = 90,filename = paste0(path,"mou_domain_nucle_hy.pdf"))

pheatmap::pheatmap(mar_cor[intersect(nucle2,rownames(mar_cor)),],cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("#FAFAFA",rev(viridis_pal(option = "G")(10))[1:5]))(100),breaks = seq(from = 0, to = 0.8, length.out = 100),
                   border_color = "white",angle_col = 90,filename = paste0(path,"mar_domain_nucle_hy.pdf"))

pheatmap::pheatmap(mac_cor[intersect(nucle2,rownames(mac_cor)),],cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("#FAFAFA",rev(viridis_pal(option = "G")(10))[1:5]))(100),breaks = seq(from = 0, to = 0.8, length.out = 100),
                   border_color = "white",angle_col = 90,filename = paste0(path,"mac_domain_nucle_hy.pdf"))







plot_nucle_and_domain_compare_mat <- function(dat,spe){
  #meta <- dat@meta.data
  #meta <- meta[!str_detect(meta$domian,"HY"),]
  
  a <- t(prop.table(x = table(as.character(dat$domian_modif2),dat$Assign_Simple2),margin = 1))
  ranking_mat <- function(data){
    library(zoo)
    window <- 1
    step <- 1
    dat_ts <- data.frame(data)
    dat2_ts <- dat_ts[order(apply(t(rollapply(t(dat_ts), width=window, by=step, FUN=mean)), 1, which.max)), ]
    return(dat2_ts)
  }
  a2 <- melt(a)
  a2 <- spread(a2,key = Var2,value = value)
  head(a2)
  rownames(a2) <- a2$Var1
  a2 <- a2[,-1]
  a2 <- a2[,paste0("D",1:22)]
  a2 <- ranking_mat(a2)
  return(a2)
}
mou_cor <- plot_nucle_and_domain_compare_mat(dat = mou_meta,spe = "mouse")
mar_cor <- plot_nucle_and_domain_compare_mat(dat = mar_meta,spe = "marmoset")
mac_cor <- plot_nucle_and_domain_compare_mat(dat = mac_meta,spe = "macaque")

nucle <- unique(c(rownames(mou_cor),rownames(mar_cor),rownames(mac_cor)))
nucle <- c("POA","AHN","SO","LHA","fiber","ZI","PeVH","SCN","PVH","ARH","VMH","DMH","TU","CVO","PH","PM","MBO","SUM","HY")
nucle2 <- nucle[-c(5,19)]
nucle3 <- intersect(intersect(intersect(nucle2,rownames(mou_cor)),rownames(mar_cor)),rownames(mac_cor))

pheatmap::pheatmap(mou_cor[intersect(nucle2,rownames(mou_cor)),],cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("#FAFAFA",rev(viridis_pal(option = "G")(10))[1:5]))(50),
                   border_color = "white",angle_col = 90)

pheatmap::pheatmap(mar_cor[intersect(nucle2,rownames(mar_cor)),],cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("#FAFAFA",rev(viridis_pal(option = "G")(10))[1:5]))(50),
                   border_color = "white",angle_col = 90)

pheatmap::pheatmap(mac_cor[intersect(nucle2,rownames(mac_cor)),],cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("#FAFAFA",rev(viridis_pal(option = "G")(10))[1:5]))(50),
                   border_color = "white",angle_col = 90)


mou_cor2 <- mou_cor[intersect(nucle3,rownames(mou_cor)),]
mar_cor2 <- mar_cor[intersect(nucle3,rownames(mar_cor)),]
mac_cor2 <- mac_cor[intersect(nucle3,rownames(mac_cor)),]

library(Hmisc)
library(psych)
cor.result <- psych::corr.test(mou_cor2,mar_cor2,method="pearson",adjust="holm")
result.r <- cor.result$r
result.p <- cor.result$p

pheatmap::pheatmap(cor.result$p,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("#FAFAFA",rev(viridis_pal(option = "G")(10))[1:5]))(50),
                   border_color = "white",angle_col = 90)



a <- cor(mou_cor2,mar_cor2)
a_test <- rcorr(data.frame(mou_cor2),data.frame(mar_cor2))
a1 <- cor(mou_cor2,mac_cor2)
a_test1 <- rcorr(mou_cor2,mac_cor2)
a2 <- cor(mar_cor2,mac_cor2)
a_test2 <- rcorr(mar_cor2,mac_cor2)

pheatmap::pheatmap(a,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("#FAFAFA",rev(viridis_pal(option = "G")(10))[1:5]))(50),
                   border_color = "white",angle_col = 90)

pheatmap::pheatmap(a1,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("#FAFAFA",rev(viridis_pal(option = "G")(10))[1:5]))(50),
                   border_color = "white",angle_col = 90)

pheatmap::pheatmap(a2,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("#FAFAFA",rev(viridis_pal(option = "G")(10))[1:5]))(50),
                   border_color = "white",angle_col = 90)



################################################# 5 Domian features 
library(COSG)
Idents(int_bin100) <- int_bin100$domain
int_bin100
marker_ske <- cosg(
  int_bin100,
  groups='all',
  assay="RNA",
  slot="data",
  mu=1,
  n_genes_user=200)

save(marker_ske,file = "")














