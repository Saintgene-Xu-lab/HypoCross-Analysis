if(T){
  library(ggplot2)
  library(ggpubr)
  library(Seurat)
  library(reshape2)
  library(cowplot)
  library(qs)
  library(stringr)
}

deg <- readRDS("deg_in_single_specie.rds")
deg$macaque$cluster <- as.character(deg$macaque$cluster)
deg$mouse$cluster <- as.character(deg$mouse$cluster)
deg$marmoset$cluster <- as.character(deg$marmoset$cluster)

get_sig_order <- function(deg){
  deg2 <- deg[deg$p_val_adj < 0.05 & deg$avg_log2FC > 0.5,]
  deg2$cluster <- factor(deg2$cluster,levels = c(paste0("EN",1:30),paste0("IN",1:31),"CHAT"))
  deg2 <- deg2[order(deg2$cluster),]
  return(deg2)
}

mou_dat_deg2 <- get_sig_order(deg = deg$mouse)
mar_dat_deg2 <- get_sig_order(deg = deg$marmoset)
mac_dat_deg2 <- get_sig_order(deg = deg$macaque)

deg <- list(mac = mac_dat_deg2,mar = mar_dat_deg2, mou = mou_dat_deg2)
cells <- c(paste0("EN",1:30),paste0("IN",1:31),"CHAT")
###### plot heatmap -------

cell_common2 <- list()
for (i in 1:length(cells)) {
  a <- intersect(deg$mou[deg$mou$cluster %in% cells[i],]$gene,deg$mar[deg$mar$cluster %in% cells[i],]$gene)
  a <- intersect(a,deg$mac[deg$mac$cluster %in% cells[i],]$gene)
  cell_common2[[cells[i]]] <- a
  
}

primate_up2 <- list()
for (i in 1:length(cells)) {
  a <- intersect(deg$mac[deg$mac$cluster %in% cells[i],]$gene,deg$mar[deg$mar$cluster %in% cells[i],]$gene)
  a <- setdiff(a,deg$mou[deg$mou$cluster %in% cells[i],]$gene)

  primate_up2[[cells[i]]] <- a
  
}

mac_up2 <- list()
for (i in 1:length(cells)) {
  a <- unique(c(deg$mou[deg$mou$cluster %in% cells[i],]$gene,deg$mar[deg$mar$cluster %in% cells[i],]$gene))
  a <- setdiff(deg$mac[deg$mac$cluster %in% cells[i],]$gene,a)

  mac_up2[[cells[i]]] <- a
}

mar_up2 <- list()
for (i in 1:length(cells)) {
  a <- unique(c(deg$mou[deg$mou$cluster %in% cells[i],]$gene,deg$mac[deg$mac$cluster %in% cells[i],]$gene))
  a <- setdiff(deg$mar[deg$mar$cluster %in% cells[i],]$gene,a)

  mar_up2[[cells[i]]]  <- a
}

mou_up2 <- list()
for (i in 1:length(cells)) {
  a <- unique(c(deg$mar[deg$mar$cluster %in% cells[i],]$gene,deg$mac[deg$mac$cluster %in% cells[i],]$gene))
  a <- setdiff(deg$mou[deg$mou$cluster %in% cells[i],]$gene,a)

  mou_up2[[cells[i]]]  <- a
}

common2 <- unique(unlist(cell_common2))
primate_up2 <- unique(unlist(primate_up2))
mac_up2 <- unique(unlist(mac_up2))
mar_up2 <- unique(unlist(mar_up2))
mou_up2 <- unique(unlist(mou_up2))

dat_plot <- data.frame(group = c("common","primate","mac","mar","mou"),num = c(length(common2),length(primate_up2),length(mac_up2),length(mar_up2),length(mou_up2)) )

library(ggpubr)
p <- ggdotchart(dat_plot, x="group", y="num", color = "group",dot.size =7,
           add = "segments", ggtheme = theme_pubr())
path <- ""
ggsave(p,filename = paste0(path,"dot_num.pdf"),width = 3,height = 4)

