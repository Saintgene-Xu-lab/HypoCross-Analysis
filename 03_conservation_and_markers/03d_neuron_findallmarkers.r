
if(T){
  library(Seurat)
  library(tidyverse)
  library(reshape2)
  library(stringr)
  library(clustree)
  library(ape)
  library(ggplot2)
  library(tidyverse)
  library(reshape2)
  library(qs)
  library(ggtree)
  library(ape)
  library(ggplot2)
  gc()
  library(COSG)
}

#########
dat <- qread("")
dat$anno_neu <- factor(dat$anno_neu,levels = c(paste0("EN",1:30),paste0("IN",1:31),"CHAT"))
Idents(dat) <- dat$anno_neu

DefaultAssay(dat) <- "RNA"
dat <- NormalizeData(dat)

################################################### Deg with in specie  -------------------------------------------------
deg <- list()
specie <- c("macaque","marmoset","mouse")
Idents(dat) <- dat$anno_neu
DefaultAssay(dat) <- "RNA"
dat <- NormalizeData(dat)
dat2 <- SplitObject(object = dat,split.by = "specie")
for (i in 1:length(specie)){
  dat2[[specie[i]]] <- NormalizeData(dat2[[specie[i]]])
  a <- FindAllMarkers(object = dat2[[specie[i]]],logfc.threshold = 0.25,min.pct = 0.25,only.pos = T)
  a <- a[a$avg_log2FC > 0.5 & a$p_val_adj < 0.05,]
  deg[[specie[i]]] <- a
}

deg_beifen <- deg
#deg <- deg_beifen

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

path <- ""
writexl::write_xlsx(mac_dat_deg2,path = paste0(path,"mac_neu_deg_fc0.5_padj0.05.xlsx"))
writexl::write_xlsx(mar_dat_deg2,path = paste0(path,"mar_neu_deg_fc0.5_padj0.05.xlsx"))
writexl::write_xlsx(mou_dat_deg2,path = paste0(path,"mou_neu_deg_fc0.5_padj0.05.xlsx"))
saveRDS(deg_beifen,file = paste0(path,"deg_in_single_specie.rds"))







