if(T){
  library(Seurat)
  library(ggplot2)
  library(tidyverse)
  library(COSG)
  library(qs)
  library(stringr)
  library(patchwork)
  library(cowplot)
}

############
whitePurple <- c("#F5F5F5" ,"#ffe3e3", "#ffc9c9", "#ffa8a8", "#ff8787", "#ff6b6b", "#c92a2a")

Nebulosa::plot_density(object = mou_t321,features = "Agrp",reduction = "umap",size = 2)+scale_color_gradientn(colours = colorRampPalette(whitePurple)(100))
FeaturePlot(mou_t321,features = c("Agrp"),reduction = "umap2",label = T)+scale_color_gradientn(colours = colorRampPalette(whitePurple)(100))

horizonExtra <- c("#000436","#021EA9", "#1632FB", "#6E34FC" ,"#C732D5" ,"#FD619D" ,"#FF9965", "#FFD32B" ,"#FFFC5A" )
Nebulosa::plot_density(object = mou_t321,features = "Agrp",reduction = "umap",size = 2)+scale_color_gradientn(colours = colorRampPalette(horizonExtra)(100))

solarExtra <- c("#3361A5","#248AF3" ,"#14B3FF", "#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D" )
Nebulosa::plot_density(object = mou_t321,features = "Agrp",reduction = "umap",size = 2)+scale_color_gradientn(colours = colorRampPalette(solarExtra)(100))

whitePurple <- c("#F5F5F5" ,"#ffe3e3", "#ffc9c9", "#ffa8a8",  "#ff6b6b", "#c92a2a")
mar <- qread("")
mac <- qread("")
mou <- qread("")
mou <- NormalizeData(mou)
mar <- NormalizeData(mar)
mac <- NormalizeData(mac)


p <- Nebulosa::plot_density(object = mou_t321,features = "Agrp",reduction = "umap",size = 2)+scale_color_gradientn(colours = colorRampPalette(whitePurple)(100))
ggsave(p,filename = paste0(save_plot,"legend",".pdf"),width = 4,height = 4)

use_genes <- c("Agrp","Pomc","Npy","Ghrh","Kiss1","Tbx3","Cited1")
save_plot <- ""

plot_st_gene2 <- function(seu,chip,gene){
  dat <- subset(seu,idents = chip)
  p <- Nebulosa::plot_density(object = dat,features = gene,reduction = "umap",size = 2)+scale_color_gradientn(colours = colorRampPalette(whitePurple)(100))+
    theme_nothing()
  return(p)
}

### mouse
Idents(mou) <- mou$chip
for (i in 1:length(use_genes)) {
  chip <- unique(mou$chip)
  for (j in 1:length(chip)) {
    result <- tryCatch({
      p <- plot_st_gene2(seu = mou,chip = chip[j],gene = use_genes[i])
      ggsave(p,filename = paste0(save_plot,"mou_",use_genes[i],"_",chip[j],".pdf"),width = 3,height = 3)
      ggsave(p,filename = paste0(save_plot,"mou_",use_genes[i],"_",chip[j],".png"),width = 4,height = 4)
    }, error = function(e) {
      print(paste0(use_genes[i],"_",chip[j]))
    })
    
    print(result)
  }
  print(paste0("Done","_",use_genes[i]))
}

### mar
Idents(mar) <- mar$chip
use_genes <- str_to_upper(use_genes)
for (i in 1:length(use_genes)) {
  chip <- unique(mar$chip)
  for (j in 1:length(chip)) {
    result <- tryCatch({
      p <- plot_st_gene2(seu = mar,chip = chip[j],gene = use_genes[i])
      ggsave(p,filename = paste0(save_plot,"mar_",use_genes[i],"_",chip[j],".pdf"),width = 4,height = 4)
      ggsave(p,filename = paste0(save_plot,"mar_",use_genes[i],"_",chip[j],".png"),width = 4,height = 4)
    }, error = function(e) {
      print(paste0(use_genes[i],"_",chip[j]))
    })
    
    print(result)
  }
  print(paste0("Done","_",use_genes[i]))
}

### mac
Idents(mac) <- mac$chip
use_genes <- str_to_upper(use_genes)
for (i in 1:length(use_genes)) {
  chip <- unique(mac$chip)
  for (j in 1:length(chip)) {
    result <- tryCatch({

      p <- plot_st_gene2(seu = mac,chip = chip[j],gene = use_genes[i])
      ggsave(p,filename = paste0(save_plot,"mac_",use_genes[i],"_",chip[j],".pdf"),width = 4,height = 4)
      ggsave(p,filename = paste0(save_plot,"mac_",use_genes[i],"_",chip[j],".png"),width = 4,height = 4)
    }, error = function(e) {

      print(paste0(use_genes[i],"_",chip[j]))
    })
    
    print(result)
  }
  print(paste0("Done","_",use_genes[i]))
}

