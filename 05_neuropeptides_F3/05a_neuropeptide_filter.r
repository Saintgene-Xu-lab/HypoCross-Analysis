
if(T){
  library(SeuratWrappers)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ggplot2)
  library(patchwork)
  library(COSG)
  library(qs)
  options(future.globals.maxSize = 3e+09)
  options(Seurat.object.assay.version = "v5")
}


load("regulatory_database.Rdata")
neu <- qread("")
DefaultAssay(neu) <- "RNA"
neu <- NormalizeData(neu)
Idents(neu) <- neu$anno6

peti_use <- qread("use_peti.qs")
path <- ""

p <- DotPlot(neu,features = peti_use)
p_pep_stat <- p$data
p_pep_stat_use <- p$data[p$data$pct.exp > 10 & p$data$avg.exp.scaled > 0.5,]
qsave(p_pep_stat,file = paste0(path,"/p_pep_stat_all_spe.qs"))
qsave(p_pep_stat_use,file = paste0(path,"/p_pep_stat_use_all_spe.qs"))

###
specie <- names(table(neu$specie))
neu <- SplitObject(neu,split.by = "specie")
neu[["mouse"]] <- subset(neu[["mouse"]],cells = rownames(neu[["mouse"]]@meta.data[neu[["mouse"]]$sample_area %in% "HY",]))
p_pep_stat <- list()
p_pep_stat_use <- list()

for (i in 1:length(neu)) {
  p <- DotPlot(neu[[specie[i]]],features = peti_use)+theme(axis.text.x = element_text(angle = 90,hjust = 1))+coord_flip()
  ggsave(p,filename = paste0(path,"/",specie[i],"_","pepi.png"),width = 10,height = 6)
  ggsave(p,filename = paste0(path,"/",specie[i],"_","pepi.pdf"),width = 10,height = 6)
  p_pep_stat[[specie[i]]] <- p$data
  p_pep_stat[[specie[i]]]$specie <- rep(specie[i],nrow(p$data))
  p_pep_stat_use[[specie[i]]] <- p$data[p$data$pct.exp > 10 & p$data$avg.exp.scaled > 0.5,]
  p_pep_stat_use[[specie[i]]]$specie <- rep(specie[i],nrow(p_pep_stat_use[[specie[i]]]))
  
}
qsave(p_pep_stat,file = paste0(path,"/p_pep_stat.qs"))
qsave(p_pep_stat_use,file = paste0(path,"/p_pep_stat_use.qs"))

#### get gaba glu cell id
data_glu <- list()
data_gaba <- list()
data_chat <- list()
data_th <- list()
data_Slc18a2 <- list()
data_Slc6a3 <- list()
data_Gad1 <- list()
data_Gad2 <- list()

for (i in 1:length(specie)){
  data_glu[[specie[i]]] <- WhichCells(object = neu[[specie[i]]],expression = Slc17a6 > 0.5)
  data_gaba[[specie[i]]] <- WhichCells(object = neu[[specie[i]]],expression = Slc32a1 > 0.5)
  data_chat[[specie[i]]] <- WhichCells(object = neu[[specie[i]]],expression = Chat > 0.5)
  data_th[[specie[i]]] <- WhichCells(object = neu[[specie[i]]],expression = Th > 0.5)
  data_Slc18a2[[specie[i]]] <- WhichCells(object = neu[[specie[i]]],expression = Th > 0.5)
  data_Slc6a3[[specie[i]]] <- WhichCells(object = neu[[specie[i]]],expression = Th > 0.5)
  data_Gad1[[specie[i]]] <- WhichCells(object = neu[[specie[i]]],expression = Th > 0.5)
  data_Gad2[[specie[i]]] <- WhichCells(object = neu[[specie[i]]],expression = Th > 0.5)
}
qsave(data_glu,file = paste0(path,"/data_glu.qs"))
qsave(data_gaba,file = paste0(path,"/data_gaba.qs"))
qsave(data_chat,file = paste0(path,"/data_chat.qs"))
qsave(data_th,file = paste0(path,"/data_th.qs"))
qsave(data_Slc18a2,file = paste0(path,"/data_Slc18a2.qs"))
qsave(data_Slc6a3,file = paste0(path,"/data_Slc6a3.qs"))
qsave(data_Gad1,file = paste0(path,"/data_Gad1.qs"))
qsave(data_Gad2,file = paste0(path,"/data_Gad2.qs"))

rm(list=ls())

