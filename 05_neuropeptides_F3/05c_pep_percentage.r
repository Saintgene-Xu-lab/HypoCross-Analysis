if(T){
  library(SeuratWrappers)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(Seurat)
  library(SeuratData)
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
Idents(neu) <- neu$specie

specie <- names(table(neu$specie))
neu <- SplitObject(neu,split.by = "specie")
neu[["mouse"]] <- subset(neu[["mouse"]],cells = rownames(neu[["mouse"]]@meta.data[neu[["mouse"]]$sample_area %in% "HY",]))

####　
peti_use2 <- qread("")
peti_use <- unique(as.character(peti_use2$features.plot))

path <- ""
p_pep_stat <- list()
p_pep_stat_use <- list()


for (i in 1:length(neu)) {
  neu[[specie[i]]] <- NormalizeData(neu[[specie[i]]])
  p <- DotPlot(neu[[specie[i]]],features = peti_use)+theme(axis.text.x = element_text(angle = 90,hjust = 1))+coord_flip()
  p_pep_stat[[specie[i]]] <- p$data
  p_pep_stat[[specie[i]]]$specie <- rep(specie[i],nrow(p$data))
  p_pep_stat_use[[specie[i]]] <- p$data[p$data$pct.exp > 10 & p$data$avg.exp.scaled > 0.5,]
  p_pep_stat_use[[specie[i]]]$specie <- rep(specie[i],nrow(p_pep_stat_use[[specie[i]]]))
}
qsave(p_pep_stat,file = paste0(path,"p_pep_distri_in_three_specie.qs"))
qsave(p_pep_stat_use,file = paste0(path,"p_pep_distri_in_three_specie_use.qs"))
rm(lsit=ls()
   

