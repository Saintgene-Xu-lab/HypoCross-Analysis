if(T){
  library(miloR)
  library(Seurat)
  library(SeuratWrappers)
  library(SingleCellExperiment)
  library(scater)
  library(dplyr)
  library(patchwork)
  library(qs)
}

setwd("")
meta <- readRDS("meta_20230904_v4.1.rds")
hy_neuron <- readRDS("")
hy_neuron
use_cell <- qread("")
DefaultAssay(hy_neuron) <- "RNA3"
meta <- meta[colnames(hy_neuron),]
table(colnames(hy_neuron) == rownames(meta))
hy_neuron$celltype <- hy_neuron$anno6
hy_neuron <- subset(hy_neuron,cells = use_cell)
hy_neuron

hy_neuron[["RNA"]] <- JoinLayers(hy_neuron[["RNA"]])
hy_neuron[["RNA3"]] <- as(object = hy_neuron[["RNA"]], Class = "Assay")

DefaultAssay(hy_neuron) <- "RNA3"
hy_neuron$celltype <- hy_neuron$anno6
set.seed(12580)
DefaultAssay(hy_neuron) <- "RNA3"

############# mouse vs monkey
hy_neuron$specie2 <- as.character(hy_neuron$specie)
hy_neuron@meta.data[hy_neuron$specie %in% c("macaque","marmoset"), ]$specie2 <- "monkey"
qsave(hy_neuron,file = "./milor_seu.qs")
rm(list=ls())

#### ------------------------------------------                 ---------------------------------------------------------
#### species-----------------------------------
hy_neuron$specie2 <- as.character(hy_neuron$specie2)
hy_neuron$cell_spec <- paste0(hy_neuron$specie2,"_",hy_neuron$celltype)
Idents(hy_neuron) <- hy_neuron$cell_spec


#### 
sce <- SingleCellExperiment(hy_neuron@assays$RNA3@counts)
sce@colData <- DataFrame(hy_neuron@meta.data)
reducedDims(sce) <- list(PCA=hy_neuron@reductions$integrated.rpca@cell.embeddings, UMAP=hy_neuron@reductions$umap@cell.embeddings)
reducedDimNames(sce)

#### milor caculate ---------------------------------------
### constract objucts 
sce_milo <- miloR::Milo(sce)
sce_milo  <- buildGraph(sce_milo , k = 30, d = 30,reduced.dim = "PCA")
### Defining representative neighbourhoods on the KNN graph
sce_milo <- makeNhoods(sce_milo, prop = 0.1, k = 30, d=30, refined = TRUE, reduced_dims = "PCA")
# it’s best to have a distribution peaking between 50 and 100.Otherwise you might consider rerunning makeNhoods increasing k and/or prop
plotNhoodSizeHist(sce_milo)

### Counting cells in neighbourhoods  !! "cell_spec" : 
sce_milo  <- countCells(sce_milo, meta.data = as.data.frame(colData(sce_milo)), sample="cell_spec")

### Differential abundance testing ！！ 
neu_design <- data.frame(colData(sce_milo))[,c("cell_spec", "specie2")]
neu_design <- distinct(neu_design)
rownames(neu_design) <- neu_design$cell_spec
### Milo uses an adaptation of the Spatial FDR correction introduced by cydar, which accounts for the overlap between neighbourhoods. 
sce_milo <- calcNhoodDistance(sce_milo, d=30)

da_results <- testNhoods(sce_milo, design = ~specie2, design.df = neu_design,norm.method = "logMS")
da_results %>%
  arrange(- SpatialFDR) %>%
  head() 
#da_results <- annotateNhoods(neu_milo, da_results, coldata_col = "cell_spec") 
da_results <- annotateNhoods(sce_milo, da_results, coldata_col = "celltype")
head(da_results)

sce_milo <- buildNhoodGraph(sce_milo)

plotUMAP(sce_milo)
plotNhoodGraphDA(sce_milo, da_results, alpha=0.05) +
  plot_layout(guides="collect") # 4 5
plotNhoodGraphDA(sce_milo, da_results, alpha=0.05) +
  plot_layout(guides="collect")+
  scale_fill_gradientn(colours = c("#F44336","#EF5350","white","#4FC3F7","#03A9F4"))

saveRDS(sce_milo,file = paste0(path,"neu_milo_mou_mankey.rds"))
saveRDS(da_results,file = paste0(path,"da_results_mou_mankey.rds"))
qsave(sce_milo,file = paste0(path,"neu_milo_mou_mankey.qs"))


p <- plotNhoodGraphDA(sce_milo, da_results, alpha=0.05) +
  plot_layout(guides="collect") # 4 5
ggsave(p,filename = paste0(path,"mou_mankey.png"),width = 10,height = 8)
ggsave(p,filename = paste0(path,"mou_mankey.pdf"),width = 10,height = 8)

p_data <- p$data
qsave(p_data,file = paste0(path,"p.qs"))
p2 <- plotDAbeeswarm(da_results, group.by = "celltype")
p_data <- p2$data
qsave(p_data,file = paste0(path,"p2.qs"))

rm(list=ls())

