
library(ggplot2)
library(RColorBrewer)
library(scales)
library(ggsci)
library(harmony)
library(dplyr)
library(spatstat)
library(Seurat)
cor1 <- pal_jco()(10)
cor2 <- pal_nejm()(8)
cor3 <- pal_ucscgb()(10)
cor4 <- pal_d3()(10)
cor5 <- pal_simpsons()(15)
cor6 <- pal_simpsons()(12)
cor7 <- pal_startrek()(15)

library(monocle3)
library(dplyr)
cor <- c(cor1,cor2,cor3,cor4,cor5,cor6,cor7,brewer.pal(9, "Set1"),brewer.pal(8, "Set2"),brewer.pal(12, "Paired"),brewer.pal(8, "Dark2"))
cor <- unique(cor)
oligo <- readRDS('./OD_20250305.rds')

specie_sel <- 'marmoset'
pbmc <- subset(oligo,subset=specie==specie_sel)
data <- GetAssayData(pbmc, assay ='RNA', slot = 'counts')

cell_metadata <- pbmc@meta.data # %>% select(Fetus,location,final_cluster, annotation)
gene_annotation <-data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- gene_annotation$gene_short_name


unique(pbmc@meta.data$specie)


data <- data[rownames(gene_annotation),rownames(cell_metadata)]
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds)
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="annotation",group_label_size=5) + ggtitle('cds.umap')



