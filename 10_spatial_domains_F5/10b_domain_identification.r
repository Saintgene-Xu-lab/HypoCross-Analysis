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

options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 10000000000000000 * 1024^2)


path <- ""
rm_cell <- qread("rm_cells.qs")

data <- qread("three_specie_merge_data.qs")
data <- subset(data,cells = intersect(colnames(data),rm_cell),invert =T)
DefaultAssay(data) <- "RNA"

data <- NormalizeData(data)
data[["RNA"]] <- split(data[["RNA"]], f = data$specie)
data <- FindVariableFeatures(data, verbose = FALSE,nfeatures = 3000)
data <- SketchData(object = data, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")

#### Perform integration on the sketched bins across samples
DefaultAssay(data) <- "sketch"
data <- FindVariableFeatures(data, verbose = F,nfeatures = 3000)
data <- ScaleData(data, verbose = F)
data <- RunPCA(data, verbose = F)
#### Perform integration on the sketched bins across samples
data <- IntegrateLayers(data, method = RPCAIntegration, orig = "pca", new.reduction = "integrated.rpca",
                        dims = 1:50, k.anchor = 50,
                        verbose = F)
# cluster the integrated data
data <- FindNeighbors(data, reduction = "integrated.rpca", dims = 1:50)
data <- FindClusters(data, resolution = 2)
data <- RunUMAP(data, reduction = "integrated.rpca", dims = 1:50, return.model = T, verbose = F,reduction.name = "umap.rpca")

########## mapping full 
data <- ProjectIntegration(object = data, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")
data <- ProjectData(object = data, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.rpca.full",
                    full.reduction = "integrated.rpca.full", dims = 1:50, refdata = list(celltype.full = "seurat_clusters"))
data <- RunUMAP(data, reduction = "integrated.rpca.full", dims = 1:50, reduction.name = "umap.full",
                reduction.key = "UMAP_full_")
## 
qsave(data,file = "int_specie_bin100.qs")

