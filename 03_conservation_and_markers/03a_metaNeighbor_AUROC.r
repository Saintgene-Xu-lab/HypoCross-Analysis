
##################### F1 
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
  library(COSG)
  library(tidyverse)
  library(RColorBrewer)
  library(COSG)
}
all_cell <- qread("all_cell.qs")

#### pse bulk
all_cell$cell_spe <- paste0(all_cell$anno2,"_",all_cell$specie)
Idents(all_cell) <- all_cell$cell_spe
pse_bulk <- AverageExpression(object = all_cell)
pse_bulk2 <- data.frame(pse_bulk$integrated)
meta <- data.frame(cell = unlist(lapply(str_split(colnames(pse_bulk2),"_"), function(x){x <- x[1]})),
                   specie =unlist(lapply(str_split(colnames(pse_bulk2),"_"), function(x){x <- x[2]})) )
rownames(meta) <- colnames(pse_bulk2)


############## auroc start ------------------------------------------------------------------------------------
##########  -------------------------------------------

cells_spe <- data.frame(table(all_cell$cell_spe))
Idents(all_cell) <- all_cell$cell_spe
get_cell_id <- function(data,cells,less_n){
  cell_sub <- as.character(cells[cells$Freq >= less_n,]$Var1)
  cell_rem <- setdiff(levels(cells$Var1),cell_sub)
  
  cell_sub_id <- list()
  for (i in 1:length(cell_sub)) {
    set.seed(12580)
    cell_sub_id[[cell_sub[i]]] <- sample(x = WhichCells(object = data,idents = cell_sub[i]),size = less_n,replace = F)
  }
  # cell_rem_id <- list()
  # for (i in 1:length(cell_rem)) {
  #   cell_rem_id[[cell_rem[i]]] <- WhichCells(object = data,idents = cell_rem[i])
  # }
  cells_id <- c(as.character(unlist(cell_sub_id))) # ,as.character(unlist(cell_rem_id))
  return(cells_id)
}
cell_id <- get_cell_id(data = all_cell,cells = cells_spe,less_n = 50)

seu <- subset(all_cell,cells = cell_id)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu,nfeatures = 1000)
hvg2 <- VariableFeatures(seu)
seu <- as.SingleCellExperiment(seu)

Idents(all_cell) <- all_cell$anno2
cosg_marker2 <- cosg(
  all_cell,
  groups='all',
  assay='SCT',
  slot='data',
  mu=1,
  n_genes_user=100)
Idents(all_cell) <- all_cell$cell_spe
pse_bulk  <- AverageExpression(object = all_cell)
pse_bulk2 <- data.frame(pse_bulk$RNA)
top_gene <- unique(as.character(as.matrix(cosg_marker2$names[1:5,])))
# all_cell_merge <- FindVariableFeatures(all_cell_merge,nfeatures = 1000)
# hvgs <- VariableFeatures(all_cell_merge)

celltype_NV = MetaNeighborUS(var_genes=top_gene, # hvg2[1:200]
                             dat=seu,
                             study_id=seu$specie,
                             cell_type=seu$anno2)
pheatmap::pheatmap(celltype_NV)

cell_anno <- c("AS_macaque","AS_marmoset","AS_mouse","CHAT_macaque","CHAT_marmoset","CHAT_mouse","END_macaque","END_marmoset","END_mouse",
               "EpC_macaque","EpC_marmoset","EpC_mouse","GABA_macaque","GABA_marmoset","GABA_mouse","GLU_macaque","GLU_marmoset","GLU_mouse","MG_macaque",
               "MG_marmoset","MG_mouse","OD_macaque","OD_marmoset","OD_mouse","OPC_macaque","OPC_marmoset","OPC_mouse","PTC_macaque","PTC_marmoset",
               "PTC_mouse","TC_macaque","TC_marmoset","TC_mouse","VLMC_macaque","VLMC_marmoset","VLMC_mouse")


colnames(celltype_NV) <- cell_anno
rownames(celltype_NV) <- cell_anno

pheatmap::pheatmap(celltype_NV[cell_order,cell_order])
pheatmap::pheatmap(celltype_NV[cell_order,cell_order],color = colorRampPalette(brewer.pal(9,'RdPu')[c(1:7)])(100))


group <- data.frame(cell = unlist(lapply(str_split(colnames(pse_bulk2),"_"), function(x){x <- x[1]})),
                    specie = unlist(lapply(str_split(colnames(pse_bulk2),"_"), function(x){x <- x[2]})))
rownames(group) <- colnames(pse_bulk2)
celltype_NV <- celltype_NV[cell_order,cell_order]

group$cell <- factor(group$cell,levels = c("GABA","GLU","CHAT","AS","OPC","OD","TC","EpC","MG","END","VLMC","PTC"))
group$specie <- factor(group$specie,levels = c("macaque", "marmoset", "mouse"))
color <- c("#13a983","#44c1f0","#3f60aa","#e8490f","#f18800","#ffa631","#9ec417","#f47983","#fcbad3","#9b5de5","#cca4e3","#7bcfa6")
names(color) <- c("GABA","GLU","CHAT","AS","OPC","OD","TC","EpC","MG","END","VLMC","PTC")
color3  <- c("#e20612","#F57F17","#00b0eb")
names(color3) <-  c("macaque", "marmoset", "mouse")
p_auroc <- ComplexHeatmap::pheatmap(mat = celltype_NV,border_color = "white",annotation_col = as.matrix(group),
                                    annotation_row = as.matrix(group),
                                    annotation_colors = list(specie = color3,cell = color),
                               color = colorRampPalette(c("#FAFAFA",brewer.pal(9,'GnBu')[c(1:5,8)]))(100),
                               cluster_rows = F,cluster_cols = F)
p_auroc
dot_group <- group 
dot_group$name <- rownames(dot_group)

melt(p_auroc@matrix) %>% ggplot(aes(x = Var1,y = Var2,fill = value))+
  geom_tile(color = "white",size = 0.5)+
  scale_size(range = c(0.1,5),limits = c(0.1,0.5))+
  scale_fill_gradientn(colours = colorRampPalette(c("#FAFAFA",brewer.pal(9,'GnBu')[c(1:5,8)]))(100))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        axis.text.y = element_blank()) + labs(y = NULL)-> p_auroc2
p_auroc2

############## auroc end ------------------------------------------------------------------------------------


###################### ------------------------------------------------ cor  start --------------------------------------------------------------------
Idents(all_cell) <- all_cell$anno2
cosg_marker2 <- cosg(
  all_cell,
  groups='all',
  assay='SCT',
  slot='data',
  mu=1,
  n_genes_user=100)
Idents(all_cell) <- all_cell$cell_spe
pse_bulk  <- AverageExpression(object = all_cell)
pse_bulk2 <- data.frame(pse_bulk$RNA)
top_gene <- unique(as.character(as.matrix(cosg_marker2$names[1:30,])))

pse_bulk2 <- na.omit(pse_bulk2[top_gene,])
#pse_bulk2 <- pse_bulk2[hvgs,]
seu_ave2_cor <- cor(pse_bulk2,method = "pearson" )

pheatmap::pheatmap(seu_ave2_cor[cell_order,cell_order],border_color = "white",color =  colorRampPalette(brewer.pal(9,'GnBu')[1:8])(100)) # YlGnBu
pheatmap::pheatmap(seu_ave2_cor[cell_order,cell_order],
                   border_color = "white",cluster_rows = F,cluster_cols = F,
                   color =  colorRampPalette(brewer.pal(9,'GnBu')[1:8])(100)) # YlGnBu


group <- data.frame(cell = unlist(lapply(str_split(colnames(pse_bulk2),"_"), function(x){x <- x[1]})),
                    specie = unlist(lapply(str_split(colnames(pse_bulk2),"_"), function(x){x <- x[2]})))
rownames(group) <- colnames(pse_bulk2)

group$cell <- factor(group$cell,levels = c("GLU","GABA","CHAT","PTC","OPC","OD","TC","AS","VLMC","MG","END","EpC"))
group$specie <- factor(group$specie,levels = c("macaque", "marmoset", "mouse"))
color <- c("#13a983","#44c1f0","#3f60aa","#e8490f","#f18800","#ffa631","#9ec417","#f47983","#fcbad3","#9b5de5","#cca4e3","#7bcfa6")
names(color) <- c("GABA","GLU","CHAT","AS","OPC","OD","TC","EpC","MG","END","VLMC","PTC")
color3  <- c("#e20612","#F57F17","#00b0eb")
names(color3) <-  c("macaque", "marmoset", "mouse")

library(aplot)#拼图

p2 <- ComplexHeatmap::pheatmap(seu_ave2_cor,annotation_colors = list(cell = color,specie = color3),
                              border_color = "white",annotation_row= group,annotation_col = group,
                              color =  colorRampPalette(brewer.pal(9,'GnBu')[1:8])(100))

dot_group <- group 
dot_group$name <- rownames(dot_group)

dat <- melt(p2@matrix)
melt(p2@matrix) %>% ggplot(aes(x = Var1,y = Var2,color = value))+
  geom_point(aes(size = value))+
  scale_size(range = c(1,5))+
  # scale_color_gradientn(colours = colorRampPalette(brewer.pal(9,'GnBu')[c(1:8)])(100))+ # 
  scale_color_gradientn(colours = colorRampPalette(brewer.pal(9,'RdPu')[c(1:7)])(100))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        axis.text.y = element_blank()) + labs(y = NULL)-> p3
p3
library(patchwork)
p3 %>%
  insert_left(p_left,width=.3) %>%
  insert_top(p_up,height = .1)

anno_spe <- dot_group %>% as.data.frame() %>% mutate(p="")%>%
  ggplot(aes(p,name,fill=specie))+
  geom_tile() + 
  scale_y_discrete(position="right") +
  scale_fill_manual(values = color3)+
  theme_minimal()+xlab(NULL) + ylab(NULL)+
  labs(fill = "Type_B")
anno_spe
anno_cell <- dot_group %>% as.data.frame() %>% mutate(p="")%>%
  ggplot(aes(p,name,fill=cell))+
  geom_tile() + 
  scale_y_discrete(position="right") +
  scale_fill_manual(values = color)+
  theme_minimal()+xlab(NULL) + ylab(NULL)+
  labs(fill = "Type_B")
anno_cell

p3 %>%
  insert_left(p_left,width=.3) %>%
  insert_top(p_up,height = .1) -> p3_plot # 8 10

p3 %>%
  insert_left(p_left,width=.3) %>%
  insert_top(p_up,height = .1) %>%
  insert_left(anno_spe,width=.02) %>% 
  insert_left(anno_cell,width=.02) -> p3_plot2
##################################################   cor   end --------------------------------------------------------------------------

################################################## umap start ---------------------------------------------------------------------------------

group
test <- CreateSeuratObject(counts = data.frame(pse_bulk$RNA[top_gene,]),meta.data = group)
test <- NormalizeData(test)
test <- FindVariableFeatures(test)
test <- ScaleData(test)
test <- RunPCA(test,npcs = 20)
test <- RunUMAP(test,dims = 1:20,min.dist = 0.5,n.neighbors = 10)# 6 15

library(Seurat)
library(ggplot2)
library(ggunchull)
library(ggrepel) # geom_label_repel()
library(tidyverse)
library(ggrastr)
library(ggplot2)

plotData <- as.data.frame(test[["umap"]]@cell.embeddings)
plotData$cluster <- test$cell
plotData$group <- test$specie
colnames(plotData)[1:2] <- c("UMAP1","UMAP2")

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)
# 
set.seed(0816)
col <- color

# 
identity_label_pos <- plotData %>%
  group_by(cluster) %>%
  summarise(
    x = median(UMAP1),
    y = median(UMAP2)
  )

# 
set.seed(111)
library(tidyquant)  
library(ggforce)
p_umap <- ggplot(plotData, aes(x = UMAP1, y = UMAP2, fill = cluster, color = cluster)) +
  
  geom_point(size = 7,aes(shape = group))+
  geom_mark_hull(aes(fill = cluster,label = cluster),
                 expand = unit(3, "mm"), radius = unit(1, "mm"),
                 concavity = 2,con.type = "elbow",con.cap = 0)+
  scale_color_manual(values = color)+
  scale_fill_manual(values = color)+
  theme_bw()
p_umap
p_umap <- ggplot(plotData, aes(x = UMAP1, y = UMAP2, fill = cluster, color = cluster)) +
  geom_point(size = 3,aes(shape = group))+
  geom_mark_hull(aes(fill = cluster,label = cluster),expand = 0.03)+ 
  scale_fill_manual(values = color)+
  scale_color_manual(values = color)+
  geom_label_repel(data = identity_label_pos,
                   aes(label = cluster, x = x, y = y),
                   fontface = "italic",
                   point.padding = unit(1, "lines"),
                   inherit.aes = F)+
  theme_bw()
p_umap # 4 5.5


library(tidyquant)
p <- ggplot(plotData, aes(x = UMAP1, y = UMAP2, fill = cluster, color = cluster)) +
  stat_unchull(alpha = 1, linewidth = 0.25, lty = 1,nbin = 100,sfac = 1.3) +
  geom_point(size = 2) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none")+
  #axis.line = element_line(arrow = arrow(type = "closed"),linewidth = 0.5),
  #axis.title = element_text(hjust = 0.05, face = "italic")) +
  #guides(color = "none", x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  scale_fill_manual(values = col) +
  scale_color_manual(values = col) +
  # 看这里——这是新增的代码，用于增加亚群标签
  geom_label_repel(data = identity_label_pos,
                   aes(label = cluster, x = x, y = y),
                   fontface = "italic",
                   point.padding = unit(0.2, "lines"),
                   inherit.aes = F)
p

saveRDS(test,file = "/umap.rds")
##############################################    umap end --------------------------------------------------------------------------
rm(list=ls())
ls()


colnames(celltype_NV) <- colnames(pse_bulk2)
rownames(celltype_NV) <- colnames(pse_bulk2)
group <- data.frame(cell = unlist(lapply(str_split(colnames(pse_bulk2),"_"), function(x){x <- x[1]})),
                    specie = unlist(lapply(str_split(colnames(pse_bulk2),"_"), function(x){x <- x[2]})))
rownames(group) <- colnames(pse_bulk2)
celltype_NV <- celltype_NV[cell_order,cell_order]

group$cell <- factor(group$cell,levels = c("GABA","GLU","CHAT","AS","OPC","OD","TC","EpC","MG","END","VLMC","PTC"))
group$specie <- factor(group$specie,levels = c("macaque", "marmoset", "mouse"))
color <- c("#13a983","#44c1f0","#3f60aa","#e8490f","#f18800","#ffa631","#9ec417","#f47983","#fcbad3","#9b5de5","#cca4e3","#7bcfa6")
names(color) <- levels(group$cell)

pheatmap::pheatmap(celltype_NV,cluster_cols = F,cluster_rows = F,color = colorRampPalette(brewer.pal(9,'RdPu')[c(1,3,9)])(20))

p2 <- ComplexHeatmap::pheatmap(mat = celltype_NV,border_color = "white",
                               annotation_legend = group,
                               cluster_rows = F,cluster_cols = F)
p2
dot_group <- group 
dot_group$name <- rownames(dot_group)

melt(p2@matrix) %>% ggplot(aes(x = Var1,y = Var2,color = value))+
  geom_point(aes(size = value))+
  scale_size(range = c(1,5))+
  # scale_color_gradientn(colours = colorRampPalette(brewer.pal(9,'GnBu')[c(1:8)])(100))+ # 
  scale_color_gradientn(colours = colorRampPalette(brewer.pal(9,'RdPu')[c(1:7)])(100))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        axis.text.y = element_blank()) + labs(y = NULL)-> p2
p2




############## ----------------------------------
Idents(all_cell) <- all_cell$cell_spe
all_cell <- BuildClusterTree(object = all_cell,dims = 1:30)
library(ggtree)
tree <- all_cell@tools$BuildClusterTree
p <- ggtree(tree)+geom_tiplab()
p

## 提取tips 的labels
cell_order <- get_taxa_name(p)

