

specie_sel <- 'macaque'
library(spatstat)
library(ggplot2)
library(dplyr)
library(Seurat)
library(maSigPro)

genes <- read.csv(paste(specie_sel,"Trajectory_genes.sc.csv",sep=''),header=T)

oligo <- readRDS('./OD_OPC_20250609.rds')
pbmc <- subset(oligo,subset=specie==specie_sel)
cds <- readRDS(paste( './monocle3_out/',specie_sel,'_monocle3/', specie_sel,'_monocle3.cds.rds',sep=''))


pseudotime <- pseudotime(cds,reduction_method = 'UMAP')
pseudotime <-pseudotime[rownames(pbmc@meta.data)]
pbmc$pseudotime <- pseudotime
p8 = FeaturePlot(pbmc, reduction ="umap", features = "pseudotime")
ggsave(paste(specie_sel,"Pseudotime_Seurat.sc.pdf",sep='_'),plot = p8, width = 8, height = 6)
saveRDS(pbmc, file =paste(specie_sel,"sco_pseudotime.sc.rds",sep='_'))
#######################################################################################################################################

pbmc <- readRDS(paste(specie_sel,"sco_pseudotime.sc.rds",sep='_'))
seurat.child <- pbmc
Idents(seurat.child) <- 'annotation'
markers <- FindAllMarkers(object = seurat.child, slot = 'data',only.pos = TRUE, logfc.threshold = 0.5)
markers <- markers[markers$p_val_adj<0.05,]
table(markers$cluster)
saveRDS(markers,file=paste(specie_sel,"OD_markers.sc.rds",sep='_') )
mac_markers <- readRDS('./macaque_monocle3/macaque_OD_markers.sc.rds')

#####___________________________________________________________

specie_sel <- 'marmoset'

pbmc <- readRDS(paste(specie_sel,"sco_pseudotime.sc.rds",sep='_'))
seurat.child <- pbmc
Idents(seurat.child) <- 'annotation'
markers <- FindAllMarkers(object = seurat.child, slot = 'data',only.pos = TRUE, logfc.threshold = 0.5)
markers <- markers[markers$p_val_adj<0.05,]
table(markers$cluster)
saveRDS(markers,file=paste(specie_sel,"OD_markers.sc.rds",sep='_') )
mar_markers <- readRDS('./marmoset_OD_markers.sc.rds')

#####___________________________________________________________

specie_sel <- 'mouse'

pbmc <- readRDS(paste(specie_sel,"sco_pseudotime.sc.rds",sep='_'))
seurat.child <- pbmc
Idents(seurat.child) <- 'annotation'
markers <- FindAllMarkers(object = seurat.child, slot = 'data',only.pos = TRUE, logfc.threshold = 0.5)
markers <- markers[markers$p_val_adj<0.05,]
table(markers$cluster)
saveRDS(markers,file=paste(specie_sel,"OD_markers.sc.rds",sep='_') )

mou_markers <- markers
#####___________________________________________________________________________________________________________________________#####


output_mou <- readRDS('./masigpro/res.cluster_mouse.rds')
output_mou <- as.data.frame(output_mou$cut)
output_mou$gene <- row.names(output_mou)
names(output_mou) <- c('cluster','gene')

output_mou_filt <- output_mou[output_mou$gene %in% mou_markers$gene,]


output_mac <- readRDS('./masigpro/res.cluster_macaque.rds')
output_mac <- as.data.frame(output_mac$cut)
output_mac$gene <- row.names(output_mac)
names(output_mac) <- c('cluster','gene')
output_mac_filt <- output_mac[output_mac$gene %in% mac_markers$gene,]

output_mar <- readRDS('./masigpro/res.cluster_marmoset.rds')
output_mar <- as.data.frame(output_mar$cut)
output_mar$gene <- row.names(output_mar)
names(output_mar) <- c('cluster','gene')
output_mar_filt <- output_mar[output_mar$gene %in% mar_markers$gene,]



# ________________________



output_mou_dn <- output_mou[output_mou$cluster %in% c('1','2'),]
output_mou_up <- output_mou[output_mou$cluster %in% c('3','4','5'),]


output_mac_dn <- output_mac[output_mac$cluster %in% c('1','2'),]
output_mac_up <- output_mac[output_mac$cluster %in% c('3','4','5'),]


output_mar_dn <- output_mar[output_mar$cluster %in% c('1','2','3','5'),]
output_mar_up <- output_mar[output_mar$cluster %in% c('4'),]

inter_dn <- intersect(intersect(output_mac_dn$gene,output_mar_dn$gene),output_mou_dn$gene)

inter_up <- intersect(intersect(output_mac_up$gene,output_mar_up$gene),output_mou_up$gene)


inter_dn_only_NHP <- setdiff(intersect(output_mac_dn$gene,output_mar_dn$gene),output_mou_dn$gene)
inter_up_only_NHP <- setdiff(intersect(output_mac_up$gene,output_mar_up$gene),output_mou_up$gene)


inter_dn_only_mouse <- setdiff(output_mou_dn$gene, intersect(output_mac_dn$gene,output_mar_dn$gene))
inter_up_only_mouse <- setdiff(output_mou_up$gene, intersect(output_mac_up$gene,output_mar_up$gene))




# ________________________


dput(setdiff(inter_dn_only_NHP ,inter_dn_only_NHP[inter_dn_only_NHP %in% mou_markers[mou_markers$cluster %in% c('COPC','OPC','DOPC'),]$gene])) # NHP dn 

dput(setdiff(inter_up_only_NHP ,inter_up_only_NHP[inter_up_only_NHP %in% mou_markers[mou_markers$cluster %in% c('OD3','OD2','OD1'),]$gene])) # NHP up

# ________________________



library(monocle3)
library(dplyr)
library(ggplot2)


nhp_dn <- c("Acvr2a", "Evi5", "Ttc28", "Gphn", "Stim1", "Rbfox2", "Nr3c2",
"Stxbp5l", "Sgip1", "Ncoa1", "Pcnx", "Ppp3ca", "Sptan1", "Slc8a1",
"Cacna2d1", "Fnbp4", "Erc1", "Arhgef12", "Fign", "Agap1", "Dgki",
"Ppp1r12a", "Hivep1", "Setbp1", "Sh3d19", "Plekha2")

nhp_up <- c("Pak2", "Cep170", "Ankrd12", "Serpini1", "Mtmr10", "Itprid2",
"Synm", "Gm12248", "Arhgef37", "Cep85")


#####################################  mouse up  ########


library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

pbmc$pseudotime_bin <- cut(
  pbmc$pseudotime,
  breaks = 100,
  labels = paste0("Bin_", 1:100)
)


table(pbmc$pseudotime_bin)

gene_set <- nhp_up

bin_avg_exp <- AverageExpression(
  pbmc,
  features = gene_set,
  group.by = "pseudotime_bin",
  assays = "RNA" 
)$RNA


bin_avg_exp_scaled <- t(scale(t(bin_avg_exp)))
pdf('../mouse_bin_avg_exp_scaled_1_up.pdf',h=2,w=7)
Heatmap(
  bin_avg_exp_scaled,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("steelblue", "white", "firebrick")),
  cluster_rows = FALSE,
  cluster_columns = FALSE, 
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_title = "Pseudotime Bins",
  row_title = "Gene Set",
  show_row_names = TRUE,
  show_column_names = FALSE,
  heatmap_legend_param = list(
    title_position = "leftcenter-rot"
  )
)
dev.off()
###########  mouse dn  ########

gene_set <- nhp_dn

bin_avg_exp <- AverageExpression(
  pbmc,
  features = gene_set,
  group.by = "pseudotime_bin",
  assays = "RNA" 
)$RNA


bin_avg_exp_scaled <- t(scale(t(bin_avg_exp)))
pdf('../mouse_bin_avg_exp_scaled_1_dn.pdf',h=2,w=7)
Heatmap(
  bin_avg_exp_scaled,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("steelblue", "white", "firebrick")),
  cluster_rows = FALSE,
  cluster_columns = FALSE, 
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_title = "Pseudotime Bins",
  row_title = "Gene Set",
  show_row_names = TRUE,
  show_column_names = FALSE,
  heatmap_legend_param = list(
    title_position = "leftcenter-rot"
  )
)
dev.off()




#####################################  macaque up  ########


specie_sel <- 'macaque'

pbmc <- readRDS(paste(specie_sel,"sco_pseudotime.sc.rds",sep='_'))
pbmc$pseudotime_bin <- cut(
  pbmc$pseudotime,
  breaks = 100,
  labels = paste0("Bin_", 1:100)
)

table(pbmc$pseudotime_bin)


gene_set <- nhp_up

bin_avg_exp <- AverageExpression(
  pbmc,
  features = gene_set,
  group.by = "pseudotime_bin",
  assays = "RNA" 
)$RNA


bin_avg_exp_scaled <- t(scale(t(bin_avg_exp)))
pdf('../macaque_bin_avg_exp_scaled_1_up.pdf',h=2,w=7)
Heatmap(
  bin_avg_exp_scaled,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("steelblue", "white", "firebrick")),
  cluster_rows = FALSE,
  cluster_columns = FALSE, 
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_title = "Pseudotime Bins",
  row_title = "Gene Set",
  show_row_names = TRUE,
  show_column_names = FALSE,
  heatmap_legend_param = list(
    title_position = "leftcenter-rot"
  )
)
dev.off()
#####################  macaque dn  ########


gene_set <- nhp_dn

bin_avg_exp <- AverageExpression(
  pbmc,
  features = gene_set,
  group.by = "pseudotime_bin",
  assays = "RNA"
)$RNA

bin_avg_exp_scaled <- t(scale(t(bin_avg_exp)))
pdf('../macaque_bin_avg_exp_scaled_1_dn.pdf',h=2,w=7)
Heatmap(
  bin_avg_exp_scaled,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("steelblue", "white", "firebrick")),
  cluster_rows = FALSE,
  cluster_columns = FALSE, 
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_title = "Pseudotime Bins",
  row_title = "Gene Set",
  show_row_names = TRUE,
  show_column_names = FALSE,
  heatmap_legend_param = list(
    title_position = "leftcenter-rot"
  )
)
dev.off()



#####################################  marmoset up  ########


specie_sel <- 'marmoset'

# genes <- read.csv(paste(specie_sel,"Trajectory_genes.sc.csv",sep=''),header=T)

pbmc <- readRDS(paste(specie_sel,"sco_pseudotime.sc.rds",sep='_'))
pbmc$pseudotime_bin <- cut(
  pbmc$pseudotime,
  breaks = 100,
  labels = paste0("Bin_", 1:100)
)

table(pbmc$pseudotime_bin)


gene_set <- nhp_up

bin_avg_exp <- AverageExpression(
  pbmc,
  features = gene_set,
  group.by = "pseudotime_bin",
  assays = "RNA" 
)$RNA


bin_avg_exp_scaled <- t(scale(t(bin_avg_exp)))
pdf('../marmoset_bin_avg_exp_scaled_1_up.pdf',h=2,w=7)
Heatmap(
  bin_avg_exp_scaled,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("steelblue", "white", "firebrick")),
  cluster_rows = FALSE,
  cluster_columns = FALSE, 
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_title = "Pseudotime Bins",
  row_title = "Gene Set",
  show_row_names = TRUE,
  show_column_names = FALSE,
  heatmap_legend_param = list(
    title_position = "leftcenter-rot"
  )
)
dev.off()
#####################  marmoset dn  ########

gene_set <- nhp_dn

bin_avg_exp <- AverageExpression(
  pbmc,
  features = gene_set,
  group.by = "pseudotime_bin",
  assays = "RNA" 
)$RNA


bin_avg_exp_scaled <- t(scale(t(bin_avg_exp)))
pdf('../marmoset_bin_avg_exp_scaled_1_dn.pdf',h=2,w=7)
Heatmap(
  bin_avg_exp_scaled,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("steelblue", "white", "firebrick")),
  cluster_rows = FALSE,
  cluster_columns = FALSE, 
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_title = "Pseudotime Bins",
  row_title = "Gene Set",
  show_row_names = TRUE,
  show_column_names = FALSE,
  heatmap_legend_param = list(
    title_position = "leftcenter-rot"
  )
)
dev.off()

