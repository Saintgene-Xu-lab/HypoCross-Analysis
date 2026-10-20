if(T){
  library(Seurat)
  library(ggplot2)
  library(reshape2)
  library(qs)
  library(reshape2)
  library(clusterProfiler)
  library(DOSE)
  library(enrichplot)
  library(ggraph)
  library(cowplot)
}
#### color code
if(T){
  color3  <- c("#e20612","#F57F17","#00b0eb")
  names(color3) <- c("macaque","marmoset","mouse")
  
  colors_4 <-c("#FF9966","#00CCCC","#6666CC","#99CC33")
  names(colors_4) <- paste0("MG",1:4)
}

mg <- qread(file = "")
Idents(mg) <- mg$cell

#####
################# ---------------4 AD mg ----
######## AD MG 
if(F){
  mg_ad <- readRDS("")
  mg_meta <- readRDS("")
  human_mouse_gene <- read.csv("Human_mouse_gene.csv")
  human_mouse_gene <- human_mouse_gene[human_mouse_gene$Symbol.y %in% rownames(mg_ad),]
  human_mouse_gene <- human_mouse_gene[!duplicated(human_mouse_gene$Symbol.y),]
  human_mouse_gene <- human_mouse_gene[!duplicated(human_mouse_gene$Symbol.x),]
  rownames(human_mouse_gene) <- human_mouse_gene$Symbol.y
  mg_ad <- mg_ad[human_mouse_gene$Symbol.y,]
  table(rownames(mg_ad) == rownames(human_mouse_gene))
  rownames(mg_ad) <- human_mouse_gene$Symbol.x
  mg_ad <- CreateSeuratObject(counts = mg_ad,meta.data = mg_meta)
  
  mg_ad <- NormalizeData(mg_ad)
  mg_ad$anno <- paste0("MG",mg_ad$seurat_clusters)
  
  Idents(mg_ad) <- mg_ad$anno
  mg_ad <- subset(mg_ad,idents = paste0("MG",c(0,1:8,10:12)))
  gc()
  mg_ad <- NormalizeData(mg_ad)
  gc()
  mg_ad[["RNA"]] <- as(mg_ad[["RNA"]],Class = "Assay")
  library(COSG)
  mg_marker <- cosg(
    mg_ad,
    groups='all',
    assay="RNA",
    slot="data",
    mu=1,
    n_genes_user=500)
  
  get_top_genes <- function(top,cos,data){
    marker <-  data$names
    score <- data$scores
    cells <- colnames(marker)
    marker_cell_list <- list()
    for (i in 1:length(cells)){
      if (length(which(score[,cells[i]] > cos)) >= top) {
        marker_cell_list[[cells[i]]] <- marker[,cells[i]][1:top]
      }else{
        marker_cell_list[[cells[i]]] <- marker[,cells[i]][which(score[,cells[i]] >= cos)]
      }
    }
    return(marker_cell_list )
  }
  mg_feature <- get_top_genes(top =100,cos = 0.1,data = mg_marker)

  DefaultAssay(mg_ad) <- "RNA"
  deg_ad <- FindAllMarkers(mg_ad,min.pct = 0.25,logfc.threshold = 0.25,only.pos = T)
  deg_ad2 <- deg_ad[deg_ad$avg_log2FC > 0.5& deg_ad$p_val_adj < 0.01,]
  deg_ad2$cluster <- as.character(deg_ad2$cluster)
  table(deg_ad2$cluster )
  deg_ad2 <- deg_ad2[!deg_ad2$cluster %in% "MG1",]
  DefaultAssay(mg) <- "RNA"
  
  path
  save(deg_ad2,p_scale_data,file = paste0(path,"deg_ad2.Rdata"))
}
path <- ""
load(paste0(path,"deg_ad2.Rdata"))

mg_class <- names(table(deg_ad2$cluster))
for (i in 1:length(mg_class)) {
  mg <- AddModuleScore(object = mg,features = list(intersect(deg_ad2[deg_ad2$cluster %in% mg_class[i],]$gene,rownames(mg))),name = paste0(mg_class[i]))
}
color3  <- c("#e20612","#F57F17","#00b0eb")
names(color3) <- c("macaque","marmoset","mouse")
DefaultAssay(mg) <- "RNA"
Idents(mg) <- mg$cell
p_scale <- DotPlot(mg,features = paste0(mg_class,"1"))
p_scale$data[,c("avg.exp.scaled","features.plot","id")] %>% spread(key = features.plot,value = avg.exp.scaled) -> p_scale_data
colnames(p_scale_data) <- c("Group","MG0:Homeostatic","MG10:Inflammatory III+++","MG11:Antiviral","MG12:Cycling",
                            "MG2:InflammatoryI++","MG3:Ribosome Biogenesis","MG4:Lipid Processing","MG5:Phagocytic",
                            "MG6:Stress Signature","MG7:Glycolytic","MG8:Inflammatory II+")

plot_use <- melt(p_scale_data) 
plot_use <- plot_use[!(plot_use$variable %in% c("MG12:Cycling")),]
plot_use$variable <- factor(plot_use$variable,levels = c("MG0:Homeostatic","MG3:Ribosome Biogenesis","MG4:Lipid Processing",
                                                         "MG7:Glycolytic","MG5:Phagocytic","MG11:Antiviral",
                                                         "MG8:Inflammatory II+","MG2:InflammatoryI++",
                                                         "MG10:Inflammatory III+++","MG6:Stress Signature"))


plot_use %>% ggplot(aes(x = variable, y = value,group = Group,color = Group))+
  geom_point(size = 3)+
  geom_line(linetype = "dashed",linewidth = 0.5)+
  theme_bw()+
  scale_color_manual(values = colors_4)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.3))-> p
p
path <- ""
ggsave(p,filename = paste0(path,"mg_ad_state_anno_2.pdf"),width = 6,height = 4)




