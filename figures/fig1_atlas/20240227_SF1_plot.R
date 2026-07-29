if(T){
  library(Seurat)
  library(qs)
  library(ggplot2)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(reshape2)
  library(RColorBrewer)
}


path <- ""

meta <- qread("")


color <- c("#13a983","#44c1f0","#3f60aa","#e8490f","#f18800","#ffa631","#9ec417","#f47983","#fcbad3","#9b5de5","#cca4e3","#7bcfa6")
names(color) <- c("GABA","GLU","CHAT","AS","OPC","OD","TC","EpC","MG","END","VLMC","PTC")
meta$major_class <- factor(meta$major_class,levels = c("GABA","GLU","CHAT","AS","OPC","OD","TC","EpC","MG","END","VLMC","PTC"))

mou_meta <- meta[meta$specie %in% "mouse",]
mac_meta <- meta[meta$specie %in% "macaque",]
mar_meta <- meta[meta$specie %in% "marmoset",]

color3  <- c("#e20612","#F57F17","#00b0eb")
names(color3) <- c("macaque","marmoset","mouse")
### umi 
plot_violin <- function(data,color){
  p <- data %>% ggplot(aes(x = major_class,y = log10(nCount_RNA),fill = major_class,color =major_class))+
    geom_violin()+
    scale_color_manual(values = color)+
    scale_fill_manual(values = color)+
    theme_bw()+
    theme(legend.position = "none")
  return(p)
}
umi_mac <- plot_violin(data = mac_meta,color = color)
umi_mar <- plot_violin(data = mar_meta,color = color)
umi_mou <- plot_violin(data = mou_meta,color = color)
umi_mac/umi_mar/umi_mou

umi_spe <- meta %>% ggplot(aes(x = specie,y = log10(nCount_RNA),fill = specie,color = specie))+
  geom_violin()+
  scale_color_manual(values = color3)+
  scale_fill_manual(values = color3)+
  theme_bw()+
  theme(legend.position = "none")
umi_spe

### gene 
plot_violin <- function(data,color){
  p <- data %>% ggplot(aes(x = major_class,y = log10(nFeature_RNA),fill = major_class,color =major_class))+
    geom_violin()+
    scale_color_manual(values = color)+
    scale_fill_manual(values = color)+
    theme_bw()+
    theme(legend.position = "none")
  return(p)
}
gene_mac <- plot_violin(data = mac_meta,color = color)
gene_mar <- plot_violin(data = mar_meta,color = color)
gene_mou <- plot_violin(data = mou_meta,color = color)
gene_mac/gene_mar/gene_mou

gene_spe <- meta %>% ggplot(aes(x = specie,y = log10(nFeature_RNA),fill = specie,color = specie))+
  geom_violin()+
  scale_color_manual(values = color3)+
  scale_fill_manual(values = color3)+
  theme_bw()+
  theme(legend.position = "none")
gene_spe

### 
plot_violin <- function(data,color){
  p <- data %>% ggplot(aes(x = major_class,y = mt,fill = major_class,color =major_class))+
    geom_violin()+
    scale_color_manual(values = color)+
    scale_fill_manual(values = color)+
    theme_bw()
  return(p)
}
mt_mac <- plot_violin(data = mac_meta,color = color)
mt_mar <- plot_violin(data = mar_meta,color = color)
mt_mou <- plot_violin(data = mou_meta,color = color)
mt_mac/mt_mar/mt_mou

p <- umi_spe+umi_mac+umi_mar+umi_mou+
  plot_layout(widths = c(0.6,4,4,4))
p
ggsave(p,filename = paste0(path,"umi2.pdf"),width = 20,height = 2.5)
ggsave(p,filename = paste0(path,"umi.pdf"),width = 20,height = 3)

p1 <- gene_spe+gene_mac+gene_mar+gene_mou+
  plot_layout(widths = c(0.6,4,4,4))
p1
ggsave(p1,filename = paste0(path,"gene2.pdf"),width = 20,height = 2.5)



### 
meta$major_class <- factor(meta$major_class,levels = c("GABA","GLU","CHAT","AS","OPC","OD","TC","EpC","MG","END","VLMC","PTC"))
meta2 <- meta[meta$sample_area %in% "HY",]
set.seed(12580)
mac <- sample(rownames(meta2[meta2$specie %in% "macaque",]),size = 80000,replace = F)
mar <- sample(rownames(meta2[meta2$specie %in% "marmoset",]),size = 80000,replace = F)
mou <- sample(rownames(meta2[meta2$specie %in% "mouse",]),size = 80000,replace = F)
meta3 <- meta2[c(mac,mar,mou),]

Cellratio <- prop.table(table(meta3$major_class, meta3$specie), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
p2 <- ggplot(data = Cellratio, aes(x =Var1, y = Freq, fill =  Var2)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=color3) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
p2

p1 <- ggplot(data = Cellratio, aes(x =Var2, y = Freq, fill =  Var1)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
p1

p3 <- ggplot(data = Cellratio, aes(x =Var2, y = Freq, fill =  Var1)) +
  geom_bar(stat = "identity", width=0.8,position=position_dodge(0.75))+
  scale_fill_manual(values=color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
p3

p4 <- ggplot(data = Cellratio, aes(x =Var1, y = Freq, fill =  Var2)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=color3) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
p2


### label transfer 
mac_to_mou <- qread("")
mar_to_mou <- qread("")

mac_to_mou$major_class <- paste0("mou_",mac_to_mou$major_class)
mac_to_mou$major_class <- factor(mac_to_mou$major_class,levels = paste0("mou_",c("GABA","GLU","CHAT","AS","OPC","OD","TC","EpC","MG","END","VLMC","PTC")))
mac_to_mou$predicted.id <- factor(mac_to_mou$predicted.id,levels = c("GABA","GLU","CHAT","AS","OPC","OD","TC","EpC","MG","END","VLMC","PTC"))
Cellratio <- prop.table(x = table(mac_to_mou$major_class,mac_to_mou$predicted.id),margin = 2)
col <- colorRampPalette(c("#FAFAFA","#EEEEEE","#757575","#212121"))(100)
pheatmap::pheatmap(mat = Cellratio,
                   cluster_rows = F,
                   cluster_cols = F,
                   border_color = "white",color = col ) # 4 5


mar_to_mou$major_class <- paste0("mou_",mar_to_mou$major_class)
mar_to_mou$major_class <- factor(mar_to_mou$major_class,levels = paste0("mou_",c("GABA","GLU","CHAT","AS","OPC","OD","TC","EpC","MG","END","VLMC","PTC")))
mar_to_mou$predicted.id <- factor(mar_to_mou$predicted.id,levels = c("GABA","GLU","CHAT","AS","OPC","OD","TC","EpC","MG","END","VLMC","PTC"))
Cellratio <- prop.table(x = table(mar_to_mou$major_class,mar_to_mou$predicted.id),margin = 2)
col <- colorRampPalette(c("#FAFAFA","#EEEEEE","#757575","#212121"))(100)
pheatmap::pheatmap(mat = Cellratio,
                   cluster_rows = F,
                   cluster_cols = F,
                   border_color = "white",color = col ) # 4 5



