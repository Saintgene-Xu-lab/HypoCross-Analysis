if(T){
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(reshape2)
  library(stringr)
  library(scales)
  library(RColorBrewer)
  library(patchwork)
  library(ggtree)
  library(ggtreeExtra)
  library(COSG)
  library(qs)
  library(stringr)
  library(ggstar)
  library(tidyverse)
  library(cowplot)
  library(ape)
  library(ggtreeExtra)
  library(ggstar)
  library(ggplot2)
  library(ggtree)
  library(treeio)
  library(ggnewscale)
  library(Seurat)
  library(reshape2)
  library(RColorBrewer)
  library(scales)
  library(ggridges)
  library(igraph)
}

########  tree_data -----------------------
tree_data <- qread("tree.graph.qs")
tree_data_layout <- qread("tree.layout.qs")
a <- as_data_frame(tree_data)
colnames(a)[1:2] <- c("parent","node")
p <- ggtree(as.phylo(a,label = "to_clust"), layout='circular')+
  geom_label(aes(label = label), fill='#aa96da',label.size = 0.01)  
p
p <- ggtree(as.phylo(a,label = "to_clust"), layout='circular')+
  geom_label(aes(x=branch,label = label), fill='#aa96da',label.size = 0.01)  
p

p <- ggtree(as.phylo(a,label = "to_clust"), layout='circular')+
  geom_nodelab(aes(label = label))+
  geom_tiplab()
p 
####### meta data ----------------------------
all_cell_meta <- qread("meta.qs")
all_cell_meta <- data.frame(all_cell_meta)

color3  <- c("#e20612","#F57F17","#00b0eb")
names(color3) <- c("macaque","marmoset","mouse")
show_col(color3)

col <- c(colorRampPalette(brewer.pal(12,'Paired'))(32),colorRampPalette(brewer.pal(8,'Set1'))(30))

#meta <- all_cell@meta.data
meta <- all_cell_meta
meta2 <- meta[!(meta$sample_area %in% c("Arc","ArcC","LH","MM","POA","PVN","SCN","VMH","ZI")),]
### 
table(meta2$specie)
set.seed(12580)
mac_cell <- sample(rownames(meta2[meta2$specie %in% "macaque",]),size = 52463 ,replace = F )
set.seed(12580)
mar_cell <- sample(rownames(meta2[meta2$specie %in% "marmoset",]),size = 52463 ,replace = F )
mou_cell <- rownames(meta2[meta2$specie %in% "mouse",])
meta3 <- meta2[c(mou_cell,mac_cell,mar_cell),]

CR_species2 <- function(meta,cell,species){
  Cellratio <- prop.table(table(meta[,species],meta[,cell]), margin = 2)
  Cellratio <- as.data.frame(Cellratio)
  p <- ggplot(Cellratio) + 
    geom_bar(aes(x =Var1, y= Freq, fill = Var2),stat = "identity",width = 0.7,linewidth = 0.5,colour = 'white',position = 'fill')+ # position = 'fill'
    labs(x='Sample',y = 'Ratio')+
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90,hjust = 1))
  return(p)
} 
p2 <- CR_species2(meta = meta3,cell ="L6" ,species = "specie") # 3 12
p2
##  
cell_num <- data.frame(cell = as.numeric(table(meta3$L6)))
rownames(cell_num) <- names(table(meta3$L6))
cell_num$radio <- cell_num$cell / sum(cell_num$cell)
cell_num$cell <- rownames(cell_num)
##  
tr <- meta %>% group_by(L6) %>% summarise(umi = mean(nCount_RNA),
                                          rna = mean(nFeature_RNA))
tr <- data.frame(tr)
colnames(tr)[1] <- "cell"

###
tr <- merge(tr,cell_num,by = "cell")
tr
### marker 

all_cell <- qread("all_cell.qs")
DefaultAssay(all_cell) <- "RNA"

meta3 <- meta[intersect(colnames(all_cell),rownames(meta)),]
all_cell <- subset(all_cell,cells = rownames(meta3))
table(rownames(meta3) == colnames(all_cell))


all_cell$L1 <- meta3$L1
all_cell$L2 <- meta3$L2
all_cell$L3 <- meta3$L3
all_cell$L4 <- meta3$L4
all_cell$L5 <- meta3$L5
all_cell$L6 <- meta3$L6
all_cell$L0 <- meta3$L0
## ----------------
DefaultAssay(all_cell) <- "RNA"
all_cell <- NormalizeData(all_cell)

Idents(all_cell) <- all_cell$L6
p_dot <- DotPlot(all_cell,cluster.idents = T,
                 features = c("Pitx2","Col3a1","Cldn5","Foxj1","Pdgfra","Col23a1","Agt","Csf1r","Mog","Chat","Slc17a6","Slc32a1"))
head(p_dot$data)
marker1 <- p_dot$data[,c("avg.exp.scaled","id","features.plot")] %>%  spread(key = features.plot,value = avg.exp.scaled)
marker2 <- p_dot$data[,c("pct.exp","id","features.plot")] %>%  spread(key = features.plot,value = pct.exp)
colnames(marker2) <- paste0(colnames(marker2),"_pt")
colnames(marker2)[1] <- "cell"
colnames(marker1)[1] <- "cell"

marker <- merge(marker1,marker2,by = "cell")

tr <- merge(tr,marker,by = "cell")
## group 
if(T){
  tr$group <- "GLU"
  tr[tr$cell %in% paste0("L6_",c("TCα","TCβ")),]$group <- "TC"
  tr[tr$cell %in% paste0("L6_",c("OPC","DOPC","COPC")),]$group <- "OPC"
  tr[tr$cell %in% paste0("L6_",c("OD3","OD1","OD2")),]$group <- "OD"
  tr[tr$cell %in% paste0("L6_",c("MG3","MG1","MG2","MG4")),]$group <- "MG"
  tr[tr$cell %in% paste0("L6_AS",c(1:3)),]$group <- "AS"
  tr[tr$cell %in% c("L6_VLMC"),]$group <- "VLMC"
  tr[tr$cell %in% c("L6_PTC"),]$group <- "PTC"
  tr[tr$cell %in% c("L6_END1","L6_END2","L6_END3","L6_END4"),]$group <- "END"
  tr[tr$cell %in% c("L6_EpC1","L6_EpC2"),]$group <- "EpC"
  tr[str_detect(tr$cell,"IN"),]$group <- "GABA"
  tr[str_detect(tr$cell,"CHAT"),]$group <- "CHAT"
}

tr$group2 <- "NEU"
tr[tr$group %in% c("AS","OD","OPC","TC","MG","EpC"),]$group2 <- "Glia"
tr[tr$group %in% c("PTC","END", "VLMC"),]$group2 <- "Other"

rownames(tr) <- tr$cell

color <- c("#13a983","#44c1f0","#3f60aa","#e8490f","#f18800","#ffa631","#9ec417","#f47983","#fcbad3","#9b5de5","#cca4e3","#7bcfa6")
tr$group <- factor(tr$group,levels = c("GABA","GLU","CHAT","AS","OPC","OD","TC","EpC","MG","END","VLMC","PTC"))
names(color) <- levels(tr$group)

################ tree circle
p <- ggtree(as.phylo(a,label = "to_clust"))
new_tree <- p$data

p_tree_l6 <- ggtree(new_tree, aes(x, y),layout='circular',branch.length="none")+xlim(-10, NA)
p_tree_l6 

whitePurple <- c("#f7fcfd" ,"#e0ecf4", "#bfd3e6", "#9ebcda", "#8c96c6", "#8c6bb1", "#88419d", "#810f7c", "#4d004b")
p2_tree <- p_tree_l6 + 
  geom_fruit(
    data=tr,
    geom=geom_star,size = 2,
    mapping=aes(y=cell, fill=group,starshape=group2),#size=Length, 
    position="identity",color ="white",
    starstroke=0.5
  ) + 
  scale_fill_manual(
    values=color,
    guide="none" 
  ) + 
  scale_size_continuous(
    range=c(1, 2), # the range of size.
    guide=guide_legend(
      keywidth=0.5, 
      keyheight=0.5,
      override.aes=list(starshape=15),
      order=2
    )
  ) +
  scale_starshape_manual(
    values=c(1,15,23),
    guide=guide_legend(
      keywidth=0.5,
      keyheight=0.5,
      order=1
    )
  )
p2_tree

p3_tree <- p2_tree + 
  new_scale_fill() + 
  geom_fruit(
    data=tr,
    geom=geom_tile,
    mapping=aes(y=cell, fill=group),
    offset=0.15,   # The distance between external layers, default is 0.03 times of x range of tree.
    pwidth=0.3 # width of the external layer, default is 0.2 times of x range of tree.
  ) + 
  scale_fill_manual(
    values=color,
    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=3)
  ) 
p3_tree

p4_tree <- p3_tree + 
  new_scale_fill() +
  geom_fruit(offset = 0.2,
             data=tr,
             geom=geom_tile,
             mapping=aes(y=cell,fill=log10(umi)),
             pwidth=0.15,
             axis.params=list(
               axis="x", # add axis text of the layer.
               text.angle=-45, # the text angle of x-axis.
               hjust=0  # adjust the horizontal position of text of axis.
             )
  )+scale_fill_gradientn(colours = colorRampPalette(brewer.pal(9,'GnBu')[1:8])(100))
p4_tree

p5_tree <-  p4_tree + 
  new_scale_fill() +
  geom_fruit(offset = 0.2,
             data=tr,
             geom=geom_tile,
             mapping=aes(y=cell,fill=log10(rna)),
             pwidth=0.15,
             axis.params=list(
               axis="x", # add axis text of the layer.
               text.angle=-45, # the text angle of x-axis.
               hjust=0  # adjust the horizontal position of text of axis.
             )
  )+scale_fill_gradientn(colours = colorRampPalette(brewer.pal(9,'RdPu')[1:8])(100))

p5_tree


p7_tree <-  p5_tree + 
  new_scale_fill() +
  geom_fruit(offset = 0.12,
             data=p2$data,
             geom=geom_bar,
             mapping=aes(x =Freq, y= Var2, fill = Var1), # mapping=aes(x =Var1, y= Freq, fill = Var2),
             orientation = "y",
             stat = 'identity',
             pwidth=0.5,
  )+scale_fill_manual(values =  c("#e20612","#F57F17","#00b0eb"))
p7_tree


p8_tree <-  p7_tree + 
  new_scale_fill() +
  geom_fruit(offset = 0.02,
             data=tr,
             geom=geom_bar,
             mapping=aes(y = cell,x = radio,fill = group), # mapping=aes(x =Var1, y= Freq, fill = Var2),
             stat = 'identity',
             orientation="y",
             pwidth=1.5
  )+scale_fill_manual(
    values=color,
    guide="none" 
  ) 
p8_tree

### 
p8_tree2 <-  p7_tree + 
  new_scale_fill() +
  geom_fruit(offset = 0.02,
             data=tr,
             geom=geom_bar,
             mapping=aes(y = cell,x = -log2(radio),fill = group), # mapping=aes(x =Var1, y= Freq, fill = Var2),
             stat = 'identity',
             orientation="y",
             pwidth=1.5
  )+scale_fill_manual(
    values=color,
    guide="none" 
  ) 
p8_tree2

## 20x20
path <- ""

if(T){
  p8_tree+geom_hilight(node=c(242), fill=color[5], color='white', alpha=0.2, extend=16.5)+# non neu
    geom_hilight(node=c(243), fill=color[3], color='white', alpha=0.2, extend=16.5)+ # chat
    geom_hilight(node=c(259,265,254,256,253,272,268,271,260,248,250), #gaba
                 fill=color[1], color='white', alpha=0.2, extend=16.5)+
    geom_hilight(node=c(246,264), # glu
                 fill=color[2], color='white', alpha=0.2, extend=16.5) -> plot_1
  
  p8_tree2+geom_hilight(node=c(242), fill=color[5], color='white', alpha=0.2, extend=16.5)+# non neu
    geom_hilight(node=c(243), fill=color[3], color='white', alpha=0.2, extend=16.5)+ # chat
    geom_hilight(node=c(259,265,254,256,253,272,268,271,260,248,250), #gaba
                 fill=color[1], color='white', alpha=0.2, extend=16.5)+
    geom_hilight(node=c(246,264), # glu
                 fill=color[2], color='white', alpha=0.2, extend=16.5) -> plot_1_v2
  
  p8_tree+geom_hilight(node=c(242), fill=color[5], color='white', alpha=0.2, extend=16.5)+# non neu
    geom_hilight(node=c(243), fill=color[3], color='white', alpha=0.2, extend=16.5)+ # chat
    geom_hilight(node=c(259,265,254,256,253,272,268,271,260,248,250), #gaba
                 fill=color[1], color='white', alpha=0.2, extend=16.5)+
    geom_hilight(node=c(246,264), # glu
                 fill=color[2], color='white', alpha=0.2, extend=16.5)+
    geom_tiplab(offset = 7.5,color = "white",size =3.5) -> plot_2
  
  
  p8_tree2+geom_hilight(node=c(242), fill=color[5], color='white', alpha=0.2, extend=16.5)+# non neu
    geom_hilight(node=c(243), fill=color[3], color='white', alpha=0.2, extend=16.5)+ # chat
    geom_hilight(node=c(259,265,254,256,253,272,268,271,260,248,250), #gaba
                 fill=color[1], color='white', alpha=0.2, extend=16.5)+
    geom_hilight(node=c(246,264), # glu
                 fill=color[2], color='white', alpha=0.2, extend=16.5)+
    geom_tiplab(offset = 7.5,color = "white",size =3.5) -> plot_2_v2
  
  p8_tree+geom_hilight(node=c(242), fill=color[5], color='white', alpha=0.2, extend=16.5)+# non neu
    geom_hilight(node=c(243), fill=color[3], color='white', alpha=0.2, extend=16.5)+ # chat
    geom_hilight(node=c(259,265,254,256,253,272,268,271,260,248,250), #gaba
                 fill=color[1], color='white', alpha=0.2, extend=16.5)+
    geom_hilight(node=c(246,264), # glu
                 fill=color[2], color='white', alpha=0.2, extend=16.5)+
    geom_tiplab(offset = 7.5,color = "black",size =3.5)-> plot_3
  
  p8_tree2+geom_hilight(node=c(242), fill=color[5], color='white', alpha=0.2, extend=16.5)+# non neu
    geom_hilight(node=c(243), fill=color[3], color='white', alpha=0.2, extend=16.5)+ # chat
    geom_hilight(node=c(259,265,254,256,253,272,268,271,260,248,250), #gaba
                 fill=color[1], color='white', alpha=0.2, extend=16.5)+
    geom_hilight(node=c(246,264), # glu
                 fill=color[2], color='white', alpha=0.2, extend=16.5)+
    geom_tiplab(offset = 7.5,color = "black",size =3.5)-> plot_3_v2
}

ggsave(plot_1,filename = paste0(path,"plot_1.pdf"),width = 20.,height = 20)
ggsave(plot_2,filename = paste0(path,"plot_2.pdf"),width = 20.,height = 20)
ggsave(plot_3,filename = paste0(path,"plot_3.pdf"),width = 20.,height = 20)
ggsave(plot_1_v2,filename = paste0(path,"plot_1_v2.pdf"),width = 20.,height = 20)
ggsave(plot_2_v2,filename = paste0(path,"plot_2_v2.pdf"),width = 20.,height = 20)
ggsave(plot_3_v2,filename = paste0(path,"plot_3_v2.pdf"),width = 20.,height = 20)

save(new_tree,tr,color,p8_tree,p8_tree2,file = "plot_ready.Rdata")
