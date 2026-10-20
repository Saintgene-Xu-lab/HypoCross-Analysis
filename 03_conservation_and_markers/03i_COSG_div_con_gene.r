if(T){
  library(Seurat)
  library(ComplexHeatmap)
  library(ggplot2)
  library(RColorBrewer)
  library(COSG)
  library(qs)
}


############ ---------
neu_sample <- qread("")
neu_sample$anno6 <- factor(neu_sample$anno6,levels = rev(c(paste0("EN",1:30),paste0("IN",1:31),"CHAT")))
Idents(neu_sample) <- neu_sample$anno6
DefaultAssay(neu_sample) <- "RNA"
dat2 <- SplitObject(object = neu_sample,split.by = "specie")
dat2 <- lapply(dat2, function(x){
  x <- NormalizeData(x)
})

mac_marker <- cosg(
  dat2$macaque,
  groups = "all",
  assay = "RNA",
  slot = "data",
  mu = 1,
  n_genes_user = 200)

mar_marker <- cosg(
  dat2$marmoset,
  groups = "all",
  assay = "RNA",
  slot = "data",
  mu = 1,
  n_genes_user = 200)

mou_marker <- cosg(
  dat2$mouse,
  groups = "all",
  assay = "RNA",
  slot = "data",
  mu = 1,
  n_genes_user = 200)

save(mac_marker,mar_marker,mou_marker,file = "cosg_deg.Rdata")
load("gene_category_20230218.Rdata")

cells <- c(paste0("EN",1:30),paste0("IN",1:31),"CHAT")
mou_gene <- list()
mar_gene<- list()
mac_gene<- list()
pri_gene<- list()
common_gene <- list()
for (i in 1:length(cells)) {
  
  ### mouse
  mou_gene[[cells[i]]] <- setdiff(mou_marker$names[1:100,cells[i]],c(mac_marker$names[1:100,cells[i]],
                                                              mar_marker$names[1:100,cells[i]]))
  
  ### marmoset
  mar_gene[[cells[i]]] <- setdiff(mar_marker$names[1:100,cells[i]],c(mac_marker$names[1:100,cells[i]],
                                                              mou_marker$names[1:100,cells[i]]))
  ### macaque
  mac_gene[[cells[i]]] <- setdiff(mac_marker$names[1:100,cells[i]],c(mar_marker$names[1:100,cells[i]],
                                                              mou_marker$names[1:100,cells[i]]))
  
  ### pri
  pri_gene[[cells[i]]] <- setdiff(intersect(mar_marker$names[1:100,cells[i]],
                                          mac_marker$names[1:100,cells[i]]),mou_marker$names[1:100,cells[i]])
  
  #### common 
  common_gene[[cells[i]]] <- intersect(intersect(mar_marker$names[1:100,cells[i]],
                                               mac_marker$names[1:100,cells[i]]),mou_marker$names[1:100,cells[i]])
  
  print(i)
  
}

#### mou

mou_gene <- unique(unlist(mou_gene))
p_mou <- DotPlot(dat2$mouse,features = mou_gene)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mou
p_mac <- DotPlot(dat2$macaque,features = mou_gene)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mac
p_mar <- DotPlot(dat2$marmoset,features = mou_gene)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mar

if(T){
  p_mar$data$gene <- rownames(p_mar$data)
  gene_cell_matrix <- p_mar$data %>%
    select(features.plot, id, avg.exp.scaled) %>%
    pivot_wider(
      names_from = id,        
      values_from = avg.exp.scaled  
    ) 
  gene_cell_matrix <- data.frame(gene_cell_matrix)
  rownames(gene_cell_matrix) <- gene_cell_matrix$features.plot
  a2 <-gene_cell_matrix[,-1]
  
  my_color <- colorRampPalette(c("white","#E1F5FE","#29B6F6"))(50)
  my_breaks <- c(seq(-4, 4, length.out = 51))  
  
  pheatmap(
    gene_cell_matrix[,-1],show_rownames = F,show_colnames = F,legend = F,
    color = my_color,
    breaks = my_breaks,     
    cluster_rows = F,
    cluster_cols = F,
    border_color = NA
  )
  
  ########
  p_mac$data$gene <- rownames(p_mac$data)
  gene_cell_matrix <- p_mac$data %>%
    select(features.plot, id, avg.exp.scaled) %>%
    pivot_wider(
      names_from = id,        
      values_from = avg.exp.scaled  
    ) 
  gene_cell_matrix <- data.frame(gene_cell_matrix)
  rownames(gene_cell_matrix) <- gene_cell_matrix$features.plot
  a2 <-gene_cell_matrix[,-1]
  
  my_color <- colorRampPalette(c("white","#E1F5FE","#29B6F6"))(50)
  my_breaks <- c(seq(-4, 4, length.out = 51))  
  
  pheatmap(
    gene_cell_matrix[,-1],show_rownames = F,show_colnames = F,legend = F,
    color = my_color,
    breaks = my_breaks,     
    cluster_rows = F,
    cluster_cols = F,
    border_color = NA
  )
  ##########
  ########
  p_mou$data$gene <- rownames(p_mou$data)
  gene_cell_matrix <- p_mou$data %>%
    select(features.plot, id, avg.exp.scaled) %>%
    pivot_wider(
      names_from = id,        
      values_from = avg.exp.scaled  
    ) 
  gene_cell_matrix <- data.frame(gene_cell_matrix)
  rownames(gene_cell_matrix) <- gene_cell_matrix$features.plot
  a2 <-gene_cell_matrix[,-1]
  
  my_color <- colorRampPalette(c("white","#E1F5FE","#29B6F6"))(50)
  my_breaks <- c(seq(-4, 4, length.out = 51))  
  
  pheatmap(
    gene_cell_matrix[,-1],show_rownames = F,show_colnames = F,legend = F,
    color = my_color,
    breaks = my_breaks,
    cluster_rows = F,
    cluster_cols = F,
    border_color = NA
  )
  
}

#### mar
mar_gene <- unique(unlist(mar_gene))
p_mou <- DotPlot(dat2$mouse,features = mar_gene)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mou
p_mac <- DotPlot(dat2$macaque,features = mar_gene)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mac
p_mar <- DotPlot(dat2$marmoset,features = mar_gene)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mar

if(T){
  p_mar$data$gene <- rownames(p_mar$data)
  gene_cell_matrix <- p_mar$data %>%
    select(features.plot, id, avg.exp.scaled) %>%
    pivot_wider(
      names_from = id,        
      values_from = avg.exp.scaled  
    ) 
  gene_cell_matrix <- data.frame(gene_cell_matrix)
  rownames(gene_cell_matrix) <- gene_cell_matrix$features.plot
  a2 <-gene_cell_matrix[,-1]
  
  my_color <- colorRampPalette(c("white", "#FFF3E0","#FF5722"))(50)
  my_breaks <- c(seq(-4, 4, length.out = 51))  
  
  pheatmap(
    gene_cell_matrix[,-1],
    color = my_color,
    breaks = my_breaks,      
    cluster_rows = F,
    cluster_cols = F,
    border_color = NA
  )
  
  ########
  p_mac$data$gene <- rownames(p_mac$data)
  gene_cell_matrix <- p_mac$data %>%
    select(features.plot, id, avg.exp.scaled) %>%
    pivot_wider(
      names_from = id,        
      values_from = avg.exp.scaled  
    ) 
  gene_cell_matrix <- data.frame(gene_cell_matrix)
  rownames(gene_cell_matrix) <- gene_cell_matrix$features.plot
  a2 <-gene_cell_matrix[,-1]
  
  my_color <- colorRampPalette(c("white", "#FFF3E0","#FF5722"))(50)
  my_breaks <- c(seq(-4, 4, length.out = 51))  
  
  pheatmap(
    gene_cell_matrix[,-1],
    color = my_color,
    breaks = my_breaks,      
    cluster_rows = F,
    cluster_cols = F,
    border_color = NA
  )
  ##########
  ########
  p_mou$data$gene <- rownames(p_mou$data)
  gene_cell_matrix <- p_mou$data %>%
    select(features.plot, id, avg.exp.scaled) %>%
    pivot_wider(
      names_from = id,        
      values_from = avg.exp.scaled  
    ) 
  gene_cell_matrix <- data.frame(gene_cell_matrix)
  rownames(gene_cell_matrix) <- gene_cell_matrix$features.plot
  a2 <-gene_cell_matrix[,-1]
  
  my_color <- colorRampPalette(c("white", "#FFF3E0","#FF5722"))(50)
  my_breaks <- c(seq(-4, 4, length.out = 51))  
  
  pheatmap(
    gene_cell_matrix[,-1],
    color = my_color,
    breaks = my_breaks,
    cluster_rows = F,
    cluster_cols = F,
    border_color = NA
  )
  
}

#### mac
mac_gene <- unique(unlist(mac_gene))
p_mou <- DotPlot(dat2$mouse,features = mac_gene)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mou
p_mac <- DotPlot(dat2$macaque,features = mac_gene)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mac
p_mar <- DotPlot(dat2$marmoset,features = mac_gene)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mar

if(T){
  p_mar$data$gene <- rownames(p_mar$data)
  gene_cell_matrix <- p_mar$data %>%
    select(features.plot, id, avg.exp.scaled) %>%
    pivot_wider(
      names_from = id,        
      values_from = avg.exp.scaled  
    ) 
  gene_cell_matrix <- data.frame(gene_cell_matrix)
  rownames(gene_cell_matrix) <- gene_cell_matrix$features.plot
  a2 <-gene_cell_matrix[,-1]
  
  my_color <- colorRampPalette(c("white", "#FFEBEE","#F44336"))(50)
  my_breaks <- c(seq(-4, 4, length.out = 51))  
  
  pheatmap(
    gene_cell_matrix[,-1],
    color = my_color,
    breaks = my_breaks,      
    cluster_rows = F,
    cluster_cols = F,
    border_color = NA
  )
  
  ########
  p_mac$data$gene <- rownames(p_mac$data)
  gene_cell_matrix <- p_mac$data %>%
    select(features.plot, id, avg.exp.scaled) %>%
    pivot_wider(
      names_from = id,        
      values_from = avg.exp.scaled  
    ) 
  gene_cell_matrix <- data.frame(gene_cell_matrix)
  rownames(gene_cell_matrix) <- gene_cell_matrix$features.plot
  a2 <-gene_cell_matrix[,-1]
  
  my_color <- colorRampPalette(c("white", "#FFEBEE","#F44336"))(50)
  my_breaks <- c(seq(-4, 4, length.out = 51))  
  
  pheatmap(
    gene_cell_matrix[,-1],
    color = my_color,
    breaks = my_breaks,      
    cluster_rows = F,
    cluster_cols = F,
    border_color = NA
  )
  ##########
  ########
  p_mou$data$gene <- rownames(p_mou$data)
  gene_cell_matrix <- p_mou$data %>%
    select(features.plot, id, avg.exp.scaled) %>%
    pivot_wider(
      names_from = id,        
      values_from = avg.exp.scaled  
    ) 
  gene_cell_matrix <- data.frame(gene_cell_matrix)
  rownames(gene_cell_matrix) <- gene_cell_matrix$features.plot
  a2 <-gene_cell_matrix[,-1]
  
  my_color <- colorRampPalette(c("white", "#FFEBEE","#F44336"))(50)
  my_breaks <- c(seq(-4, 4, length.out = 51))  
  
  pheatmap(
    gene_cell_matrix[,-1],
    color = my_color,
    breaks = my_breaks,
    cluster_rows = F,
    cluster_cols = F,
    border_color = NA
  )
  
}

#### pri 
pri_gene <- unique(unlist(pri_gene))
p_mou <- DotPlot(dat2$mouse,features = pri_gene)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mou
p_mac <- DotPlot(dat2$macaque,features = pri_gene)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mac
p_mar <- DotPlot(dat2$marmoset,features = pri_gene)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mar

if(T){
  p_mar$data$gene <- rownames(p_mar$data)
  gene_cell_matrix <- p_mar$data %>%
    select(features.plot, id, avg.exp.scaled) %>%
    pivot_wider(
      names_from = id,        
      values_from = avg.exp.scaled  
    ) 
  gene_cell_matrix <- data.frame(gene_cell_matrix)
  rownames(gene_cell_matrix) <- gene_cell_matrix$features.plot
  a2 <-gene_cell_matrix[,-1]
  
  my_color <- colorRampPalette(c("white", "#E8EAF6","#3F51B5"))(50)
  my_breaks <- c(seq(-4, 4, length.out = 51))  
  
  pheatmap(
    gene_cell_matrix[,-1],
    color = my_color,
    breaks = my_breaks,      
    cluster_rows = F,
    cluster_cols = F,
    border_color = NA
  )
  
  ########
  p_mac$data$gene <- rownames(p_mac$data)
  gene_cell_matrix <- p_mac$data %>%
    select(features.plot, id, avg.exp.scaled) %>%
    pivot_wider(
      names_from = id,        
      values_from = avg.exp.scaled  
    ) 
  gene_cell_matrix <- data.frame(gene_cell_matrix)
  rownames(gene_cell_matrix) <- gene_cell_matrix$features.plot
  a2 <-gene_cell_matrix[,-1]
  
  my_color <- colorRampPalette(c("white", "#E8EAF6","#3F51B5"))(50)
  my_breaks <- c(seq(-4, 4, length.out = 51))  
  
  pheatmap(
    gene_cell_matrix[,-1],
    color = my_color,
    breaks = my_breaks,      
    cluster_rows = F,
    cluster_cols = F,
    border_color = NA
  )
  ##########
  ########
  p_mou$data$gene <- rownames(p_mou$data)
  gene_cell_matrix <- p_mou$data %>%
    select(features.plot, id, avg.exp.scaled) %>%
    pivot_wider(
      names_from = id,        
      values_from = avg.exp.scaled  
    ) 
  gene_cell_matrix <- data.frame(gene_cell_matrix)
  rownames(gene_cell_matrix) <- gene_cell_matrix$features.plot
  a2 <-gene_cell_matrix[,-1]
  
  my_color <- colorRampPalette(c("white", "#E8EAF6","#3F51B5"))(50)
  my_breaks <- c(seq(-4, 4, length.out = 51))  
  
  pheatmap(
    gene_cell_matrix[,-1],
    color = my_color,
    breaks = my_breaks,
    cluster_rows = F,
    cluster_cols = F,
    border_color = NA
  )
  
}


#### common 
common_gene <- unique(unlist(common_gene))
p_mou <- DotPlot(dat2$mouse,features = common_gene)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mou
p_mac <- DotPlot(dat2$macaque,features = common_gene)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mac
p_mar <- DotPlot(dat2$marmoset,features = common_gene)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mar

if(T){
  p_mar$data$gene <- rownames(p_mar$data)
  gene_cell_matrix <- p_mar$data %>%
    select(features.plot, id, avg.exp.scaled) %>%
    pivot_wider(
      names_from = id,        
      values_from = avg.exp.scaled  
    ) 
  gene_cell_matrix <- data.frame(gene_cell_matrix)
  rownames(gene_cell_matrix) <- gene_cell_matrix$features.plot
  a2 <-gene_cell_matrix[,-1]
  
  my_color <- colorRampPalette(c("white", "#EDE7F6","#673AB7"))(50)
  my_breaks <- c(seq(-4, 4, length.out = 51))  
  
  pheatmap(
    gene_cell_matrix[,-1],
    color = my_color,
    breaks = my_breaks,      
    cluster_rows = F,
    cluster_cols = F,
    border_color = NA
  )
  
  ########
  p_mac$data$gene <- rownames(p_mac$data)
  gene_cell_matrix <- p_mac$data %>%
    select(features.plot, id, avg.exp.scaled) %>%
    pivot_wider(
      names_from = id,        
      values_from = avg.exp.scaled  
    ) 
  gene_cell_matrix <- data.frame(gene_cell_matrix)
  rownames(gene_cell_matrix) <- gene_cell_matrix$features.plot
  a2 <-gene_cell_matrix[,-1]
  
  my_color <- colorRampPalette(c("white", "#EDE7F6","#673AB7"))(50)
  my_breaks <- c(seq(-4, 4, length.out = 51))  # 强制 -4 ~ 4
  
  pheatmap(
    gene_cell_matrix[,-1],
    color = my_color,
    breaks = my_breaks,     
    cluster_rows = F,
    cluster_cols = F,
    border_color = NA
  )
  ##########
  ########
  p_mou$data$gene <- rownames(p_mou$data)
  gene_cell_matrix <- p_mou$data %>%
    select(features.plot, id, avg.exp.scaled) %>%
    pivot_wider(
      names_from = id,        
      values_from = avg.exp.scaled  
    ) 
  gene_cell_matrix <- data.frame(gene_cell_matrix)
  rownames(gene_cell_matrix) <- gene_cell_matrix$features.plot
  a2 <-gene_cell_matrix[,-1]
  
  my_color <- colorRampPalette(c("white", "#EDE7F6","#673AB7"))(50)
  my_breaks <- c(seq(-4, 4, length.out = 51))  
  
  pheatmap(
    gene_cell_matrix[,-1],
    color = my_color,
    breaks = my_breaks,
    cluster_rows = F,
    cluster_cols = F,
    border_color = NA
  )
  
}






