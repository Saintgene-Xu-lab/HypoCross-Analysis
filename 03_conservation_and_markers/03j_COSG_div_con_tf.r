if(T){
  library(Seurat)
  library(ComplexHeatmap)
  library(ggplot2)
  library(RColorBrewer)
  library(COSG)
  library(qs)
}



neu_sample <- qread()
neu_sample$anno6 <- factor(neu_sample$anno6,levels = rev( c(paste0("EN",1:30),paste0("IN",1:31),"CHAT")))
Idents(neu_sample) <- neu_sample$anno6
DefaultAssay(neu_sample) <- "RNA"
dat2 <- SplitObject(object = neu_sample,split.by = "specie")
dat2 <- lapply(dat2, function(x){
  x <- NormalizeData(x)
})


load("cosg_deg.Rdata")
load("gene_category_20230218.Rdata")

cells <- c(paste0("EN",1:30),paste0("IN",1:31),"CHAT")
mou_tf <- list()
mar_tf<- list()
mac_tf<- list()
pri_tf<- list()
common_tf <- list()
for (i in 1:length(cells)) {
  
  ### mouse
  mou_tf[[cells[i]]] <- setdiff(mou_marker$names[,cells[i]],c(mac_marker$names[,cells[i]],
                                                                                mar_marker$names[,cells[i]]))
  
  mou_tf[[cells[i]]] <- intersect(mou_tf[[cells[i]]],gene_category$tf$mouse_tf_list)
  
  ### marmoset
  mar_tf[[cells[i]]] <- setdiff(mar_marker$names[,cells[i]],c(mac_marker$names[,cells[i]],
                                                                     mou_marker$names[,cells[i]]))
  
  mar_tf[[cells[i]]] <- intersect(mar_tf[[cells[i]]],gene_category$tf$mouse_tf_list)
  
  
  ### macaque
  mac_tf[[cells[i]]] <- setdiff(mac_marker$names[,cells[i]],c(mar_marker$names[,cells[i]],
                                                                     mou_marker$names[,cells[i]]))
  
  mac_tf[[cells[i]]] <- intersect(mac_tf[[cells[i]]],gene_category$tf$mouse_tf_list)
  
  
  ### pri
  pri_tf[[cells[i]]] <- setdiff(intersect(mar_marker$names[,cells[i]],
                                         mac_marker$names[,cells[i]]),mou_marker$names[,cells[i]])
  
  pri_tf[[cells[i]]] <- intersect(pri_tf[[cells[i]]],gene_category$tf$mouse_tf_list)
  
  #### common 
  common_tf[[cells[i]]] <- intersect(intersect(mar_marker$names[,cells[i]],
                                             mac_marker$names[,cells[i]]),mou_marker$names[,cells[i]])
  
  common_tf[[cells[i]]] <- intersect(common_tf[[cells[i]]],gene_category$tf$mouse_tf_list)
  
  print(i)
  
}

#### mou

mou_tf <- unique(unlist(mou_tf))
p_mou <- DotPlot(dat2$mouse,features = mou_tf)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mou
p_mac <- DotPlot(dat2$macaque,features = mou_tf)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mac
p_mar <- DotPlot(dat2$marmoset,features = mou_tf)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
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
  
  my_color <- colorRampPalette(c("white","#E1F5FE","#29B6F6"))(50)
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
  
  my_color <- colorRampPalette(c("white","#E1F5FE","#29B6F6"))(50)
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

#### mar
mar_tf <- unique(unlist(mar_tf))
p_mou <- DotPlot(dat2$mouse,features = mar_tf)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mou
p_mac <- DotPlot(dat2$macaque,features = mar_tf)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mac
p_mar <- DotPlot(dat2$marmoset,features = mar_tf)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
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
mac_tf <- unique(unlist(mac_tf))
p_mou <- DotPlot(dat2$mouse,features = mac_tf)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mou
p_mac <- DotPlot(dat2$macaque,features = mac_tf)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mac
p_mar <- DotPlot(dat2$marmoset,features = mac_tf)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
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
pri_tf <- unique(unlist(pri_tf))
p_mou <- DotPlot(dat2$mouse,features = pri_tf)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mou
p_mac <- DotPlot(dat2$macaque,features = pri_tf)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mac
p_mar <- DotPlot(dat2$marmoset,features = pri_tf)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
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
common_tf <- unique(unlist(common_tf))
p_mou <- DotPlot(dat2$mouse,features = common_tf)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mou
p_mac <- DotPlot(dat2$macaque,features = common_tf)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p_mac
p_mar <- DotPlot(dat2$marmoset,features = common_tf)+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
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







