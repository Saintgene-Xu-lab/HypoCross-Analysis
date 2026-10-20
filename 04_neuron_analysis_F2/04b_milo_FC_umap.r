tbUMAP <- readRDS("assign_umap_and_cell.rds")
szDot = 1.5

library(ggplot2)

Rg = round(range(-tbUMAP$logFC));
tbUMAP_EN7 <- tbUMAP[tbUMAP$celltype == 'EN7', ]

logFC <- tbUMAP$logFC
# 
cellType <- unique(tbUMAP$celltype)
nCC <- length(cellType)

# 
for (nC in 1:nCC) {
  lgSel <- tbUMAP$celltype == cellType[nC]
  logFC_Sel <- logFC[lgSel]
  CV_r <- abs(mean(logFC_Sel) / sd(logFC_Sel))
  tbUMAP$logFC[lgSel] <- logFC_Sel * CV_r
}

#
#
tbUMAP_other <- tbUMAP[tbUMAP$celltype != 'EN7', ]
tbUMAP_other <- tbUMAP_other[order(tbUMAP_other$logFC),]
tbUMAP_EN7 <- tbUMAP[tbUMAP$celltype == 'EN7', ]
tbUMAP_EN7 <- tbUMAP_EN7[order(-tbUMAP_EN7$logFC),]
Rg = round(range(-tbUMAP$logFC));

ggplot() +
  geom_point(data = tbUMAP_other, aes(x = x, y = y, color = -logFC), size = szDot, alpha = 0.8) +
  geom_point(data = tbUMAP_EN7, aes(x = x, y = y, color = -logFC), size = szDot, alpha = 0.8) +
  scale_color_gradientn(colors = c("#4575b4", "#f0f0f0", "#d73027"), values = scales::rescale(c(Rg[1], 0, Rg[2])), limits = c(Rg[1], Rg[2]))+
  annotate("text", x = mean(tbUMAP_EN7$x), y = mean(tbUMAP_EN7$y), label = "EN7", hjust = -0.2, vjust = -0.2, color = "black", size = 5) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
  guides(color = guide_colorbar(title = "Log FC"))

library(cowplot)
szDot = 0.5

ggplot() +
  geom_point(data = tbUMAP_other, aes(x = x, y = y, color = -logFC), size = szDot, alpha = 0.8) +
  geom_point(data = tbUMAP_EN7, aes(x = x, y = y, color = -logFC), size = szDot, alpha = 0.8) +
  scale_color_gradientn(colors = c("#4575b4", "#f0f0f0", "#d73027"), values = scales::rescale(c(Rg[1], 0, Rg[2])), limits = c(Rg[1], Rg[2]))+
  annotate("text", x = mean(tbUMAP_EN7$x), y = mean(tbUMAP_EN7$y), label = "EN7", hjust = -0.2, vjust = -0.2, color = "black", size = 5) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
  guides(color = guide_colorbar(title = "Log FC")) -> p 
p
ggsave(p,filename = "milor_umap2.pdf",width = 4,height = 4)

tbUMAP2 <- rbind(tbUMAP_EN7,tbUMAP_other)

tbUMAP2
add_EN30 <- data.frame(matrix(ncol = 7,nrow = 3,data = 0))
colnames(add_EN30) <- colnames(tbUMAP2)
add_EN30$celltype <- "EN30" ### recover cell during filtering
add_EN30$logFC <- 0
add_EN30 <- rbind(tbUMAP2,add_EN30)
add_EN30$celltype <- factor(add_EN30$celltype,levels = c(paste0("EN",1:30),paste0("IN",1:31),"CHAT"))
anno <- data.frame(table(add_EN30$celltype))
anno$anno <- rownames(anno)
anno$seq <- 1:nrow(anno)
cells <- anno$Var1
add_EN30$X <- "NA"
for (i in 1:length(cells)) {
  add_EN30[add_EN30$celltype %in% cells[i],]$X <- anno[anno$Var1 %in% cells[i],]$seq
}
add_EN30 <- add_EN30[order(add_EN30$X),]
add_EN30$X <- factor(add_EN30$X,levels = 1:62)


Rg = round(range(-add_EN30$logFC));
ggplot(data = add_EN30,aes(x = X,y = -logFC, color = -logFC))+
  geom_quasirandom(alpha = 1) +
  theme_bw(base_size = 22) + 
  scale_color_gradientn(colors = c("#4575b4", "#f0f0f0", "#d73027"))+
  scale_x_discrete(labels = c(paste0("EN",1:30),paste0("IN",1:31),"CHAT"))+
  theme(strip.text.y = element_text(angle = 0),axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))->p

p <- p/p

path <- ""
ggsave(p,filename = paste0(path,"dot_milor2.pdf"),width = 20,height = 5)

