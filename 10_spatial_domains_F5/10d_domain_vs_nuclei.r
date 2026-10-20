if(T){
  library(Seurat)
  library(ggplot2)
  library(tidyverse)
  library(reshape2)
  library(qs)
  library(RColorBrewer)
  library(cowplot)
  library(viridis)
  display.brewer.all()
}

mou <- qread("mou.qs")
mar <- qread("mar.qs")
mac <- qread("mac.qs")

##### mac --------------------------------------
meta3 <- mac@meta.data
meta3 <- meta3[!meta3$domian_modif2 == "HY",]
meta3 <- meta3[!meta3$Assign_Simple2 == "HY",]
ann1 <- names(table(meta3$domian_modif2))
ann2 <- names(table(meta3$Assign_Simple2))
ann1_result <- list()
ann2_result<- list()
for (i in 1:length(ann1)){
  ann2_result<- list()
  for (j in 1:length(ann2)) {
    overlap_result <- intersect(rownames(meta3[meta3$domian_modif2 %in% ann1[i],]),rownames(meta3[meta3$Assign_Simple2 %in% ann2[j],]))
    union_result <- unique(c(rownames(meta3[meta3$domian_modif2clu6_202512 %in% ann1[i],]),rownames(meta3[meta3$Assign_Simple2 %in% ann2[j],])))
    ann2_result[[ann2[j]]] <- length(overlap_result)/length(union_result)
  }
  ann1_result[[ann1[i]]] <- as.numeric(ann2_result)
}

result <- data.frame(ann1_result)
rownames(result ) <- ann2
mac_result <- result
nucle <- c("POA","AHN","SO","LHA","fiber","ZI","PeVH","SCN","PVH","ARH","VMH","DMH","TU","CVO","PH","PM","MBO","SUM")
nucle2 <- nucle[-c(5)]
mac_result <- mac_result[nucle2,]
pheatmap::pheatmap(mac_result,cluster_rows = F,
                   cluster_cols = F,border_color = "white",
                   color = colorRampPalette(c("#FAFAFA",rev(viridis_pal(option = "G")(10))[1:5]))(50))

###### mar --------------------------------
meta3 <- mar@meta.data
meta3 <- meta3[!meta3$Assign_Simple2 == "HY",]
meta3 <- meta3[!meta3$domian_modif2 == "HY",]
ann1 <- names(table(meta3$domian_modif2))
ann2 <- names(table(meta3$Assign_Simple2))
ann1_result <- list()
ann2_result<- list()
for (i in 1:length(ann1)){
  ann2_result<- list()
  for (j in 1:length(ann2)) {
    overlap_result <- intersect(rownames(meta3[meta3$domian_modif2 %in% ann1[i],]),rownames(meta3[meta3$Assign_Simple2 %in% ann2[j],]))
    union_result <- unique(c(rownames(meta3[meta3$domian_modif2 %in% ann1[i],]),rownames(meta3[meta3$Assign_Simple2 %in% ann2[j],])))
    ann2_result[[ann2[j]]] <- length(overlap_result)/length(union_result)
  }
  ann1_result[[ann1[i]]] <- as.numeric(ann2_result)
}

result <- data.frame(ann1_result)
rownames(result ) <- ann2
mar_result <- result
nucle <- c("POA","AHN","SO","LHA","ZI","PeVH","SCN","PVH","ARH","VMH","DMH","TU","PH","PM","MBO","SUM")
mar_result <- mar_result[nucle,]
pheatmap::pheatmap(mar_result,cluster_rows = F,
                   cluster_cols = F,border_color = "white",
                   color = colorRampPalette(c("#FAFAFA",rev(viridis_pal(option = "G")(10))[1:5]))(50))

###### mou --------------------------------
meta3 <- mou@meta.data
meta3 <- meta3[!meta3$Assign_Simple2 == "HY",]
meta3 <- meta3[!meta3$domian_modif2 == "HY",]
ann1 <- names(table(meta3$domian_modif2))
ann2 <- names(table(meta3$Assign_Simple2))
ann1_result <- list()
ann2_result<- list()
for (i in 1:length(ann1)){
  ann2_result<- list()
  for (j in 1:length(ann2)) {
    overlap_result <- intersect(rownames(meta3[meta3$domian_modif2 %in% ann1[i],]),rownames(meta3[meta3$Assign_Simple2 %in% ann2[j],]))
    union_result <- unique(c(rownames(meta3[meta3$domian_modif2 %in% ann1[i],]),rownames(meta3[meta3$Assign_Simple2 %in% ann2[j],])))
    ann2_result[[ann2[j]]] <- length(overlap_result)/length(union_result)
  }
  ann1_result[[ann1[i]]] <- as.numeric(ann2_result)
}
result <- data.frame(ann1_result)
rownames(result ) <- ann2
mou_result <- result
nucle <- c("POA","AHN","SO","LHA","ZI","PeVH","SCN","PVH","ARH","VMH","DMH","TU","CVO","PH","PM","MBO","SUM")
mou_result <- mou_result[nucle,]
pheatmap::pheatmap(mou_result,cluster_rows = F,
                   cluster_cols = F,border_color = "white",
                   color = colorRampPalette(c("#FAFAFA",rev(viridis_pal(option = "G")(10))[1:5]))(50))


library("ComplexHeatmap")
library("zoo")
ranking_mat <- function(data){
  library(zoo)
  window <- 1
  step <- 1
  dat_ts <- data.frame(data)
  dat2_ts <- dat_ts[order(apply(t(rollapply(t(dat_ts), width=window, by=step, FUN=mean)), 1, which.max)), ]
  return(dat2_ts)
}

nucle <- c("ARH","VMH","PVH","SO","TU","SCN","MBO","SUM","PeVH","CVO","ZI","PM","PH","POA","DMH","AHN","LHA")
nucle2 <- intersect(nucle,rownames(mar_result))

domain <- paste0("D",c(7,10,18,19,22,17,15,12,20,13,16,2,4,14,6,3,5,11,8,1,9,21))

mou_result3 <- ranking_mat(mou_result[,domain])
mar_result3 <- ranking_mat(mar_result[,domain])
mac_result3 <- ranking_mat(mac_result[,domain])


path <- ""
pheatmap::pheatmap(t(mou_result3),cluster_rows = F,breaks = seq(from = 0, to = 0.5, length.out = 100),
                   cluster_cols = F,border_color = "white",angle_col = c("90"),filename = paste0(path,"mou_t.pdf"),
                   color = colorRampPalette(c("#FAFAFA",rev(viridis_pal(option = "G")(10))[1:5]))(100))
pheatmap::pheatmap(t(mar_result3),cluster_rows = F,breaks = seq(from = 0, to = 0.5, length.out = 100),
                   cluster_cols = F,border_color = "white",angle_col = c("90"),filename = paste0(path,"mar_t.pdf"),
                   color = colorRampPalette(c("#FAFAFA",rev(viridis_pal(option = "G")(10))[1:5]))(100))
pheatmap::pheatmap(t(mac_result3),cluster_rows = F,breaks = seq(from = 0, to = 0.5, length.out = 100),
                   cluster_cols = F,border_color = "white",angle_col = c("90"),filename = paste0(path,"mac_t.pdf"),
                   color = colorRampPalette(c("#FAFAFA",rev(viridis_pal(option = "G")(10))[1:5]))(100))

save(mac_result3,mar_result3,mou_result3,file ="plot.Rdata" )

pheatmap::pheatmap(t(mac_result3),cluster_rows = F,breaks = seq(from = 0, to = 0.5, length.out = 100),
                   cluster_cols = F,border_color = "white",angle_col = c("90"),
                   color = colorRampPalette(c("#FAFAFA",rev(viridis_pal(option = "G")(10))[1:5]))(100))





