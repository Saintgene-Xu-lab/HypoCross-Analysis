
if(T){
  library(Seurat)
  library(qs)
  library(reshape2)
  library(RColorBrewer)
  library(scales)
  library(tidyverse)
  library(stringr)
}

mac <- qread("")
mar <- qread("")
mou <- qread("")
############ -----

##### mac --------------------------------------
meta3 <- mac@meta.data
#meta3 <- meta3[!meta3$Assign_Simple2 == "HY",]
meta3 <- meta3[!meta3$banksy_clu6_202512 == "HY",]
meta3 <- meta3[!meta3$Assign_Simple2 == "HY",]
ann1 <- names(table(meta3$banksy_clu6_202512))
ann2 <- names(table(meta3$Assign_Simple2))
ann1_result <- list()
ann2_result<- list()
for (i in 1:length(ann1)){
  ann2_result<- list()
  for (j in 1:length(ann2)) {
    overlap_result <- intersect(rownames(meta3[meta3$banksy_clu6_202512 %in% ann1[i],]),rownames(meta3[meta3$Assign_Simple2 %in% ann2[j],]))
    union_result <- unique(c(rownames(meta3[meta3$banksy_clu6_202512 %in% ann1[i],]),rownames(meta3[meta3$Assign_Simple2 %in% ann2[j],])))
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
                   color = colorRampPalette(c("#FAFAFA","#E1F5FE","#81D4FA","#29B6F6","#039BE5"))(100))


###### mar --------------------------------
meta3 <- mar@meta.data
meta3 <- meta3[!meta3$Assign_Simple2 == "HY",]
meta3 <- meta3[!meta3$banksy_clu6_202512 == "HY",]
ann1 <- names(table(meta3$banksy_clu6_202512))
ann2 <- names(table(meta3$Assign_Simple2))
ann1_result <- list()
ann2_result<- list()
for (i in 1:length(ann1)){
  ann2_result<- list()
  for (j in 1:length(ann2)) {
    overlap_result <- intersect(rownames(meta3[meta3$banksy_clu6_202512 %in% ann1[i],]),rownames(meta3[meta3$Assign_Simple2 %in% ann2[j],]))
    union_result <- unique(c(rownames(meta3[meta3$banksy_clu6_202512 %in% ann1[i],]),rownames(meta3[meta3$Assign_Simple2 %in% ann2[j],])))
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
                   color = colorRampPalette(c("#FAFAFA","#E1F5FE","#81D4FA","#29B6F6","#039BE5"))(100))

###### mar --------------------------------
meta3 <- mou@meta.data
meta3 <- meta3[!meta3$Assign_Simple2 == "HY",]
meta3 <- meta3[!meta3$banksy_clu6_20251202 == "HY",]
ann1 <- names(table(meta3$banksy_clu6_20251202))
ann2 <- names(table(meta3$Assign_Simple2))
ann1_result <- list()
ann2_result<- list()
for (i in 1:length(ann1)){
  ann2_result<- list()
  for (j in 1:length(ann2)) {
    overlap_result <- intersect(rownames(meta3[meta3$banksy_clu6_20251202 %in% ann1[i],]),rownames(meta3[meta3$Assign_Simple2 %in% ann2[j],]))
    union_result <- unique(c(rownames(meta3[meta3$banksy_clu6_20251202 %in% ann1[i],]),rownames(meta3[meta3$Assign_Simple2 %in% ann2[j],])))
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
                   color = colorRampPalette(c("#FAFAFA","#E1F5FE","#81D4FA","#29B6F6","#039BE5"))(100))

add_col <- data.frame(rep(0,17))
colnames(add_col) <- "MarMG"
rownames(add_col) <- rownames(mou_result)
mou_result2 <- cbind(mou_result,add_col)
mac_result2 <- cbind(mac_result,add_col)
mar_result2 <- mar_result

nucle2 <- c("ARH","VMH","PVH","SO","TU","SCH","MBO","Pe","ME","BM","ZI","PM","PH","POA","DMH","AHN","LHA","MarMG")
mou_result3 <- ranking_mat(mou_result2[,nucle2])
mar_result3 <- ranking_mat(mar_result2[,nucle2])
mac_result3 <- ranking_mat(mac_result2[,nucle2])

colnames(mou_result3 ) <- paste0("st:",colnames(mou_result3))
colnames(mar_result3 ) <- paste0("st:",colnames(mar_result3))
colnames(mac_result3 ) <- paste0("st:",colnames(mac_result3))

pheatmap::pheatmap(t(mou_result3),cluster_rows = F,breaks = seq(from = 0, to = 0.5, length.out = 100),
                   cluster_cols = F,border_color = "white",filename = paste0(path,"mou_result.pdf"),
                   color = colorRampPalette(c("#FAFAFA","#E1F5FE","#81D4FA","#29B6F6","#039BE5"))(100))
pheatmap::pheatmap(t(mar_result3),cluster_rows = F,breaks = seq(from = 0, to = 0.5, length.out = 100),
                   cluster_cols = F,border_color = "white",filename = paste0(path,"mar_result.pdf"),
                   color = colorRampPalette(c("#FAFAFA","#E1F5FE","#81D4FA","#29B6F6","#039BE5"))(100))
pheatmap::pheatmap(t(mac_result3),cluster_rows = F,breaks = seq(from = 0, to = 0.5, length.out = 100),
                   cluster_cols = F,border_color = "white",filename = paste0(path,"mac_result.pdf"),
                   color = colorRampPalette(c("#FAFAFA","#E1F5FE","#81D4FA","#29B6F6","#039BE5"))(100))




