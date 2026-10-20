if(T){
  library(Seurat)
  library(ggplot2)
  library(tidyverse)
  library(reshape2)
  library(stringr)
  library(qs)
}

###### mac1
int <- qread("")

############### -----
meta <- int@meta.data
meta$major_anno <- "aa"
meta$near_anno <- "aa"

mac <- read.csv("macaque_bin100_majority_251202.csv")
mac$id <- paste0("mac_",mac$Unnamed..0)
rownames(mac) <- mac$id
mac <- mac[rownames(meta[meta$specie %in% "macaque",]),]
mac$clu6 <- str_remove(mac$clu6,"[0-9]")

mar <- read.csv("marmoset_bin100_majority_251202.csv")
mar$id <- paste0("mar_",mar$Unnamed..0)
rownames(mar) <- mar$id
mar <- mar[rownames(meta[meta$specie %in% "marmoset",]),]
mar$clu6 <- str_remove(mar$clu6,"[0-9]")

mou <- read.csv("mouse_clu6_bin100_majority_251202.csv")
mou$id <- paste0("mou_",mou$name)
rownames(mou) <- mou$id
mou <- mou[rownames(meta[meta$specie %in% "mouse",]),]
mou$clu6 <- str_remove(mou$clu6,"[0-9]")

meta[meta$specie %in% "macaque",]$major_anno <- mac$clu6
meta[meta$specie %in% "marmoset",]$major_anno <- mar$clu6
meta[meta$specie %in% "mouse",]$major_anno <- mou$clu6

int$bansky_202501202 <- meta$major_anno

##########
mou <- qread("")
mar <- qread("")
mac <- qread("")

mou_meta <- mou@meta.data
rownames(mou_meta) <- paste0("mou_",rownames(mou_meta))
mar_meta <- mar@meta.data
rownames(mar_meta) <- paste0("mar_",rownames(mar_meta))
mac_meta <- mac@meta.data
rownames(mac_meta) <- paste0("mac_",rownames(mac_meta))

meta <- int@meta.data
meta$domian_modif2 <- "HY"
meta[intersect(rownames(meta),rownames(mou_meta)),]$domian_modif2  <- mou_meta[intersect(rownames(meta),rownames(mou_meta)),]$domian_modif2
meta[intersect(rownames(meta),rownames(mar_meta)),]$domian_modif2  <- mar_meta[intersect(rownames(meta),rownames(mar_meta)),]$domian_modif2
meta[intersect(rownames(meta),rownames(mac_meta)),]$domian_modif2  <- mac_meta[intersect(rownames(meta),rownames(mac_meta)),]$domian_modif2

int$domian_modif2 <- meta$domian_modif2


id <- names(table(meta$bansky_202501202))[-1]
meta2 <- meta[meta$bansky_202501202 %in% id,]
meta2 <- meta2[!meta2$bansky_202501202 == "HY",]
meta2 <- meta2[!meta2$domian_modif2 == "HY",]


mou$bansky_20251202 <- "HY"
id <- intersect(colnames(mou),str_remove(colnames(int),"mou_"))
mou@meta.data[id,]$bansky_20251202 <- int@meta.data[paste0("mou_",id),]$bansky_202501202

mar$bansky_20251202 <- "HY"
id <- intersect(colnames(mar),str_remove(colnames(int),"mar_"))
mar@meta.data[id,]$bansky_20251202 <- int@meta.data[paste0("mar_",id),]$bansky_202501202

mac$bansky_20251202 <- "HY"
id <- intersect(colnames(mac),str_remove(colnames(int),"mac_"))
mac@meta.data[id,]$bansky_20251202 <- int@meta.data[paste0("mac_",id),]$bansky_202501202

Idents(mou) <- mou$bansky_20251202
Idents(mar) <- mar$bansky_20251202
Idents(mac) <- mac$bansky_20251202


if(T){
  #### macaque 
  meta3 <- meta2[meta2$specie %in% "macaque",]
  ann1 <- names(table(meta3$domian_modif2))
  ann2 <- names(table(meta3$bansky_202501202))
  ann1_result <- list()
  ann2_result<- list()
  for (i in 1:length(ann1)){
    ann2_result<- list()
    for (j in 1:length(ann2)) {
      overlap_result <- intersect(rownames(meta3[meta3$domian_modif2 %in% ann1[i],]),rownames(meta3[meta3$bansky_202501202 %in% ann2[j],]))
      union_result <- unique(c(rownames(meta3[meta3$domian_modif2 %in% ann1[i],]),rownames(meta3[meta3$bansky_202501202 %in% ann2[j],])))
      ann2_result[[ann2[j]]] <- length(overlap_result)/length(union_result)
    }
    ann1_result[[ann1[i]]] <- as.numeric(ann2_result)
  }
  
  result <- data.frame(ann1_result)
  rownames(result ) <- ann2
  mac_result <- result
  pheatmap::pheatmap(mac_result,cluster_rows = F,
                     cluster_cols = F,border_color = "white",
                     color = colorRampPalette(c("#EEEEEE","#FFEBEE","#EF9A9A","#F44336","#D32F2F"))(100))
  #### marmoset
  meta3 <- meta2[meta2$specie %in% "marmoset",]
  ann1 <- names(table(meta3$domian_modif2))
  ann2 <- names(table(meta3$bansky_202501202))
  ann1_result <- list()
  ann2_result<- list()
  for (i in 1:length(ann1)){
    ann2_result<- list()
    for (j in 1:length(ann2)) {
      overlap_result <- intersect(rownames(meta3[meta3$domian_modif2 %in% ann1[i],]),rownames(meta3[meta3$bansky_202501202 %in% ann2[j],]))
      union_result <- unique(c(rownames(meta3[meta3$domian_modif2%in% ann1[i],]),rownames(meta3[meta3$bansky_202501202 %in% ann2[j],])))
      ann2_result[[ann2[j]]] <- length(overlap_result)/length(union_result)
    }
    ann1_result[[ann1[i]]] <- as.numeric(ann2_result)
  }
  
  result <- data.frame(ann1_result)
  rownames(result ) <- ann2
  mar_result <- result
  pheatmap::pheatmap(mar_result,cluster_rows = F,
                     cluster_cols = F,border_color = "white",
                     color = colorRampPalette(c("#EEEEEE","#FFEBEE","#EF9A9A","#F44336","#D32F2F"))(100))
  
  #### macaque 
  meta3 <- meta2[meta2$specie %in% "mouse",]
  ann1 <- names(table(meta3$domian_modif2))
  ann2 <- names(table(meta3$bansky_202501202))
  ann1_result <- list()
  ann2_result<- list()
  for (i in 1:length(ann1)){
    ann2_result<- list()
    for (j in 1:length(ann2)) {
      overlap_result <- intersect(rownames(meta3[meta3$domian_modif2 %in% ann1[i],]),rownames(meta3[meta3$bansky_202501202 %in% ann2[j],]))
      union_result <- unique(c(rownames(meta3[meta3$domian_modif2 %in% ann1[i],]),rownames(meta3[meta3$bansky_202501202 %in% ann2[j],])))
      ann2_result[[ann2[j]]] <- length(overlap_result)/length(union_result)
    }
    ann1_result[[ann1[i]]] <- as.numeric(ann2_result)
  }
  
  result <- data.frame(ann1_result)
  rownames(result ) <- ann2
  mou_result <- result
  pheatmap::pheatmap(mou_result,cluster_rows = F,
                     cluster_cols = F,border_color = "white",
                     color = colorRampPalette(c("#EEEEEE","#FFEBEE","#EF9A9A","#F44336","#D32F2F"))(100))
}


###### arrange heatmap
path <- ""

nucle1 <- c("ARH","VMH","PVH","SO","TU","SCH","MBO","Pe","ME","BM","ZI","PM","PH","POA","DMH","AHN","LHA")
nucle2 <- c("ARH","VMH","PVH","SO","TU","SCH","MBO","Pe","ME","BM","ZI","PM","PH","POA","DMH","AHN","LHA","MarMG")

add_col <- data.frame(t(data.frame(rep(0,22))))
rownames(add_col) <- "MarMG"
colnames(add_col) <- colnames(mac_result)
mac_result2 <- rbind(mac_result,add_col)
mou_result2 <- rbind(mou_result,add_col)
domain <- paste0("D",c(7,10,18,19,22,17,15,12,20,13,16,2,4,14,6,3,5,11,8,1,9,21))

mac_result2 <- mac_result2[nucle2,domain]
mar_result2 <- mar_result[nucle2,domain]
mou_result2 <- mou_result2[nucle2,domain]

rownames(mac_result2) <- paste0("st:",rownames(mac_result2))
rownames(mar_result2) <- paste0("st:",rownames(mar_result2))
rownames(mou_result2) <- paste0("st:",rownames(mou_result2))

pheatmap::pheatmap(mac_result2,cluster_rows = F,breaks = seq(from = 0, to = 0.5, length.out = 100),
                   cluster_cols = F,border_color = "white",filename = paste0(path,"mac_result.pdf"),
                   color = colorRampPalette(c("#F5F5F5","#FFEBEE","#EF9A9A","#F44336","#D32F2F"))(100))
pheatmap::pheatmap(mar_result2,cluster_rows = F,breaks = seq(from = 0, to = 0.5, length.out = 100),
                   cluster_cols = F,border_color = "white",filename = paste0(path,"mar_result.pdf"),
                   color = colorRampPalette(c("#F5F5F5","#FFEBEE","#EF9A9A","#F44336","#D32F2F"))(100))
pheatmap::pheatmap(mou_result2,cluster_rows = F,breaks = seq(from = 0, to = 0.5, length.out = 100),
                   cluster_cols = F,border_color = "white",filename = paste0(path,"mou_result.pdf"),
                   color = colorRampPalette(c("#F5F5F5","#FFEBEE","#EF9A9A","#F44336","#D32F2F"))(100))

