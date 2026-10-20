if(T){
  library(Seurat)
  library(qs)
  library(viridis)
  library(ComplexHeatmap)
  library(ggplot2)
  library(RColorBrewer)
  library(COSG)
}
mar <- qread("mar.qs")
mou <- qread("mou.qs")
mac <- qread("mac.qs")

mou <- NormalizeData(mou)
mar <- NormalizeData(mar)
mac <- NormalizeData(mac)

out_path <- ""

n <- 0.5
p1 <- DimPlot(mou,cells =  rownames(mou_meta[mou_meta$chip %in% "T316",]),cols = "#EEEEEE",pt.size = 1,
              cells.highlight = WhichCells(mou,expression =  Sst > n & Slc17a6 > n),reduction = "umap2",raster = F,sizes.highlight = 1)+theme_nothing()
p2 <- DimPlot(mar,cells =  rownames(mar_meta[mar_meta$chip2 %in% "T016",]),cols = "#EEEEEE",pt.size = 1,
              cells.highlight = WhichCells(mar,expression =  SST > n & SLC17A6 > n ),reduction = "umap2",raster = F,sizes.highlight = 1)+theme_nothing()
p3 <- DimPlot(mac,cells =  rownames(mac_meta[mac_meta$chip %in% "T792",]),cols = "#EEEEEE",
              cells.highlight = WhichCells(mac,expression =  SST > n & SLC17A6 > n ),reduction = "umap2",raster = F,sizes.highlight = 1)+theme_nothing()
p_gene <- p3|p2|p1
p_gene

ggsave(p_gene,filename = paste0(out_path,"En7_sst_Slc17a6.pdf"),width = 12,height = 4)
ggsave(p1,filename = paste0(out_path,"En7_sst_slc17a6_mou.pdf"),width = 4,height = 4)
ggsave(p2,filename = paste0(out_path,"En7_sst_slc17a6_mar.pdf"),width = 4,height = 4)
ggsave(p3,filename = paste0(out_path,"En7_sst_slc17a6_mac.pdf"),width = 4,height = 4)

####LB
Idents(mou) <- mou$laber_transfer_1109
Idents(mac) <- mac$laber_transfer_1109
Idents(mar) <- mar$laber_transfer_1109
p1 <- DimPlot(mou,cells =  rownames(mou_meta[mou_meta$chip %in% "T316",]),sizes.highlight = 1,cols = "#EEEEEE",cols.highlight = col[["EN7"]],pt.size = 1,
              cells.highlight = intersect(rownames(mou_meta[mou_meta$chip %in% "T316",]),WhichCells(mou,idents = "EN7_Sst/Esr1")))+theme_nothing()
p2 <- DimPlot(mar,cells =  rownames(mar_meta[mar_meta$chip2 %in% "T016",]),sizes.highlight =1,cols = "#EEEEEE",cols.highlight = col[["EN7"]],pt.size = 1,
              cells.highlight = intersect(rownames(mar_meta[mar_meta$chip2 %in% "T016",]),WhichCells(mar,idents = "EN7_Sst/Esr1")))+theme_nothing()
p3 <- DimPlot(mac,cells =  rownames(mac_meta[mac_meta$chip %in% "T792",]),sizes.highlight = 1,cols = "#EEEEEE",cols.highlight = col[["EN7"]],
              cells.highlight = intersect(rownames(mac_meta[mac_meta$chip %in% "T792",]),WhichCells(mac,idents = "EN7_Sst/Esr1")))+theme_nothing()

p_lb <- p3|p2|p1
p_lb

ggsave(p_lb,filename = paste0(out_path,"En7_lb.pdf"),width = 12,height = 4)
ggsave(p1,filename = paste0(out_path,"En7_lb_mou.pdf"),width = 4,height = 4)
ggsave(p2,filename = paste0(out_path,"En7_lb_mar.pdf"),width = 4,height = 4)
ggsave(p3,filename = paste0(out_path,"En7_lb_mac.pdf"),width = 4,height = 4)





