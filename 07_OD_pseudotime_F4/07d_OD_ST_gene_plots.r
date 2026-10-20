library(ggplot2)

speci = 'macaque'
sel_chip <- "T790" #  marmoset:"T015_T965"  mouse:T315 macaque:T790
chip_id <- sel_chip
print("Beging to run mouse chip ...")
print(chip_id)
this_chip <- subset(mac_st, subset = chip == chip_id )



highlight_types <- c("COPC","OPC","DOPC", "OD3", "OD2", "OD1")

df <- this_chip@meta.data
df$color_group <- ifelse(
  df$celltype_pred %in% highlight_types,
  df$celltype_pred,
  "Other"
)

df$point_size <- ifelse(df$color_group == "Other", 0.3, 0.3)
df$point_alpha <- ifelse(df$color_group == "Other", 0.3, 1)

p <- ggplot(df, aes(x = rx, y = -ry)) +
  geom_point(
    data = subset(df, color_group == "Other"),
    aes(color = color_group),
    size = unique(subset(df, color_group == "Other")$point_size),
    alpha = unique(subset(df, color_group == "Other")$point_alpha)
  ) +
  geom_point(
    data = subset(df, color_group != "Other"),
    aes(color = color_group),
    size = unique(subset(df, color_group != "Other")$point_size),
    alpha = unique(subset(df, color_group != "Other")$point_alpha)
  ) +
  scale_color_manual(


    values =  c("OD1" = "#66C2A5","OD2" = "#FC8D62", "OD3" = "#8DA0CB","DOPC" = "#E78AC3",
         "OPC" = "#A6D854","COPC" = "#FFD92F","Other" = "#F0F0F0"),
    name = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank()
    # plot.margin = unit(c(-0.5, -0.5, -0.5, -0.5), "cm")
    # aspect.ratio = 1
  )

pdf(paste('ST_geneEXP', speci, 'chip', chip_id, "_cts_0928_2.pdf", sep='_'), h=2, w=2)
print(p)
dev.off()

######## ==============================================================================
#              Fig
######## ==============================================================================

nhp_dn <- c("Acvr2a", "Evi5", "Ttc28", "Gphn", "Stim1", "Rbfox2", "Nr3c2",
"Stxbp5l", "Sgip1", "Ncoa1", "Pcnx", "Ppp3ca", "Sptan1", "Slc8a1",
"Cacna2d1", "Fnbp4", "Erc1", "Arhgef12", "Fign", "Agap1", "Dgki",
"Ppp1r12a", "Hivep1", "Setbp1", "Sh3d19", "Plekha2")
nhp_dn_upper <- dput(toupper(nhp_dn))
nhp_up <- c("Pak2", "Cep170", "Ankrd12", "Serpini1", "Mtmr10", "Itprid2",
"Synm", "Gm12248", "Arhgef37", "Cep85")

nhp_up_upper <- c("PAK2", "CEP170", "ANKRD12", "SERPINI1", "MTMR10", "ITPRID2",
"STNM", "GM12248", "ARHEF37", "CEP85")
speci = 'mouse'
sel_chip <- "T315" # #  marmoset:"T015_T965"  mouse:T315 macaque:T790
chip_id <- sel_chip
# onechip <- "T295"
print("Beging to run mouse chip ...")
print(chip_id)
this_chip <- subset(mou_st, subset = chip == chip_id )

x<- -(this_chip@meta.data$rx)
y<- -(this_chip@meta.data$ry) # y<- -(this_chip@meta.data$ry)
# y<- this_chip@meta.data$ry # y<- -(this_chip@meta.data$ry)


library(ggplot2)
library(patchwork)
for(gene_name in nhp_up) {
    # gene_name <- 'Cck'
    if(gene_name %in% row.names(this_chip@assays$RNA@data)){
        print(gene_name)
        gene_expression<-this_chip@assays$RNA@data[gene_name,]
        data<-data.frame(x=x,y=y,gene_expression=gene_expression)
        data <- data[order(data$gene_expression), ]

        p_scatter <- ggplot(data, aes(x = x, y = y, color = gene_expression)) +
            geom_point(size = 1.5,  alpha = 1) +
            scale_color_viridis_c(option = "plasma", begin = 0.2, end = 1) +
            labs(
                title = gene_name
            ) +
          theme_minimal() +
          theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),

            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),  
            legend.position = "top",  
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 12),
            panel.grid.major = element_blank(),  
            panel.grid.minor = element_blank(),  
            panel.border = element_rect(colour = "black", fill = NA, size = 0.5)  
          )
        
        data_normalized <- data %>%
        mutate(x_normalized = (x - min(x)) / (max(x) - min(x)))  
        p_trend <-  ggplot(data_normalized, aes(x = x_normalized, y = gene_expression)) +
            stat_summary_bin(
                fun = "mean",
                geom = "line",
                color = "red",
                linewidth = 1,
                bins = 30  
            ) +
            scale_x_continuous(
                limits = c(0, 1),  
                breaks = seq(0, 1, by = 0.2), 
                #labels = scales::percent_format()  
            ) +
            labs(y = "Mean Expression", x = "Normalized Spatial X-coordinate") +
            theme_minimal(base_size = 12) +
            theme(
                axis.title = element_text(size = 10, face = "bold"),
                axis.line = element_line(color = "black", linewidth = 0.5),
                panel.grid = element_blank()
            )

        p_combined <- p_scatter / p_trend +
        plot_layout(heights = c(3, 1))  
        print(p_combined)
  
    }
}


for(gene_name in nhp_dn) {
    print(gene_name)
    # gene_name <- 'Cck'
    
    if(gene_name %in% row.names(this_chip@assays$RNA@data)){
        gene_expression<-this_chip@assays$RNA@data[gene_name,]
        data<-data.frame(x=x,y=y,gene_expression=gene_expression)
        data <- data[order(data$gene_expression), ]
        p_scatter <- ggplot(data, aes(x = x, y = y, color = gene_expression)) +
            geom_point(size = 1.5,  alpha = 1) +
            scale_color_viridis_c(option = "plasma", begin = 0.2, end = 1) +
            labs(
                title = gene_name
            ) +
          theme_minimal() +
          theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),

            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  
            legend.position = "top",  
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 14),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),  
            panel.border = element_rect(colour = "black", fill = NA, size = 0.5)  
          )
     
        p_trend <- ggplot(data, aes(x = x, y = gene_expression)) +
            stat_summary_bin(
                fun = "mean",
                geom = "line",
                color = "red",
                linewidth = 1,
                bins = 30  
            ) +
            labs(y = "Mean Expression", x = "Spatial X-coordinate (bins)") +
            theme_minimal(base_size = 16) +  
            theme(
                axis.title = element_text(size = 16, face = "bold", color = "black"),
                axis.text.x =  element_blank(),
                axis.text.y = element_text(size = 16, color = "black"),
                axis.line = element_line(color = "black", linewidth = 0.5),  
                legend.title = element_text(face = "bold", size = 22),
                legend.text = element_text(size = 22),
                legend.position = "none", 
                panel.grid.major = element_blank(),  
                panel.grid.minor = element_blank(),
                plot.margin = unit(c(5, 5, 5, 5), "mm")  
            )
        p_combined <- p_scatter / p_trend +
            plot_layout(heights = c(3, 1))  

        pdf(paste('./dn_overla_mou/ST_geneEXP' , gene_name , speci , 'chip', chip_id ,"_dn_overlay_0729.pdf",sep='_'))
        # print(p)
        print(p_combined)
        dev.off()
    }
}




######## ==============================================================================
#              Fig legend
######## ==============================================================================



cell_od <- c("OD1","OD2","OD3","DOPC","OPC","COPC")
color_od <- c("#FFE0B2","#FFE082","#FFC107","#FFB74D","#FB8C00","#E65100")



library(ggplot2)
library(cowplot) 
library(grid)

use_cor <- c("#FFE0B2","#FFE082","#FFC107","#FFB74D","#FB8C00","#E65100")
names(use_cor) <- c("OD1","OD2","OD3","DOPC","OPC","COPC")

dummy_plot <- ggplot(data.frame(group = names(use_cor),
                               value = 1:6,
                               color = use_cor)) +
  geom_point(aes(x = 0, y = 0, color = group),
             size = 0, alpha = 0) +
  scale_color_manual(
    values = use_cor,
    name = NULL,
    guide = guide_legend(
      override.aes = list(size = 5, shape = 19, alpha = 1),
      keywidth = unit(12, "mm")
    )
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.spacing.y = unit(0.5, "cm"),
    legend.text = element_text(size = 8, margin = margin(l = 8))
  )

legend_only <- get_legend(dummy_plot)
save_plot("legend.pdf", legend_only, base_width = 1, base_height = 2)

pdf("legend_direct.pdf", width = 1, height = 2)
grid.newpage()
grid.draw(legend_only)
dev.off()

legend_grob <- ggplotGrob(ggplot() + annotation_custom(legend_only) + theme_void())
ggsave("legend_wrapped.pdf", legend_grob, width = 1, height = 2)



######## ==============================================================================
#              Fig Heatmap
######## ==============================================================================

library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize)




nhp_dn <- c("Acvr2a",  "Ncoa1", "Pcnx", "Sptan1", "Slc8a1",
"Cacna2d1", "Fnbp4", "Erc1",  "Hivep1", "Setbp1")

nhp_up <- c("Pak2", "Cep170", "Ankrd12", "Serpini1", "Mtmr10", "Itprid2",
"Synm", "Gm12248", "Arhgef37", "Cep85")
specie_sel <- 'mouse'
setwd(paste('./monocle3_out/' ,specie_sel,'_monocle3',sep=''))
# genes <- read.csv(paste(specie_sel,"Trajectory_genes.sc.csv",sep=''),header=T)
#####################################  mouse up  ########
pbmc <- readRDS(paste(specie_sel,"sco_pseudotime.sc.rds",sep='_'))

pbmc$pseudotime_bin <- cut(
  pbmc$pseudotime,
  breaks = 100,
  labels = paste0("Bin_", 1:100)
)
table(pbmc$pseudotime_bin)

gene_set <- nhp_up

bin_avg_exp <- AverageExpression(
  pbmc,
  features = gene_set,
  group.by = "pseudotime_bin",
  assays = "RNA"
)$RNA


bin_avg_exp_scaled <- t(scale(t(bin_avg_exp)))
pdf(paste("./OD_Pseudotime_Heatmap/" ,specie_sel, "_bin_avg_exp_scaled_1_up.pdf", sep=''),h=1.5,w=7)
# 创建热图
Heatmap(
  bin_avg_exp_scaled,
  name = "Average Expression\n (Z-score)",
  col = colorRamp2(c(-2, 0, 2), c("steelblue", "white", "firebrick")),
  cluster_rows = FALSE,
  cluster_columns = FALSE, 
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_title = paste(specie_sel," Pseudotime Bins (up)", sep=''),
  row_title = "",
  show_row_names = TRUE,
  show_column_names = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 6, fontface = "bold"),  
    title_position = "leftcenter-rot",
    labels_gp = gpar(fontsize = 8)  
  )
)
dev.off()

pdf(paste("./OD_Pseudotime_Heatmap/" ,specie_sel, "_bin_avg_exp_scaled_1_up_noGeneName.pdf", sep=''),h=1,w=7)
Heatmap(
  bin_avg_exp_scaled,
  name = "Average Expression\n (Z-score)",
  col = colorRamp2(c(-2, 0, 2), c("steelblue", "white", "firebrick")),
  cluster_rows = FALSE,
  cluster_columns = FALSE, 
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_title = paste(specie_sel," Pseudotime Bins (up)", sep=''),
  row_title = "",
  show_row_names = FALSE,
  show_column_names = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 6, fontface = "bold"),  
    title_position = "leftcenter-rot",
    labels_gp = gpar(fontsize = 8)  
  )
)
dev.off()
###########  mouse dn  ########

gene_set <- nhp_dn

bin_avg_exp <- AverageExpression(
  pbmc,
  features = gene_set,
  group.by = "pseudotime_bin",
  assays = "RNA" 
)$RNA

bin_avg_exp_scaled <- t(scale(t(bin_avg_exp)))
pdf(paste("./OD_Pseudotime_Heatmap/" ,specie_sel, "_bin_avg_exp_scaled_1_dn.pdf", sep=''),h=1.5,w=7)
Heatmap(
  bin_avg_exp_scaled,
  name = "Average Expression\n (Z-score)",
  col = colorRamp2(c(-2, 0, 2), c("steelblue", "white", "firebrick")),
  cluster_rows = FALSE,
  cluster_columns = FALSE, 
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_title = paste(specie_sel," Pseudotime Bins  (down)", sep=''),
  row_title = " ",
  show_row_names = TRUE,
  show_column_names = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 6, fontface = "bold"),  
    title_position = "leftcenter-rot",
    labels_gp = gpar(fontsize = 8)  
  )
)
dev.off()

pdf(paste("./OD_Pseudotime_Heatmap/" ,specie_sel, "_bin_avg_exp_scaled_1_dn_noGeneName.pdf", sep=''),h=1,w=7)
Heatmap(
  bin_avg_exp_scaled,
  name = "Average Expression\n (Z-score)",
  col = colorRamp2(c(-2, 0, 2), c("steelblue", "white", "firebrick")),
  cluster_rows = FALSE,
  cluster_columns = FALSE, 
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_title = paste(specie_sel," Pseudotime Bins  (down)", sep=''),
  row_title = " ",
  show_row_names = FALSE,
  show_column_names = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 6, fontface = "bold"),  
    title_position = "leftcenter-rot",
    labels_gp = gpar(fontsize = 8)  
  )
)
dev.off()

######## ==============================================================================
#              Fig OD umap
######## ==============================================================================
specie_sel <- 'all'
library(spatstat)
library(ggplot2)
library(dplyr)
library(Seurat)
library(maSigPro)

# genes <- read.csv(paste(specie_sel,"Trajectory_genes.sc.csv",sep=''),header=T)
library(qs)
oligo <- qread("./data/od_opc_update_umap.qs")

pbmc <- oligo


pbmc@meta.data$annotation[is.na(pbmc@meta.data$annotation)] <- "Unknown"

pbmc@meta.data$annotation <- factor(
  pbmc@meta.data$annotation,
  levels = names(which(table(pbmc@meta.data$annotation) > 0))
)

table(pbmc@meta.data$annotation)

pbmc@meta.data$annotation <- factor(pbmc@meta.data$annotation, levels=c("COPC","OPC" ,"DOPC", "OD3" ,"OD2" ,"OD1" ))
levels(pbmc@meta.data$annotation)

cor <- c("OD1" = "#FFE0B2",
         "OD2" = "#FFE082",
         "OD3" = "#FFC107",
         "DOPC" = "#FFB74D",
         "OPC" = "#FB8C00",
         "COPC" = "#FB8C00"
         )


cor <- c("OD1" = "#66C2A5",
         "OD2" = "#FC8D62",
         "OD3" = "#8DA0CB",
         "DOPC" = "#E78AC3",
         "OPC" = "#A6D854",
         "COPC" = "#FFD92F"
         )

p <- DimPlot(pbmc, group.by='annotation',pt.size = 0.1,

        cols = cor)+

  theme(axis.text = element_text(size=8,colour = 'black'),

        axis.title = element_text(size=8,colour = 'black'),

        panel.background=element_rect(fill='transparent', color='black',size = 1),

        legend.key=element_blank(),legend.key.size  = unit(0.5,'cm'),

        legend.key.height = unit(0.5,'cm'),legend.key.width  = unit(0.5,'cm'),

        legend.text = element_text(size=8,colour = 'black'),

        plot.margin = unit(c(0,0,-0.1,-0.1),'cm'))

ggsave("./OD_reUMAP_0923_3spe_spli_3_1.pdf",p,width =7,height = 6,units = 'cm')


p <- DimPlot(pbmc, group.by='annotation',split.by = 'specie',pt.size = 0.1,

        cols = cor)+

  theme(axis.text = element_text(size=8,colour = 'black'),

        axis.title = element_text(size=8,colour = 'black'),
        # axis.title = element_blank(),
        panel.background=element_rect(fill='transparent', color='black',size = 1),

        legend.key=element_blank(),legend.key.size  = unit(0.5,'cm'),

        legend.key.height = unit(0.1,'cm'),legend.key.width  = unit(0.5,'cm'),

        legend.text = element_text(size=7,colour = 'black'),

        plot.margin = unit(c(0,0,-0.1,-0.1),'cm'),
            plot.title = element_blank() 
  )

ggsave("./OD_reUMAP_0923_3spe_spli_2.1_1.pdf",p,width = 14,height = 5,units = 'cm')

