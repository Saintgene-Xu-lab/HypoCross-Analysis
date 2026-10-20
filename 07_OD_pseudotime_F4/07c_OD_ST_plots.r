library(pheatmap)
library(tidyr)
library(ggpubr)
library(ggalluvial)
library(doParallel)
library(foreach)
library(ggExtra)
library(ggalluvial)
library(ggplot2)
library(ggrepel)
library(grDevices)
library(iterators)
library(methods)
library(parallel)
library(purrr)  
library(dplyr)  
library(readxl)
library(patchwork)

#####
speci <- 'macaque'

folder_path <- "./macaque_non20x_40xOD_5xTC_0.5xEN1_0.2xIN1_IN13mk_85c_csv/"
csv_files <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)

mou_st <- readRDS('./spatial_rds_20250328/macaque.rds')

mou_st@meta.data <- mou_st@meta.data %>%
  mutate(cell = row.names(mou_st@meta.data))
mou_st@meta.data <- mou_st@meta.data %>%
  mutate(cell = row.names(mou_st@meta.data))
csv_files <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)

merged_csv_data <- csv_files %>%
  map(~ read.csv(.x) %>% dplyr::select("cell", "celltype_pred")) %>%
  reduce(bind_rows) %>%
  distinct(cell, .keep_all = TRUE)
mou_st@meta.data <- mou_st@meta.data %>%
  left_join(merged_csv_data, by = "cell")

head(mou_st@meta.data)

rownames(mou_st@meta.data) <- mou_st@meta.data$cell
##############

highlight_types <- c("COPC","OPC","DOPC", "OD3", "OD2", "OD1")
df <- this_chip@meta.data
df$color_group <- ifelse(
  df$celltype_pred %in% highlight_types, 
  df$celltype_pred, 
  "Other"
)

pdf(paste('ST_geneEXP' ,  speci , 'chip', chip_id ,"_cts_0729_3.pdf",sep='_'),h=4,w=4)
ggplot(df, aes(x = rx, y = -ry, color = color_group)) +
  geom_point(size = 0.4) +
  scale_color_manual(
    values = c( "COPC" = "#eeb401","OPC" = "#942d8d", "DOPC" = "#2e409a", "OD3" = "#d75427", "OD2" = "#00abf0", "OD1" = "#3cb346",  "Other" = "gray80"),
    name = "Cell Type"
  ) +
  theme_minimal(base_size = 22) +  
  theme(
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.text.x = element_text(
      angle = 90, 
      hjust = 1,  
      vjust = 1, 
      size = 16,
      color = "black"
    ),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.5),  
    legend.title = element_text(face = "bold", size = 22),
    legend.text = element_text(size = 22),
    legend.position = "none", 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    plot.margin = unit(c(5, 5, 5, 5), "mm")  
  ) 

dev.off()

nhp_dn <- c("Acvr2a", "Evi5", "Ttc28", "Gphn", "Stim1", "Rbfox2", "Nr3c2",
"Stxbp5l", "Sgip1", "Ncoa1", "Pcnx", "Ppp3ca", "Sptan1", "Slc8a1",
"Cacna2d1", "Fnbp4", "Erc1", "Arhgef12", "Fign", "Agap1", "Dgki",
"Ppp1r12a", "Hivep1", "Setbp1", "Sh3d19", "Plekha2")
nhp_dn_upper <- dput(toupper(nhp_dn))

nhp_up <- c("PAK2", "CEP170", "ANKRD12", "SERPINI1", "MTMR10", "ITPRID2",
"STNM", "GM12248", "ARHEF37", "CEP85")

sel_chip <- "T790"
chip_id <- sel_chip
# onechip <- "T295"
print("Beging to run mouse chip ...")
print(chip_id)
this_chip <- subset(mou_st, subset = chip == chip_id )
this_chip <- subset(this_chip, subset = celltype_pred %in% highlight_types )

x<-this_chip@meta.data$rx
y<- -(this_chip@meta.data$ry) 


for(gene_name in nhp_up) {
	# gene_name <- 'Cck'
	if(gene_name %in% row.names(this_chip@assays$RNA@data)){
		gene_expression<-this_chip@assays$RNA@data[gene_name,]
		data<-data.frame(x=x,y=y,gene_expression=gene_expression)
		data <- data[order(data$gene_expression), ]

		p_scatter <- ggplot(data, aes(x = x, y = y, color = gene_expression)) + 
		  geom_point(size = 1.5) + 
		  scale_color_viridis_c(option = "plasma", begin = 0.2, end = 1) +
		  labs(
		    title = gene_name
		  ) +
		  theme_minimal() +
		  theme(
		    axis.title.x = element_blank(),  
		    axis.text.x = element_blank(),   
		    panel.grid = element_blank(),
		    legend.position = "right"
		  )

		p_trend <- ggplot(data, aes(x = x, y = gene_expression)) +
		  stat_summary_bin(
		    fun = "mean",
		    geom = "line", 
		    color = "red",
		    linewidth = 1,
		    bins = 30  
		  ) +
		  labs(y = "Mean Expression", x = "Spatial X-coordinate") +
		  theme_minimal(base_size = 12) + 
		  theme(
		    axis.title = element_text(size = 12, face = "bold", color = "black"),
		    axis.text.x = element_text(
		      angle = 0, 
		      hjust = 1,  
		      vjust = 1, 
		      size = 20,
		      color = "black"
		    ),
		    axis.text.y = element_text(size = 20, color = "black"),
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

		pdf(paste('./up_overla_mac/ST_geneEXP' , gene_name , speci , 'chip', chip_id ,"up_overlay_0729_revY.pdf",sep='_'))
		print(p_combined)
		dev.off()
	}
}


for(gene_name in nhp_dn_upper) {
	print(gene_name)
	# gene_name <- 'Cck'
	if(gene_name %in% row.names(this_chip@assays$RNA@data)){
		gene_expression<-this_chip@assays$RNA@data[gene_name,]
		data<-data.frame(x=x,y=y,gene_expression=gene_expression)
		data <- data[order(data$gene_expression), ]
		p_scatter <- ggplot(data, aes(x = x, y = y, color = gene_expression)) + 
		  geom_point(size = 1.5) + 
		  scale_color_viridis_c(option = "plasma", begin = 0.2, end = 1) +
		  labs(
		    title = gene_name
		  ) +
		  theme_minimal() +
		  theme(
		    axis.title.x = element_blank(),  
		    axis.text.x = element_blank(),   
		    panel.grid = element_blank(),
		    legend.position = "right"
		  )

		p_trend <- ggplot(data, aes(x = x, y = gene_expression)) +
		  stat_summary_bin(
		    fun = "mean",
		    geom = "line", 
		    color = "red",
		    linewidth = 1,
		    bins = 30  
		  ) +
		  labs(y = "Mean Expression", x = "Spatial X-coordinate") +
		  theme_minimal(base_size = 12) +  
		  theme(
=		    axis.title = element_text(size = 12, face = "bold", color = "black"),
		    axis.text.x = element_text(
		      angle = 0, 
		      hjust = 1,  
		      vjust = 1, 
		      size = 20,
		      color = "black"
		    ),
		    axis.text.y = element_text(size = 20, color = "black"),
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
		# print(p_combined)

		pdf(paste('./dn_overla_mac/ST_geneEXP' , gene_name , speci , 'chip', chip_id ,"_dn_overlay_0729.pdf",sep='_'))
		# print(p)
		print(p_combined)
		dev.off()
	}
}


################################################################
#         marmoset
################################################################



#####
speci <- 'marmoset'

folder_path <- "./marmoset_non8x_10xOD_0.4xEN1_0.2xIN1_IN13mk_85c_csv/"

csv_files <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)

mou_st <- readRDS('./spatial_rds_20250328/marmoset.rds')

mou_st@meta.data <- mou_st@meta.data %>%
  mutate(cell = row.names(mou_st@meta.data))
mou_st@meta.data <- mou_st@meta.data %>%
  mutate(cell = row.names(mou_st@meta.data))
csv_files <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)



merged_csv_data <- csv_files %>%
  map(~ read.csv(.x) %>% dplyr::select("cell", "celltype_pred")) %>%
  reduce(bind_rows) %>%
  distinct(cell, .keep_all = TRUE)
mou_st@meta.data <- mou_st@meta.data %>%
  left_join(merged_csv_data, by = "cell")

head(mou_st@meta.data)

rownames(mou_st@meta.data) <- mou_st@meta.data$cell
##############


nhp_dn <- c("Acvr2a", "Evi5", "Ttc28", "Gphn", "Stim1", "Rbfox2", "Nr3c2",
"Stxbp5l", "Sgip1", "Ncoa1", "Pcnx", "Ppp3ca", "Sptan1", "Slc8a1",
"Cacna2d1", "Fnbp4", "Erc1", "Arhgef12", "Fign", "Agap1", "Dgki",
"Ppp1r12a", "Hivep1", "Setbp1", "Sh3d19", "Plekha2")
nhp_dn_upper <- dput(toupper(nhp_dn))

nhp_up <- c("PAK2", "CEP170", "ANKRD12", "SERPINI1", "MTMR10", "ITPRID2",
"STNM", "GM12248", "ARHEF37", "CEP85")

sel_chip <- "T015_T965"
chip_id <- sel_chip
# onechip <- "T295"
print("Beging to run mouse chip ...")
print(chip_id)
this_chip <- subset(mou_st, subset = chip == chip_id )
this_chip <- subset(this_chip, subset = celltype_pred %in% highlight_types )

x<- this_chip@meta.data$rx
y<- this_chip@meta.data$ry 

for(gene_name in nhp_up) {
	# gene_name <- 'Cck'
	if(gene_name %in% row.names(this_chip@assays$RNA@data)){
		gene_expression<-this_chip@assays$RNA@data[gene_name,]
		data<-data.frame(x=x,y=y,gene_expression=gene_expression)
		data <- data[order(data$gene_expression), ]

		
		p_scatter <- ggplot(data, aes(x = x, y = y, color = gene_expression)) + 
		  geom_point(size = 1.5) + 
		  scale_color_viridis_c(option = "plasma", begin = 0.2, end = 1) +
		  labs(
		    title = gene_name
		  ) +
		  theme_minimal() +
		  theme(
		    axis.title.x = element_blank(), 
		    axis.text.x = element_blank(),  
		    panel.grid = element_blank(),
		    legend.position = "right"
		  )

		p_trend <- ggplot(data, aes(x = x, y = gene_expression)) +
		  stat_summary_bin(
		    fun = "mean",
		    geom = "line", 
		    color = "red",
		    linewidth = 1,
		    bins = 30  
		  ) +
		  labs(y = "Mean Expression", x = "Spatial X-coordinate") +
		  theme_minimal(base_size = 12) + 
		  theme(
		    axis.title = element_text(size = 12, face = "bold", color = "black"),
		    axis.text.x = element_text(
		      angle = 0, 
		      hjust = 1,  
		      vjust = 1, 
		      size = 20,
		      color = "black"
		    ),
		    axis.text.y = element_text(size = 20, color = "black"),
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

		pdf(paste('./up_overla_mar/ST_geneEXP' , gene_name , speci , 'chip', chip_id ,"up_overlay_0729.pdf",sep='_'))
		# print(p)
		print(p_combined)
		dev.off()
	}
}


for(gene_name in nhp_dn_upper) {
	print(gene_name)
	# gene_name <- 'Cck'
	if(gene_name %in% row.names(this_chip@assays$RNA@data)){
		gene_expression<-this_chip@assays$RNA@data[gene_name,]
		data<-data.frame(x=x,y=y,gene_expression=gene_expression)
		data <- data[order(data$gene_expression), ]
		p_scatter <- ggplot(data, aes(x = x, y = y, color = gene_expression)) + 
		  geom_point(size = 1.5) + 
		  scale_color_viridis_c(option = "plasma", begin = 0.2, end = 1) +
		  labs(
		    title = gene_name
		  ) +
		  theme_minimal() +
		  theme(
		    axis.title.x = element_blank(), 
		    axis.text.x = element_blank(),   
		    panel.grid = element_blank(),
		    legend.position = "right"
		  )

		p_trend <- ggplot(data, aes(x = x, y = gene_expression)) +
		  stat_summary_bin(
		    fun = "mean",
		    geom = "line", 
		    color = "red",
		    linewidth = 1,
		    bins = 30  
		  ) +
		  labs(y = "Mean Expression", x = "Spatial X-coordinate") +
		  theme_minimal(base_size = 12) +  
		  theme(
		    axis.title = element_text(size = 12, face = "bold", color = "black"),
		    axis.text.x = element_text(
		      angle = 0, 
		      hjust = 1,  
		      vjust = 1, 
		      size = 20,
		      color = "black"
		    ),
		    axis.text.y = element_text(size = 20, color = "black"),
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
		# print(p_combined)

		pdf(paste('./dn_overla_mar/ST_geneEXP' , gene_name , speci , 'chip', chip_id ,"_dn_overlay_0729.pdf",sep='_'))
		print(p_combined)
		dev.off()
	}
}
################################################################
library(ggplot2)
highlight_types <- c("COPC","OPC","DOPC", "OD3", "OD2", "OD1")
df <- this_chip@meta.data
df$color_group <- ifelse(
  df$celltype_pred %in% highlight_types, 
  df$celltype_pred, 
  "Other"
)

pdf(paste('ST_geneEXP' ,  speci , 'chip', chip_id ,"_cts_0729_3.pdf",sep='_'),h=4,w=4)
ggplot(df, aes(x = rx, y = ry, color = color_group)) +
  geom_point(size = 0.4) +
  scale_color_manual(
    values = c( "COPC" = "#eeb401","OPC" = "#942d8d", "DOPC" = "#2e409a", "OD3" = "#d75427", "OD2" = "#00abf0", "OD1" = "#3cb346",  "Other" = "gray80"),
    name = "Cell Type"
  ) +
  theme_minimal(base_size = 22) +  
  theme(
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.text.x = element_text(
      angle = 90, 
      hjust = 1,  
      vjust = 1, 
      size = 16,
      color = "black"
    ),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.5),  
    legend.title = element_text(face = "bold", size = 22),
    legend.text = element_text(size = 22),
    legend.position = "none", 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    plot.margin = unit(c(5, 5, 5, 5), "mm")  
  ) 

dev.off()




################################################################
#         mouse
################################################################


#####
speci <- 'mouse'
folder_path <- "./mouse_non8x_56xOD_2xTC_0.4xEN1_0.2xIN1_IN13mk_85c_csv/"

csv_files <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)

mou_st <- readRDS('./spatial_rds_20250328/mouse.rds')

mou_st@meta.data <- mou_st@meta.data %>%
  mutate(cell = row.names(mou_st@meta.data))
mou_st@meta.data <- mou_st@meta.data %>%
  mutate(cell = row.names(mou_st@meta.data))
csv_files <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)

merged_csv_data <- csv_files %>%
  map(~ read.csv(.x) %>% dplyr::select("cell", "celltype_pred")) %>%
  reduce(bind_rows) %>%
  distinct(cell, .keep_all = TRUE)
mou_st@meta.data <- mou_st@meta.data %>%
  left_join(merged_csv_data, by = "cell")

head(mou_st@meta.data)

rownames(mou_st@meta.data) <- mou_st@meta.data$cell
##############

highlight_types <- c("COPC","OPC","DOPC", "OD3", "OD2", "OD1")
df <- this_chip@meta.data
df$color_group <- ifelse(
  df$celltype_pred %in% highlight_types, 
  df$celltype_pred, 
  "Other"
)

pdf(paste('ST_geneEXP' ,  speci , 'chip', chip_id ,"_cts_0729_3.pdf",sep='_'),h=4,w=5)
ggplot(df, aes(x = -rx, y = -ry, color = color_group)) +
  geom_point(size = 0.4) +
  scale_color_manual(
    values = c( "COPC" = "#eeb401","OPC" = "#942d8d", "DOPC" = "#2e409a", "OD3" = "#d75427", "OD2" = "#00abf0", "OD1" = "#3cb346",  "Other" = "gray80"),
    name = "Cell Type"
  ) +
  theme_minimal(base_size = 22) +  
  theme(
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.text.x = element_text(
      angle = 90, 
      hjust = 1,  
      vjust = 1, 
      size = 16,
      color = "black"
    ),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.5),  
    legend.title = element_text(face = "bold", size = 22),
    legend.text = element_text(size = 22),
    legend.position = "none",  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    plot.margin = unit(c(5, 5, 5, 5), "mm")  
  ) 

dev.off()



nhp_dn <- c("Acvr2a", "Evi5", "Ttc28", "Gphn", "Stim1", "Rbfox2", "Nr3c2",
"Stxbp5l", "Sgip1", "Ncoa1", "Pcnx", "Ppp3ca", "Sptan1", "Slc8a1",
"Cacna2d1", "Fnbp4", "Erc1", "Arhgef12", "Fign", "Agap1", "Dgki",
"Ppp1r12a", "Hivep1", "Setbp1", "Sh3d19", "Plekha2")
nhp_up <- c("Pak2", "Cep170", "Ankrd12", "Serpini1", "Mtmr10", "Itprid2",
"Synm", "Gm12248", "Arhgef37", "Cep85")


sel_chip <- "T316"
chip_id <- sel_chip
print("Beging to run mouse chip ...")
print(chip_id)
this_chip <- subset(mou_st, subset = chip == chip_id )
this_chip <- subset(this_chip, subset = celltype_pred %in% highlight_types )

x<- -(this_chip@meta.data$rx)
y<- -(this_chip@meta.data$ry) 
library(ggplot2)
library(patchwork)
for(gene_name in nhp_up) {
	# gene_name <- 'Cck'
	if(gene_name %in% row.names(this_chip@assays$RNA@data)){
		gene_expression<-this_chip@assays$RNA@data[gene_name,]
		data<-data.frame(x=x,y=y,gene_expression=gene_expression)
		data <- data[order(data$gene_expression), ]

		p_scatter <- ggplot(data, aes(x = x, y = y, color = gene_expression)) + 
		  geom_point(size = 1.5) + 
		  scale_color_viridis_c(option = "plasma", begin = 0.2, end = 1) +
		  labs(
		    title = gene_name
		  ) +
		  theme_minimal() +
		  theme(
		    axis.title.x = element_blank(),  
		    axis.text.x = element_blank(),  
		    panel.grid = element_blank(),
		    legend.position = "right"
		  )

		p_trend <- ggplot(data, aes(x = x, y = gene_expression)) +
		  stat_summary_bin(
		    fun = "mean",
		    geom = "line", 
		    color = "red",
		    linewidth = 1,
		    bins = 30  
		  ) +
		  labs(y = "Mean Expression", x = "Spatial X-coordinate") +
		  theme_minimal(base_size = 12) +  
		  theme(
		    axis.title = element_text(size = 12, face = "bold", color = "black"),
		    axis.text.x = element_text(
		      angle = 0, 
		      hjust = 1,  
		      vjust = 1, 
		      size = 20,
		      color = "black"
		    ),
		    axis.text.y = element_text(size = 20, color = "black"),
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

		pdf(paste('./up_overla_mou/ST_geneEXP' , gene_name , speci , 'chip', chip_id ,"up_overlay_0729.pdf",sep='_'),w=8,h=6)
		print(p_combined)
		dev.off()
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
		  geom_point(size = 1.5) + 
		  scale_color_viridis_c(option = "plasma", begin = 0.2, end = 1) +
		  labs(
		    title = gene_name
		  ) +
		  theme_minimal() +
		  theme(
		    axis.title.x = element_blank(),  
		    axis.text.x = element_blank(),  
		    panel.grid = element_blank(),
		    legend.position = "right"
		  )

		p_trend <- ggplot(data, aes(x = x, y = gene_expression)) +
		  stat_summary_bin(
		    fun = "mean",
		    geom = "line", 
		    color = "red",
		    linewidth = 1,
		    bins = 30  
		  ) +
		  labs(y = "Mean Expression", x = "Spatial X-coordinate") +
		  theme_minimal(base_size = 12) +  
		  theme(
		    axis.title = element_text(size = 12, face = "bold", color = "black"),
		    axis.text.x = element_text(
		      angle = 0, 
		      hjust = 1, 
		      vjust = 1, 
		      size = 20,
		      color = "black"
		    ),
		    axis.text.y = element_text(size = 20, color = "black"),
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

		pdf(paste('./dn_overla_mou/ST_geneEXP' , gene_name , speci , 'chip', chip_id ,"_dn_overlay_0729.pdf",sep='_'),w=8,h=6)
		print(p_combined)
		dev.off()
	}
}

