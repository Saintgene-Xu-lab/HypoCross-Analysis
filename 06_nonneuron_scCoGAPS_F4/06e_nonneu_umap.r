library(RColorBrewer)
library(scales)
library(ggsci)
library(Seurat)
library(dplyr)


dir <- './'
nn <- readRDS(paste(dir,"non_neu_male_69142_20250305.rds",sep=''))


head(nn@meta.data,2)

nn@meta.data <- nn@meta.data %>%
  mutate(annotation = case_when(
    annotation == "TCα" ~ "TCa",
    TRUE ~ annotation  
  ))

nn@meta.data <- nn@meta.data %>%
  mutate(annotation = case_when(
    annotation == "TCβ" ~ "TCb",
    TRUE ~ annotation  
  ))

pdf('nn_reUMAP_1014_3spe_spli.pdf',h=10,w=20)
DimPlot(nn, group.by = "annotation",split.by = 'specie', label = T, reduction = 'umap', cols = cor ,pt.size = 0.1)
dev.off()
library(ggplot2)

p <- DimPlot(nn, group.by='annotation',split.by = 'specie',label=F,pt.size = 0.1,

        cols = cor)+

  theme(axis.text = element_text(size=8,colour = 'black'),

        axis.title = element_text(size=8,colour = 'black'),
        # axis.title = element_blank(),
        panel.background=element_rect(fill='transparent', color='black',size = 1),

        legend.key=element_blank(),legend.key.size  = unit(0.5,'cm'),

        legend.key.height = unit(0.1,'cm'),legend.key.width  = unit(0.5,'cm'),

        legend.text = element_text(size=7,colour = 'black'),

        plot.margin = unit(c(0,0,-0.1,-0.1),'cm'),
            plot.title = element_blank()  # 去掉标题
  )
ggsave("./nn_reUMAP_1014_3spe_spli_2.3.pdf",p,width = 42,height = 15,units = 'cm')




#### 
# color use



my_vector <-  unique(nn@meta.data$annotation)


sorted_vector <- sort(my_vector)

print(sorted_vector)

my_cor <- cor[1:23]

my_vec_cor <- data.frame(cell_types = my_vector, color_use = my_cor)
print(my_vec_cor)

write.xlsx(my_vec_cor, "./nn_color_use_20250609.xlsx")
