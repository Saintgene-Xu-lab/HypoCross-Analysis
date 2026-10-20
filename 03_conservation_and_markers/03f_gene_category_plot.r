########### public datas -----------------------------------------------------------------
if(T){
  load("regulatory_database.Rdata")
  load("gene_category_20230218.Rdata")
  peti <- readxl::read_xlsx("petitle_all.xlsx")
  regulatory_database$ion_factor
  gene_category$tf$mouse_tf_list
  gene_category$epi$epi_gene_list_mouse
  gene_category$lr$mouse_receptor
  gene_category$lr$mouse_ligand
  
  category_list <- list(#IEGs = iges,
    Ligand = gene_category$lr$mouse_ligand,
    Receptor = gene_category$lr$mouse_receptor,
    Peptide = unique(peti$ligand),
    IOC = regulatory_database$ion_factor,
    EPI = gene_category$epi$epi_gene_list_mouse,
    TF = gene_category$tf$mouse_tf_list)
  term <- c("Ligand","Receptor","Peptide","IOC","EPI","TF")
}
color <- c("#AB47BC","#F44336","#FF9800","#03A9F4","#227c9d","#17c3b2","#ffcb77","#ffc2b4","#4DD0E1","#fb8f67","#E0E0E0")
names(color) = c("pri","mac","mar","mou","TF","EPI","Ligand","Receptor","Peptide","IOC","Other")
load("divergence_genes.Rdata")
########### function --------------
category_genes_con <- function(data,category){
  a <- data.frame(genes = data)
  rownames(a) <- a$genes
  a$group <- "Other"
  if (length(intersect(a$genes,category$Ligand)) > 0) {
    a[intersect(a$genes,category_list$Ligand),]$group <- "Ligand"
  }
  if (length(intersect(a$genes,category$Receptor)) > 0) {
    a[intersect(a$genes,category_list$Receptor),]$group <- "Receptor"
  }
  if (length(intersect(a$genes,category$Peptide)) > 0) {
    a[intersect(a$genes,category_list$Peptide),]$group <- "Peptide"
  }
  if (length(intersect(a$genes,category$IOC)) > 0) {
    a[intersect(a$genes,category_list$IOC),]$group <- "IOC"
  }
  if (length(intersect(a$genes,category$EPI)) > 0) {
    a[intersect(a$genes,category_list$EPI),]$group <- "EPI"
  }
  if (length(intersect(a$genes,category$TF)) > 0) {
    a[intersect(a$genes,category_list$TF),]$group <- "TF"
  }
  return(a)
}


######### ----
######################## conservative features -----------------------------------------------------------------------------------
single_feature <- readRDS("deg_in_single_specie.rds")
single_feature_v2 <- lapply(single_feature, function(x){
  x <- x[x$avg_log2FC > 0.5 & x$p_val_adj < 0.05,]
})


cells <- names(table(single_feature_v2$macaque$cluster))
conservate_cells <- list()
for (i in 1:length(cells)) {
  conservate_cells[[cells[i]]] <- intersect(intersect(single_feature_v2[[1]][single_feature_v2[[1]]$cluster %in% cells[i],]$gene,
                                                      single_feature_v2[[2]][single_feature_v2[[2]]$cluster %in% cells[i],]$gene),single_feature_v2[[3]][single_feature_v2[[3]]$cluster %in% cells[i],]$gene)
}

conservate_gene_category <- category_genes_con(data = unique(as.character(unlist(conservate_cells))),category = category_list)

conservate_gene_category_fre1 <- data.frame(table(conservate_gene_category$group))
conservate_gene_category_fre1$prop <- round(conservate_gene_category_fre1$Freq / sum(conservate_gene_category_fre1$Freq),digits = 4) * 100

conservate_gene_category_fre2 <- data.frame(table(conservate_gene_category[conservate_gene_category$group %in% term,]$group))
conservate_gene_category_fre2$prop <- round(conservate_gene_category_fre2$Freq / sum(conservate_gene_category_fre2$Freq),digits = 4) * 100

# Add label position
conservate_gene_category_fre1 <- conservate_gene_category_fre1 %>%
  arrange(desc(prop)) %>% 
  mutate(lab.ypos = cumsum(prop) - 0.5*prop) # 

# Add label position
conservate_gene_category_fre2 <- conservate_gene_category_fre2 %>%
  arrange(desc(prop)) %>% 
  mutate(lab.ypos = cumsum(prop) - 0.5*prop) # 

path <- ""
p <- ggplot(conservate_gene_category_fre1, aes(x = "", y = prop, fill = Var1)) +
  geom_bar(width = 1, stat = "identity", color = "white")+
  scale_fill_manual(values = color)+
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = prop), color = "white")+
  theme_void()
p
ggsave(p,filename = paste0(path,"Conservative_feature_with_other.pdf"),width = 5,height = 5)


p2 <- ggplot(conservate_gene_category_fre2, aes(x = "", y = prop, fill = Var1)) +
  geom_bar(width = 1, stat = "identity", color = "white")+
  scale_fill_manual(values = color)+
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = prop), color = "white")+
  theme_void()
p2
ggsave(p2,filename = paste0(path,"Conservative_feature.pdf"),width = 5,height = 5)

######### --------  primate

my_gene_category <- category_genes_con(data = unique(primate_up_genes$dat),category = category_list)
my_gene_category_fre1 <- data.frame(table(my_gene_category$group))
my_gene_category_fre1$prop <- round(my_gene_category_fre1$Freq / sum(my_gene_category_fre1$Freq),digits = 4) * 100

my_gene_category_fre2 <- data.frame(table(my_gene_category[my_gene_category$group %in% term,]$group))
my_gene_category_fre2$prop <- round(my_gene_category_fre2$Freq / sum(my_gene_category_fre2$Freq),digits = 4) * 100

# Add label position
my_gene_category_fre1 <- my_gene_category_fre1 %>%
  arrange(desc(prop)) %>% 
  mutate(lab.ypos = cumsum(prop) - 0.5*prop) # 

# Add label position
my_gene_category_fre2 <- my_gene_category_fre2 %>%
  arrange(desc(prop)) %>% 
  mutate(lab.ypos = cumsum(prop) - 0.5*prop) # 

path <- ""
p <- ggplot(my_gene_category_fre1, aes(x = "", y = prop, fill = Var1)) +
  geom_bar(width = 1, stat = "identity", color = "white")+
  scale_fill_manual(values = color)+
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = prop), color = "white")+
  theme_void()
p
ggsave(p,filename = paste0(path,"primate_up_feature_with_other.pdf"),width = 5,height = 5)


p2 <- ggplot(my_gene_category_fre2, aes(x = "", y = prop, fill = Var1)) +
  geom_bar(width = 1, stat = "identity", color = "white")+
  scale_fill_manual(values = color)+
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = prop), color = "white")+
  theme_void()
p2
ggsave(p2,filename = paste0(path,"primate_up_feature.pdf"),width = 5,height = 5)

############ macaque
my_gene_category <- category_genes_con(data = unique(mac_up_genes$dat),category = category_list)
my_gene_category_fre1 <- data.frame(table(my_gene_category$group))
my_gene_category_fre1$prop <- round(my_gene_category_fre1$Freq / sum(my_gene_category_fre1$Freq),digits = 4) * 100

my_gene_category_fre2 <- data.frame(table(my_gene_category[my_gene_category$group %in% term,]$group))
my_gene_category_fre2$prop <- round(my_gene_category_fre2$Freq / sum(my_gene_category_fre2$Freq),digits = 4) * 100

# Add label position
my_gene_category_fre1 <- my_gene_category_fre1 %>%
  arrange(desc(prop)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop) # 

# Add label position
my_gene_category_fre2 <- my_gene_category_fre2 %>%
  arrange(desc(prop)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop) # 

path <- ""
p <- ggplot(my_gene_category_fre1, aes(x = "", y = prop, fill = Var1)) +
  geom_bar(width = 1, stat = "identity", color = "white")+
  scale_fill_manual(values = color)+
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = prop), color = "white")+
  theme_void()
p
ggsave(p,filename = paste0(path,"mac_up_feature_with_other.pdf"),width = 5,height = 5)


p2 <- ggplot(my_gene_category_fre2, aes(x = "", y = prop, fill = Var1)) +
  geom_bar(width = 1, stat = "identity", color = "white")+
  scale_fill_manual(values = color)+
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = prop), color = "white")+
  theme_void()
p2
ggsave(p2,filename = paste0(path,"mac_up_feature.pdf"),width = 5,height = 5)

#### mar
my_gene_category <- category_genes_con(data = unique(mar_up_genes$dat),category = category_list)
my_gene_category_fre1 <- data.frame(table(my_gene_category$group))
my_gene_category_fre1$prop <- round(my_gene_category_fre1$Freq / sum(my_gene_category_fre1$Freq),digits = 4) * 100

my_gene_category_fre2 <- data.frame(table(my_gene_category[my_gene_category$group %in% term,]$group))
my_gene_category_fre2$prop <- round(my_gene_category_fre2$Freq / sum(my_gene_category_fre2$Freq),digits = 4) * 100

# Add label position
my_gene_category_fre1 <- my_gene_category_fre1 %>%
  arrange(desc(prop)) %>% 
  mutate(lab.ypos = cumsum(prop) - 0.5*prop) # 

# Add label position
my_gene_category_fre2 <- my_gene_category_fre2 %>%
  arrange(desc(prop)) %>% 
  mutate(lab.ypos = cumsum(prop) - 0.5*prop) # 

path <- ""
p <- ggplot(my_gene_category_fre1, aes(x = "", y = prop, fill = Var1)) +
  geom_bar(width = 1, stat = "identity", color = "white")+
  scale_fill_manual(values = color)+
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = prop), color = "white")+
  theme_void()
p
ggsave(p,filename = paste0(path,"mar_up_feature_with_other.pdf"),width = 5,height = 5)


p2 <- ggplot(my_gene_category_fre2, aes(x = "", y = prop, fill = Var1)) +
  geom_bar(width = 1, stat = "identity", color = "white")+
  scale_fill_manual(values = color)+
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = prop), color = "white")+
  theme_void()
p2
ggsave(p2,filename = paste0(path,"mar_up_feature.pdf"),width = 5,height = 5)

##### mou
my_gene_category <- category_genes_con(data = unique(mou_up_genes$dat),category = category_list)
my_gene_category_fre1 <- data.frame(table(my_gene_category$group))
my_gene_category_fre1$prop <- round(my_gene_category_fre1$Freq / sum(my_gene_category_fre1$Freq),digits = 4) * 100

my_gene_category_fre2 <- data.frame(table(my_gene_category[my_gene_category$group %in% term,]$group))
my_gene_category_fre2$prop <- round(my_gene_category_fre2$Freq / sum(my_gene_category_fre2$Freq),digits = 4) * 100

# Add label position
my_gene_category_fre1 <- my_gene_category_fre1 %>%
  arrange(desc(prop)) %>% 
  mutate(lab.ypos = cumsum(prop) - 0.5*prop) # 

# Add label position
my_gene_category_fre2 <- my_gene_category_fre2 %>%
  arrange(desc(prop)) %>% 
  mutate(lab.ypos = cumsum(prop) - 0.5*prop) # 

path <- ""
p <- ggplot(my_gene_category_fre1, aes(x = "", y = prop, fill = Var1)) +
  geom_bar(width = 1, stat = "identity", color = "white")+
  scale_fill_manual(values = color)+
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = prop), color = "white")+
  theme_void()
p
ggsave(p,filename = paste0(path,"mou_up_feature_with_other.pdf"),width = 5,height = 5)


p2 <- ggplot(my_gene_category_fre2, aes(x = "", y = prop, fill = Var1)) +
  geom_bar(width = 1, stat = "identity", color = "white")+
  scale_fill_manual(values = color)+
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = prop), color = "white")+
  theme_void()
p2
ggsave(p2,filename = paste0(path,"mou_up_feature.pdf"),width = 5,height = 5)









