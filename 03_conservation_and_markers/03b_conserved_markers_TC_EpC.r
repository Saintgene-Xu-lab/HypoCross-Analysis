library(spatstat)
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(cowplot)
library(RColorBrewer)
library(CoGAPS)
library(projectR)
library(pheatmap)
library(viridis)
library(patchwork)
library(purrr)

mac_st <- readRDS('./macaque_ST_transfer_result_20250407.RDS')
mar_st <- readRDS('./marmoset_ST_transfer_result_20250407.RDS')
mou_st <- readRDS('./mouse_ST_transfer_result_20250407.RDS')

dir <- './'
nn <- readRDS(paste(dir,"non_neu_male_69142_20250305.rds",sep=''))



####################################################################################


nn_TC <- subset(nn , subset= major_class == 'TC')
conserved_markers_TCa <- conserved_markers

conserved_markers_TCa <- FindConservedMarkers(
  object = nn_TC,
  ident.1 = "TCα",  
  ident.column = "annotation",
  grouping.var = "specie",
  only.pos = TRUE
)
setwd(dir)
saveRDS(conserved_markers ,'conserved_markers_TCa.RDS')
conserved_markers_TCb <- FindConservedMarkers(
  object = nn,
  ident.1 = "TCβ",  
  ident.column = "annotation",
  grouping.var = "specie",
  only.pos = TRUE
)
saveRDS(conserved_markers_TCb ,'conserved_markers_TCb.RDS')


conserved_markers <- conserved_markers_TCb


tc_a <- conserved_markers_TCa %>%
  mutate(gene = rownames(.)) %>%
  select(
    gene,
    mouse_log2FC = mouse_avg_log2FC,
    macaque_log2FC = macaque_avg_log2FC,
    marmoset_log2FC = marmoset_avg_log2FC,
    mouse_padj = mouse_p_val_adj,
    macaque_padj = macaque_p_val_adj,
    marmoset_padj = marmoset_p_val_adj
  ) %>%
  filter(mouse_padj < 0.05)  

tc_b <- conserved_markers_TCb %>%
  mutate(gene = rownames(.)) %>%
  select(
    gene,
    mouse_log2FC = mouse_avg_log2FC,
    macaque_log2FC = macaque_avg_log2FC,
    marmoset_log2FC = marmoset_avg_log2FC,
    mouse_padj = mouse_p_val_adj,
    macaque_padj = macaque_p_val_adj,
    marmoset_padj = marmoset_p_val_adj
  ) %>%
  filter(mouse_padj < 0.05)  

tc_a <- tc_a %>%
  mutate(
    direction_mouse_macaque = sign(mouse_log2FC) == sign(macaque_log2FC),
    direction_mouse_marmoset = sign(mouse_log2FC) == sign(marmoset_log2FC),
    direction_all = direction_mouse_macaque & direction_mouse_marmoset
  )

table(tc_a$direction_all)  

mouse_specific_genes <- tc_a %>%
  filter(
    (macaque_padj >= 0.05 | sign(mouse_log2FC) != sign(macaque_log2FC)) |
    (marmoset_padj >= 0.05 | sign(mouse_log2FC) != sign(marmoset_log2FC))
  )


primate_specific_genes <- conserved_markers_TCa[
  (conserved_markers_TCa$macaque_p_val_adj < 0.05 & conserved_markers_TCa$macaque_avg_log2FC > 0) &
  (conserved_markers_TCa$marmoset_p_val_adj < 0.05 & conserved_markers_TCa$marmoset_avg_log2FC > 0) &
  (conserved_markers_TCa$mouse_p_val_adj >= 0.05 | abs(conserved_markers_TCa$mouse_avg_log2FC) < 0.1),  
]


primate_genes <- rownames(primate_specific_genes)
write.csv(primate_specific_genes, "primate_specific_genes.csv")


library(pheatmap)

heatmap_data <- primate_specific_genes[, c("mouse_avg_log2FC", "macaque_avg_log2FC", "marmoset_avg_log2FC")]
rownames(heatmap_data) <- rownames(primate_specific_genes)
names(heatmap_data) <- c('mouse', 'macaque','marmoset')
pheatmap(
  heatmap_data,
  cluster_rows = TRUE,
  show_rownames = TRUE, 
  fontsize_row = 10,
  fontsize_col =15,
  main = "Log2FC of Primate-Specific Genes"
)

ma_data <- data.frame(
  gene = rownames(primate_specific_genes),
  avg_expr = (primate_specific_genes$mouse_pct.1 + primate_specific_genes$macaque_pct.1 + primate_specific_genes$marmoset_pct.1) / 3,
  log2FC_primate = (primate_specific_genes$macaque_avg_log2FC + primate_specific_genes$marmoset_avg_log2FC) / 2
)

ggplot(ma_data, aes(x = avg_expr, y = log2FC_primate)) +
  geom_point(color = ifelse(ma_data$gene %in% primate_genes, "red", "gray")) +
  labs(x = "Average Expression (Pct.1)", y = "Log2FC (Primate Avg)", title = "MA Plot") +
  theme_minimal()

library(clusterProfiler)
library(org.Mm.eg.db)  

entrez_ids <- bitr(primate_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)


go_res <- enrichGO(
  gene = entrez_ids$ENTREZID,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP"  
)
dotplot(go_res, showCategory = 20)


#################################################################### TCb ######################################################################

primate_specific_genes <- conserved_markers_TCb[
  (conserved_markers_TCb$macaque_p_val_adj < 0.05 & conserved_markers_TCb$macaque_avg_log2FC > 0) &
  (conserved_markers_TCb$marmoset_p_val_adj < 0.05 & conserved_markers_TCb$marmoset_avg_log2FC > 0) &
  (conserved_markers_TCb$mouse_p_val_adj >= 0.05 | abs(conserved_markers_TCb$mouse_avg_log2FC) < 0.1),  # mouse 不显著
]

nrow(primate_specific_genes)  

primate_genes <- rownames(primate_specific_genes)
write.csv(primate_specific_genes, "primate_specific_genes.csv")


heatmap_data <- primate_specific_genes[, c("mouse_avg_log2FC", "macaque_avg_log2FC", "marmoset_avg_log2FC")]
rownames(heatmap_data) <- rownames(primate_specific_genes)
names(heatmap_data) <- c('mouse', 'macaque','marmoset')
pheatmap(
  heatmap_data,
  cluster_rows = TRUE,
    cluster_cols = FALSE,
  show_rownames = TRUE,  
  fontsize_row = 10,
  fontsize_col =15,
  main = "Log2FC of Primate-Specific Genes"
)

ma_data <- data.frame(
  gene = rownames(primate_specific_genes),
  avg_expr = (primate_specific_genes$mouse_pct.1 + primate_specific_genes$macaque_pct.1 + primate_specific_genes$marmoset_pct.1) / 3,
  log2FC_primate = (primate_specific_genes$macaque_avg_log2FC + primate_specific_genes$marmoset_avg_log2FC) / 2
)

ggplot(ma_data, aes(x = avg_expr, y = log2FC_primate)) +
  geom_point(color = ifelse(ma_data$gene %in% primate_genes, "red", "gray")) +
  labs(x = "Average Expression (Pct.1)", y = "Log2FC (Primate Avg)", title = "MA Plot") +
  theme_minimal()

library(clusterProfiler)
library(org.Mm.eg.db)  

entrez_ids <- bitr(primate_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

go_res <- enrichGO(
  gene = entrez_ids$ENTREZID,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP"  
)
dotplot(go_res, showCategory = 20)



