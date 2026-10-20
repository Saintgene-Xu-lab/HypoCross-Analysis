
############################

result <- readRDS('./mouse_CoGAPS_result.rds')
macaque_aucMat <- readRDS('./macaque_aucMat_CoGAPS.rds')
mouse_aucMat <- readRDS('./mouse_aucMat_CoGAPS.rds')
marmoset_aucMat <- readRDS('./marmoset_aucMat_CoGAPS.rds')

library(tidyverse)
df_auc <- bind_rows(
  as.data.frame(macaque_aucMat) %>% mutate(Species = "Macaque", Pattern = rownames(macaque_aucMat)),
  as.data.frame(mouse_aucMat) %>% mutate(Species = "Mouse", Pattern = rownames(mouse_aucMat)),
  as.data.frame(marmoset_aucMat) %>% mutate(Species = "Marmoset", Pattern = rownames(marmoset_aucMat)),
) %>% 
  pivot_longer(cols = -c(Species, Pattern), names_to = "CellType", values_to = "AUC")
head(df_auc)

pattern_to_plots <- unique(df_auc$Pattern)  

library(ggplot2)
library(viridis) 

for (pattern_to_plot in pattern_to_plots){
  plot_data <- df_auc %>% 
    filter(Pattern == pattern_to_plot)
  plot_data$CellType <- gsub("TCα", "TCa", plot_data$CellType)
  plot_data$CellType <- gsub("TCβ", "TCb", plot_data$CellType)

  head(plot_data)
  p <- ggplot(plot_data, aes(x = Species, y = CellType, fill = AUC)) +
    geom_tile(color = "white", linewidth = 0.5) +
   scale_fill_gradientn(
    colors = c("white", "lightpink", "red", "darkred"),
    values = scales::rescale(c(min(plot_data$AUC), quantile(plot_data$AUC, 0.5), 
                              quantile(plot_data$AUC, 0.9), max(plot_data$AUC)))
  )+
    labs(
      title = paste("AUCell Scores for", pattern_to_plot),
      x = "Species", y = "Cell Type", fill = "AUC Score"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  pdf(paste('./' , pattern_to_plot,'.pdf',sep=''),w=3,h=5)
  print(p)
  dev.off()
}