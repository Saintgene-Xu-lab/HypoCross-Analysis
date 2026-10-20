library(spatstat)
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(cowplot)
library(RColorBrewer)
library(CoGAPS)
# library(projectR)
library(pheatmap)
library(viridis)
library(patchwork)
library(purrr)

mac_st <- readRDS('./macaque_ST_transfer_result_20250407.RDS')
mar_st <- readRDS('./marmoset_ST_transfer_result_20250407.RDS')
mou_st <- readRDS('./mouse_ST_transfer_result_20250407.RDS')


speci = "mouse"


######## ==============================================================================
#              TC EpC ST expression (point Chips)
######## ==============================================================================


plot_enhanced_spatial <- function(seurat_obj, gene, target_celltype, species, output_dir, chips_to_plot = "all") {
  plot_gene <- ifelse(species %in% c("mac", "mar"), toupper(gene), gene)
  message("\n[Start] ", species, " | ", gene, " | ", target_celltype)

  # 数据准备
  spatial_data <- tryCatch({
    df <- data.frame(
      x = if (species %in% c("mac", "mar")) seurat_obj@meta.data$rx else -seurat_obj@meta.data$rx,
      y = -seurat_obj@meta.data$ry,
      celltype = seurat_obj@meta.data$celltype_pred,
      chip = seurat_obj$chip,
      expression = as.numeric(seurat_obj@assays$RNA@data[plot_gene, ])
    )

    # 数据验证
    if (nrow(df) == 0) stop("empty")
    message("|-- PASS")
    message("|   cell number: ", nrow(df))
    message("|   X: [", paste(round(range(df$x),2), collapse = ", "), "]")
    message("|   Y: [", paste(round(range(df$y),2), collapse = ", "), "]")

    df
  }, error = function(e) {
    warning("! Fail: ", e$message)
    return(NULL)
  })

  if (is.null(spatial_data)) return(invisible(NULL))

  spatial_data <- spatial_data %>%
    mutate(
      display_group = case_when(
        celltype == target_celltype & expression > 0 ~ "target_high",
        celltype == target_celltype ~ "target_low",
        TRUE ~ "background"
      )
    )

  if (!identical(chips_to_plot, "all")) {
    available_chips <- unique(spatial_data$chip)
    matching_chips <- chips_to_plot[chips_to_plot %in% available_chips]

    if (length(matching_chips) == 0) {
      warning("! no chip\n",
              "  use: ", paste(available_chips, collapse = ", "), "\n",
              "  request: ", paste(chips_to_plot, collapse = ", "))
      return(invisible(NULL))
    }
    spatial_data <- spatial_data %>% filter(chip %in% matching_chips)
    message("|-- filt chip: ", paste(matching_chips, collapse = ", "))
  }

  target_counts <- table(spatial_data$display_group)
  message("|-- cell group:")
  message("|   target_high: ", target_counts["target_high"] %||% 0)
  message("|   target_low: ", target_counts["target_low"] %||% 0)

  if (sum(spatial_data$display_group == "target_high") == 0) {
    warning("! no target")
    return(invisible(NULL))
  }

  point_size <- ifelse(nrow(spatial_data) > 5000, 0.3,
                      ifelse(nrow(spatial_data) > 1000, 0.5, 1.2))

  for (current_chip in unique(spatial_data$chip)) {
    chip_data <- spatial_data %>% filter(chip == current_chip)

    p <- ggplot() +
      # bacground
      geom_point(
        data = filter(chip_data, display_group == "background"),
        aes(x = x, y = y), color = "#F0F0F0", size = point_size, alpha = 0.3
      ) +
      # 低表达目标细胞
      geom_point(
        data = filter(chip_data, display_group == "target_low"),
        aes(x = x, y = y), color = "#FF9999", size = point_size * 1.2
      ) +
      geom_point(
        data = filter(chip_data, display_group == "target_high"),
        aes(x = x, y = y), color = "#E41A1C", size = point_size * 1.5
      ) +
      coord_fixed() +
      theme_void() +
      labs(title = paste(gene, "in", target_celltype)) +
      theme(plot.title = element_text(size = 8, hjust = 0.5))

    # 保存图片
    filename <- sprintf("%s_%s_in_%s_chip%s.png",
                       species, gene, target_celltype, current_chip)
    tryCatch({
      ggsave(
        file.path(output_dir, filename),
        plot = p, width = 6, height = 6, units = "in", dpi = 150
      )
      message("|-- save: ", filename)
    }, error = function(e) {
      warning("! save fail: ", filename, "\n  error: ", e$message)
    })
  }
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

run_all_combinations <- function(output_dir = "spatial_plots") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("make dir: ", output_dir)
  }

  for (species in names(species_objs)) {
    sobj <- species_objs[[species]]$obj
    chips <- species_objs[[species]]$chips_to_plot

    species_dir <- file.path(output_dir, species)
    if (!dir.exists(species_dir)) dir.create(species_dir)

    message("\n", rep("=", 40))
    message("process: ", species)
    message("chip: ", paste(chips, collapse = ", "))
    message(rep("=", 40))

    for (combination in analysis_combinations) {
      celltype <- combination$celltype
      genes <- combination$genes

      message("\n[cell type] ", celltype)

      for (gene in genes) {
        plot_enhanced_spatial(
          seurat_obj = sobj,
          gene = gene,
          target_celltype = celltype,
          species = species,
          output_dir = species_dir,
          chips_to_plot = chips
        )
      }
    }
  }
  message("\n Finish: ", normalizePath(output_dir))
}

species_objs <- list(
  mou = list(
    obj = mou_st,
    chips_to_plot = c("T323")  
  ),
  mac = list(
    obj = mac_st,
    chips_to_plot = c("T787","T785","T773")  
  ),
  mar = list(
    obj = mar_st,
    chips_to_plot = c("T015_T965","T008_T970","T012_T623")   
  )
)


analysis_combinations <- list(
  list(celltype = "TCα", genes = c("Crym", "Glis3")),
  list(celltype = "TCβ", genes = c("Fat3")),

  list(celltype = "MG3", genes = c("Cd74")),
  list(celltype = "AS2", genes = c("Clu")),
  list(celltype = "AS1", genes = c("Slc7a10","Itih3")),

  list(celltype = "EpC1", genes = c("Ccdc39", "Ccdc181")),
  list(celltype = "EpC2", genes = c("Ak8"))
)

run_all_combinations()

