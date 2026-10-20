#!/usr/bin/env Rscript
# Standalone extraction of the current inte_area.plot.r entropy scatterplot.
# Run from any directory:
#   Rscript entropy.2026-10-13.r
# Optional positional arguments: input QS, palette CSV, output SVG.
# Requires: SeuratObject, qs, dplyr, ggplot2, ggrepel, svglite.
# Input QS must contain specie, clu6, and clu_2 ... clu_50 metadata.
# The QS is approximately 21 GB; allow sufficient RAM for deserialization.
# No upstream integration/clustering is rerun; existing assignments are inputs.
# Historical SVG differs from current source (includes HY and y ticks > ln(2));
# this script reproduces the CURRENT calculation, not verified historical values.

required_packages <- c("SeuratObject", "qs", "dplyr", "ggplot2", "ggrepel", "svglite")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
}

base_dir <- "/public/group/weiwu_group/majianyun/data/Hypothalamus/obj.all.spatial"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 3L) stop("Usage: Rscript entropy.2026-10-13.r [input.qs] [palette.csv] [output.svg]")
input_qs <- if (length(args) >= 1L) args[[1L]] else
  file.path(base_dir, "parcellation", "iseu.no_split.qs")
palette_csv <- if (length(args) >= 2L) args[[2L]] else
  file.path(base_dir, "banksy_colorv2.csv")
output_svg <- if (length(args) >= 3L) args[[3L]] else
  file.path(base_dir, "parcellation", "figures", "entropy.2026-10-13", "etp.v2.ba.svg")
output_csv <- paste0(tools::file_path_sans_ext(output_svg), ".regions.csv")
output_session <- paste0(tools::file_path_sans_ext(output_svg), ".sessionInfo.txt")
outputs <- c(output_svg, output_csv, output_session)
if (tolower(tools::file_ext(output_svg)) != "svg") stop("Output must have an .svg extension.")
if (!file.exists(input_qs) || !file.exists(palette_csv)) stop("Input QS or palette CSV does not exist.")
if (any(file.exists(outputs))) {
  stop("Refusing to overwrite existing outputs: ", paste(outputs[file.exists(outputs)], collapse = ", "))
}

message("Loading metadata from ", input_qs)
seu <- qs::qread(input_qs)
metadata <- seu[[]]
rm(seu)
invisible(gc())
cluster_columns <- paste0("clu_", 2:50)
required_columns <- c("specie", "clu6", cluster_columns)
if (!all(required_columns %in% colnames(metadata))) {
  stop("Missing metadata columns: ", paste(setdiff(required_columns, colnames(metadata)), collapse = ", "))
}
metadata <- metadata[, required_columns, drop = FALSE]
if (!nrow(metadata) || anyNA(metadata)) stop("Required metadata must be nonempty and contain no missing values.")
metadata$specie <- as.character(metadata$specie)
metadata$clu6 <- as.character(metadata$clu6)
metadata$clu6[metadata$clu6 == "HY1"] <- "HY"
species <- c("mouse", "macaque", "marmoset")
if (!setequal(unique(metadata$specie), species)) stop("Expected exactly mouse, macaque, and marmoset.")

# Entropy in nats. Positive pseudocounts make both probabilities strictly positive.
entropy <- function(v) {
  p <- v / sum(v)
  -sum(p * log(p))
}

# For each K, normalize pseudocount-adjusted counts across clusters within species.
# Use explicit species/cluster names to avoid factor-code indexing ambiguities.
# Accumulate scores instead of retaining two cell-by-49 matrices in memory.
mouse_nhp_sum <- numeric(nrow(metadata))
macaque_marmoset_sum <- numeric(nrow(metadata))
for (column in cluster_columns) {
  labels <- as.character(metadata[[column]])
  counts <- as.matrix(table(labels, factor(metadata$specie, levels = species))) + 1
  corrected <- sweep(counts, 2L, colSums(counts), FUN = "/")

  # Primate abundance is the mean of the two species-corrected abundances.
  # It is not the sum of raw macaque and marmoset cell counts.
  primate <- rowMeans(corrected[, c("macaque", "marmoset"), drop = FALSE])
  mouse_nhp <- apply(cbind(mouse = corrected[, "mouse"], primate = primate), 1L, entropy)
  macaque_marmoset <- apply(corrected[, c("macaque", "marmoset"), drop = FALSE], 1L, entropy)
  # entropy() internally normalizes the macaque/marmoset pair to sum to one.
  index <- match(labels, rownames(corrected))
  mouse_nhp_sum <- mouse_nhp_sum + mouse_nhp[index]
  macaque_marmoset_sum <- macaque_marmoset_sum + macaque_marmoset[index]
  message("Calculated ", column)
}
metadata$cell_mou_mon_entropy_mean <- mouse_nhp_sum / length(cluster_columns)
metadata$cell_mac_mar_entropy_mean <- macaque_marmoset_sum / length(cluster_columns)

# FIRST average entropy over K for each cell bin, THEN compute its ratio.
# Region NMER = log(median(cell-level ratios)); plotted axes are separate means.
metadata$ratio_mean <- metadata$cell_mac_mar_entropy_mean /
  metadata$cell_mou_mon_entropy_mean
regional <- dplyr::summarise(
  dplyr::group_by(metadata, clu6),
  `mouse-NHP entropy` = mean(cell_mou_mon_entropy_mean),
  `macaque-marmoset entropy` = mean(cell_mac_mar_entropy_mean),
  entropy_ratio = log(median(ratio_mean)),
  count = dplyr::n(),
  .groups = "drop"
)
regional$NMER <- regional$entropy_ratio
# HY remains in entropy count tables; exclude it only from regional plotting.
regional <- dplyr::filter(regional, clu6 != "HY")
regional <- dplyr::arrange(regional, dplyr::desc(NMER))
regional$clu6 <- factor(regional$clu6, levels = regional$clu6)
if (!nrow(regional)) stop("No annotated regions remain after excluding HY.")
values <- as.matrix(regional[, c("mouse-NHP entropy", "macaque-marmoset entropy")])
if (any(!is.finite(values)) || any(values < 0 | values > log(2) + 1e-12)) {
  stop("Entropy values must be finite and lie between zero and ln(2).")
}

palette <- utils::read.csv(palette_csv, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
if (!all(c("nucle", "color") %in% names(palette))) stop("Palette must contain nucle and color columns.")
if (anyNA(palette[, c("nucle", "color")]) || anyDuplicated(palette$nucle)) stop("Palette labels must be unique and nonmissing.")
color_map <- stats::setNames(palette$color, palette$nucle)
if (!all(as.character(regional$clu6) %in% names(color_map))) {
  stop("Missing colors: ", paste(setdiff(as.character(regional$clu6), names(color_map)), collapse = ", "))
}

# Same point/label settings as current source, plus a seed for reproducibility.
# Spatial offsets, companion bar plots, scaled entropy, and other analyses omitted.
plot <- ggplot2::ggplot(regional, ggplot2::aes(
  x = `mouse-NHP entropy`, y = `macaque-marmoset entropy`, label = clu6, colour = clu6
)) +
  ggplot2::geom_point(size = 4) +
  ggrepel::geom_text_repel(
    size = 5, force_pull = 30, min.segment.length = 0, max.overlaps = Inf,
    seed = 42, max.time = Inf, max.iter = 10000
  ) +
  ggplot2::scale_colour_manual(values = color_map) +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = "none")

# Explicit plot/device; omit undocumented fix_text_size argument from original.
if (!dir.exists(dirname(output_svg))) {
  if (!dir.create(dirname(output_svg), recursive = TRUE)) stop("Cannot create output directory.")
}
ggplot2::ggsave(output_svg, plot = plot, device = svglite::svglite,
                width = 4, height = 4, units = "in", bg = "white")
utils::write.csv(regional, output_csv, row.names = FALSE)
writeLines(c(
  paste("Input QS:", normalizePath(input_qs)),
  paste("Palette CSV:", normalizePath(palette_csv)),
  "Cluster numbers: 2:50; pseudocount: 1; logarithm: natural; label seed: 42",
  utils::capture.output(sessionInfo())
), output_session)
message("Saved SVG: ", output_svg)
message("Saved regional values: ", output_csv)
message("Saved environment details: ", output_session)