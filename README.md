# codeversion261013

Analysis code for comparative snRNA-seq and spatial transcriptomics (Stereo-seq) of the hypothalamus across three species: mouse, macaque, and marmoset.

The codebase is organized into 14 modules that follow the analysis workflow and figure logic of the paper — from preprocessing and integration to conservation analysis, neuropeptide comparison, spatial parcellation, and cross-region integration (Fig. 1–5 and supplementary figures).

> Code organization (directory structure, file naming, and documentation) was completed with AI assistance.

---

## Repository layout

| Directory | Description |
|---|---|
| [`01_snRNA_preprocess/`](#01_snrna_preprocess) | snRNA-seq preprocessing and quality control |
| [`02_snRNA_integration_annotation/`](#02_snrna_integration_annotation) | Integration, clustering, and hierarchical cell-type annotation |
| [`03_conservation_and_markers/`](#03_conservation_and_markers) | Cross-species conservation and marker gene analysis |
| [`04_neuron_analysis_F2/`](#04_neuron_analysis_f2) | Neuronal differential abundance and subtype analysis (Fig. 2) |
| [`05_neuropeptides_F3/`](#05_neuropeptides_f3) | Neuropeptide expression and co-expression (Fig. 3) |
| [`06_nonneuron_scCoGAPS_F4/`](#06_nonneuron_sccogaps_f4) | scCoGAPS pattern analysis of non-neuronal cells (Fig. 4) |
| [`07_OD_pseudotime_F4/`](#07_od_pseudotime_f4) | Oligodendrocyte pseudotime analysis (Fig. 4) |
| [`08_ST_data_build/`](#08_st_data_build) | Spatial transcriptomics object construction |
| [`09_ST_parcellation_BANKSY/`](#09_st_parcellation_banksy) | BANKSY spatial parcellation and annotation comparison |
| [`10_spatial_domains_F5/`](#10_spatial_domains_f5) | Cross-species spatial domain identification (Fig. 5) |
| [`11_regional_entropy_diversity/`](#11_regional_entropy_diversity) | Regional diversity and species-mixing entropy |
| [`12_nuclei_correlation/`](#12_nuclei_correlation) | Nuclei / chip-level transcriptional correlations |
| [`13_crossregion_scVI/`](#13_crossregion_scvi) | Cross-region scVI integration and similarity metrics |
| [`14_supplementary_plots/`](#14_supplementary_plots) | Supplementary figures and QC visualization |

---

## Module overview

### 01_snRNA_preprocess

Per-library preprocessing of snRNA-seq data for all three species, including quality control, SCTransform normalization, and doublet removal (DoubletFinder). One script per species plus shared doublet-detection utilities.

- `01a_mouse_preprocess.r` — mouse
- `01b_marmoset_preprocess.r` — marmoset
- `01c_macaque_preprocess.r` — macaque
- `01d_doublet_lib.r` — doublet-removal helpers

### 02_snRNA_integration_annotation

Reading and batch-QC of multi-species snRNA data, benchmarking of integration methods (CCA / RPCA / Harmony / fastMNN), and hierarchical cell-type annotation (L1–L6). Includes non-neuronal cluster refinement, label-transfer preparation, neuronal subtype split and re-integration, and cluster-tree construction.

- `02a` reading and batch QC
- `02b` / `02f` integration benchmarks (all cells / neurons)
- `02d` non-neuronal cluster refinement
- `02e` label-transfer preparation
- `02g` / `02h` neuronal subtype split and integration
- `02i` L1–L6 hierarchical annotation
- `02j` / `02k` cluster-tree structure and plotting

### 03_conservation_and_markers

Cross-species cell-type conservation and marker analysis. MetaNeighbor provides conservation scores (AUROC); COSG and FindAllMarkers identify conserved and species-specific genes; results are validated against the HypoMap reference.

- `03a` MetaNeighbor conservation AUROC
- `03b` / `03c` conserved markers and primate-elevated genes
- `03d`–`03h` neuronal markers, gene categories, divergence UMAPs
- `03i` / `03j` COSG conserved/divergent genes and transcription factors

### 04_neuron_analysis_F2

Cross-species comparison of neuronal populations. Milo is used for differential abundance testing; the module also produces Fig. 2 subtype UMAPs, marker gene panels (EN7, Sst, etc.), and summary diagrams.

- `04a` / `04b` Milo differential abundance and logFC UMAP
- `04c`–`04e` neuronal subtype panels
- `04f` Fig. 2 summary diagram

### 05_neuropeptides_F3

Cross-species comparison of neuropeptide expression: filtering the neuropeptide gene set, defining expressing cell populations per species, quantifying expression percentages, and analyzing co-expression patterns (Fig. 3).

- `05a` neuropeptide gene filtering
- `05b` expressing-cell definitions
- `05c` expression percentage
- `05d` co-expression analysis

### 06_nonneuron_scCoGAPS_F4

scCoGAPS pattern decomposition of non-neuronal cells (glia and related types), projectR projection of identified patterns onto cross-species data, and supervised UMAP visualization (Fig. 4).

- `06a` scCoGAPS and projectR projection
- `06b` pattern feature plots
- `06c` supervised UMAP
- `06d` / `06e` pattern heatmap and non-neuronal UMAP

### 07_OD_pseudotime_F4

Monocle3 pseudotime analysis of oligodendrocyte lineage cells (OD/OPC), identifying pseudotime-dependent gene expression changes and validating them in spatial transcriptomics data (Fig. 4).

- `07a` Monocle object preparation
- `07b` pseudotime heatmap (100 bins)
- `07c`–`07e` spatial validation panels

### 08_ST_data_build

Conversion of Stereo-seq SAW outputs into Seurat objects at cell-bin, Bin50, and Bin100 resolutions for all three species, including mouse Allen CCF-registered data. Provides shell batch scripts and R object-generation scripts.

- `08a`–`08d` batch shell scripts
- `08e`–`08i` mouse / macaque object generation
- `08k`–`08m` marmoset object generation

### 09_ST_parcellation_BANKSY

Spatial parcellation with BANKSY, including clustering, refinement, and comparison against anatomical annotations and spatial domains. Also covers microglia state annotation (MG1–MG4).

- `09a` BANKSY parcellation and refinement
- `09b` / `09c` parcellation vs. annotation / domain comparison
- `09d` / `09e` microglia domain plots and state annotation

### 10_spatial_domains_F5

Cross-species spatial domain identification using ortholog gene sets, LeverageScore sketching, and RPCA integration, followed by clustering and comparison with nuclei annotations and BANKSY parcellation (Fig. 5).

- `10a` ortholog data preparation
- `10b` spatial domain identification
- `10c` domain analysis and plots
- `10d` domain vs. nuclei comparison

### 11_regional_entropy_diversity

Regional diversity and species-mixing analyses: unsupervised spatial clustering, cross-species mixing entropy, local cell-type diversity, and nuclei-level transcriptional correlation heatmaps.

- `11a` cross-species neighborhood / embedding library
- `11b` unsupervised regional clustering
- `11c` species-mixing entropy calculation
- `11d` regional diversity and correlation plots

### 12_nuclei_correlation

Chip- and nuclei-level transcriptional correlation analysis of spatial transcriptomics data, including within-species chip comparison and cross-species nuclei comparison.

- `12a` / `12b` within-species chip correlations
- `12c` / `12d` three-species merge and cross-species correlation heatmaps

### 13_crossregion_scVI

Cross-region, cross-species integration of spatial transcriptomics data (Python notebooks): data merging, scVI integration, KL divergence, and RSM / NCM / NDM similarity metrics.

- `13a_species_merge.ipynb` — multi-region data merge
- `13b_scVI_integration.ipynb` — scVI integration and similarity analysis

### 14_supplementary_plots

Supplementary figure generation: sequencing QC panels, routine and per-donor UMAPs, chip replicate correlations, and ST gene-density plots.

- `14a` SF1 QC plots
- `14b` / `14c` routine and per-donor UMAPs
- `14d` / `14e` chip replicate correlation and gene density

---

## Requirements

- **R** (v5.0.3): Seurat, harmony, COSG, Milo, Monocle 3, scCoGAPS, projectR, uwot, ClusterR, and related packages
- **Python** (≥ 3.x): scvi-tools, scanpy, anndata
- See `library()` / `import` statements in each script for the full dependency list

---

## Notes

Code organization and documentation were completed with AI assistance. Analysis logic and parameters follow the paper's Methods section. A detailed code–Methods correspondence is available in `99_check_report.md`.
