# HypoCross-Analysis
Analysis code for **"Cross-species transcriptomic atlases reveal extensive specialization in the ancient hypothalamus"** (mouse, marmoset, macaque snRNA-seq + Stereo-seq).

The repository is organized in **two layers**:

- **`preprocessing/` (pipeline)** — shared upstream processing, ordered by stage. These build the objects that all figures depend on (per-species QC → integration/annotation → spatial object generation → parcellation → cross-region integration).
- **`figures/`** — downstream, figure-specific analysis and plotting, organized by main figure and by panel.

Shared helper functions live in **`lib/`**. Large intermediate data objects are **not** in git; see `data_objects/README.md`.

> Code was contributed by multiple authors; upstream (`preprocessing`) and figure scripts follow slightly different styles (`.rds`/`future` vs `qs`/Seurat v5). Absolute paths and blanked `""` paths from the authors' clusters remain and must be set locally before running (see *Reproducing* below).

## Repository structure

```
HypoCross-Analysis/
├── preprocessing/
|    ├── 01_snRNA_preprocessing/           # per-species QC + SCT + doublet removal
│    │    ├── mouse_preprocess.r
│    │    ├── macaque_preprocess.r
│    │    └── marmoset_preprocess.r
│    ├── 02_snRNA_integration_annotation/  # unified taxonomy
│    │    ├── nonneuronal_refine_cluster.r  #   per-class integration & annotation 
│    │    ├── transfer_snRNA_preparation.r  #   prepare per-species objects for label transfer
│    │    └── snRNA_qc_and_exploration.r    #   (scratch QC — not a required step)
│    ├── 03_spatial_generation/            # build binned Stereo-seq objects (Bin50 / Bin100)
│    │    ├── run_gen_bin50.sh              #   entry points (drive the gen_*.r scripts)
│    │    ├── run_gen_bin100.sh
│    │    ├── run_gen_marmoset2.sh
│    │    ├── run_gen_allen_ccf.sh
│    │    ├── run_gen_allen_ccf_bin100.sh
│    │    └── gen_{mouse,macaque,marmoset,marmoset2}_bin50.r , gen_mouse_allen{,_bin100}.r
│    ├── 04_spatial_parcellation/          # spatial domains
│    │    ├── banksy_refine.r               #   BANKSY spatial-domain clustering
│    │    ├── area_integration_unsupervised.r
│    │    └── area_integration_plot.r
│    ├── 05_crossregion_scvi/              # cross-region / cross-species scVI integration (Fig 6)
│    │    ├── species_merge.ipynb
│    │    └── region_integration.ipynb      #   scVI integration + KL divergence + RSM/NDM
├── figures/
│   ├── fig1_atlas/
│   ├── fig2_neuronal/
│   ├── fig3_neuropeptide/
│   └── fig5_parcellation/
├── environment/                      # TODO: add sessionInfo / conda env / renv.lock
└── MANIFEST_reorganization.md        # full old→new mapping + cleanup decisions
```

## Figure-to-code map

Fill in the input-object column as you finalize; this table is what reviewers use to trace each result.

| Panel(s) | Script | Notes |
|---|---|---|
| Fig 1C (UMAPs) | `figures/fig1_atlas/plot_umap_major_class.R`, `plot_umap_per_species.R` | major-class + per-species embeddings |
| Fig 1C (dendrogram) | `figures/fig1_atlas/tree_1_structure.R`, `tree_2_plot.R` | circular subtype dendrogram |
| Fig 1D/E (concordance) | `figures/fig1_atlas/auc_cor_umap_labeltransfer.R` | AUROC + correlation, label-transfer concordance |
| Fig 1F (ST density) | `figures/fig1_atlas/plot_st_gene_density.R` | ARH marker density |
| Fig 1 / donor2 (fig S2) | `figures/fig1_atlas/plot_donor2_reproducibility.R` | donor-to-donor reproducibility |
| fig S2B (ST correlation) | `figures/fig1_atlas/st_cor_1…4*.R` | cross-species section correlation |
| fig S1 (QC) | `figures/fig1_atlas/figS1_qc_plot.R` | UMI/gene QC |
| Fig 2A (neuronal taxonomy) | `figures/fig2_neuronal/1…4_*.R`, `plot_umap.R` | integration → subcluster → hierarchical levels |
| Fig 2B (Milo) | `figures/fig2_neuronal/milo_1_get_cells.R`, `milo_2_fc_umap.R` | differential abundance |
| Fig 2C/D (EN7 / SST) | `figures/fig2_neuronal/5_plot_en7.R`, `6_plot_sst.R` | VMH SST+/SLC17A6+ population |
| Fig 2F/G (conserved/divergent) | `figures/fig2_neuronal/7_conserved_divergent_genes.R`, `8_*_tf.R`, `ontology_1…5_*.R` | Hypo-neuronal ontology |
| Fig 3A–F (neuropeptides) | `figures/fig3_neuropeptide/1…4_*.R` | peptide id / per-species / % / co-expression |
| Fig 5 (parcellation) | `04_spatial_parcellation/banksy_refine.r` + `figures/fig5_parcellation/1…6_*.R` | BANKSY domains + analysis/plots |
| Fig 5L–R (marmoset MG niche) | `figures/fig5_parcellation/7_plot_marmoset_mg_domain.R`, `8_mg_state_annotation.R` | Mar-MPOA / MG1 + GO |
| Fig 6 (cross-region divergence) | `05_crossregion_scvi/species_merge.ipynb`, `region_integration.ipynb` | KL divergence, RSM/NDM |

**Not yet in this collection (still being gathered):** Fig 4 (non-neuronal scCoGAPS panels — the objects are built in `02_snRNA_integration_annotation/nonneuronal_refine_cluster.r`), the ligand–receptor / long-range parts of Fig 3, and the plotting for parts of Fig 6.

## Reproducing

1. Set paths: pipeline scripts contain absolute cluster paths; figure scripts contain blanked `qread("")` / `readRDS(paste0(path, ...))`. Replace with your local paths (a future `config.R` with a single `DATA_DIR` is recommended).
2. Environment: add `environment/` (see below). Two R stacks are used (Seurat v4 + SCT for the pipeline; Seurat v5 for several figure scripts) plus a Python/scVI env for `05_*`.
3. Run order: `01 → 02 → 03 → 04 → 05`, then `figures/figN/` in numeric order.
