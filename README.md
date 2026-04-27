# Final-Project-of-Data-Science-for-Biology-and-Medicine
This is a collection of codes and results of a final project for BIEN3320 in HKUST lectured by Professor Jiguang Wang
# Integrative Immune Profiling of Germ Cell Tumors Using Bulk and Spatial Transcriptomics

**Course:** BIEN 3320 Data-Science-for-Biology-and-Medicine 
**Author:** Atticus Zhou
**Repository:** https://github.com/atticuszhou288-maker/Final-Project-of-Data-Science-for-Biology-and-Medicine.git

## Overview

This repository contains the figure outputs and analysis code for an integrative immune microenvironment analysis of germ cell tumors (GCT), combining bulk RNA-seq and 10x Visium HD spatial transcriptomics. The study characterizes immune heterogeneity across GCT histological subtypes and identifies spatial immune niches using independent molecular scoring methods.

## Datasets

- **Bulk RNA-seq**: 15 GCT samples across 5 histological subtypes (Germinoma, Yolk Sac Tumor, Immature Teratoma, Mature Teratoma, Mixed GCT). Samples include intracranial, ovarian, and testicular lesions. Clinical metadata includes age, sex, lesion site, and race.
- **Spatial transcriptomics (10x Visium HD)**: 2,975 spots from an intracranial GCT tissue section. Each spot is assigned a spatial domain (`SD-Germinoma`, `SD-YST`, `SD-TLS`, `SD-Stromal`, `SD-Neuronal`) based on cell-type proportion composition of the bin (spot), derived from reference-based deconvolution and spatial neighbor analysis.

## Data Availability

**These datasets are proprietary and not yet published.** If you are interested in accessing the data, please contact the course instructor or the data owner directly.

## Methods

- **Bulk RNA-seq**: Filtering (genes with ≥10 counts in ≥3 samples), normalization to log₂(CPM+1). Immune and stromal scores computed via the `tidyestimate` package (ESTIMATE algorithm). Immune cell enrichment scores computed via `GSVA::ssGSEA` using MSigDB C8 single-cell-derived gene sets.
- **Spatial transcriptomics**: Spot-level QC (min. 5 genes per spot, min. 10 spots per gene), library-size normalization to 10,000 counts, log₁₀₆ transformation. Immune cell scores computed using a custom pure-Python implementation of ssGSEA with literature-curated immune lineage markers (≥5 genes per cell type).
- **Clustering**: Unsupervised K-means clustering (K=3) on the ssGSEA scores to define immune types (Immune_High, Immune_Mid, Immune_Low), fully independent of the pre-defined spatial domain labels.
- **Visualization**: UMAP for dimensionality reduction; stacked bar charts for domain–cluster correspondence; bar plots with Kruskal-Wallis tests for cross-domain comparisons.

## Reproducibility

- **R analysis**: Run `Bulk-RNA-seq_Immune_Analysis.R` in RStudio (requires `tidyestimate`, `GSVA`, `msigdbr`, `ggplot2`, `ggpubr`, `rstatix`).
- **Python analysis**: Open `Spatial_Analysis.ipynb` in Google Colab and execute all cells sequentially. Required packages: `scanpy`, `pandas`, `numpy`, `matplotlib`, `seaborn`, `scikit-learn`.

## License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

## Acknowledgments

We thank the course instructor Professor Jiguang Wang and TAs for providing the processed GCT datasets and for guidance throughout this project. The ESTIMATE and ssGSEA algorithms are based on the original publications by Yoshihara et al. (2013) and Barbie et al. (2009), respectively.
