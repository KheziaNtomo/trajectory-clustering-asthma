# Trajectory-Based Clustering to Identify Drug Response Subgroups

Code for the analysis described in:

> **Trajectory based clustering to identify asthma subgroups following treatment with the selective CXCR2 antagonist, AZD5069**
>
> Asamoah K, Yang F, Adcock IM, Vuckovic D, Uddin M, Chung KF, Chadeau-Hyam M

## Overview

This repository implements a trajectory-based consensus clustering framework to identify patient response subgroups from longitudinal clinical trial data. The method uses Dynamic Time Warping (DTW) with stability-enhanced consensus clustering on multivariate time series.

## Scripts

Run in order:

| Script | Description |
|--------|-------------|
| `00_config.R` | Shared configuration: libraries, constants, data corrections |
| `01_data_processing.R` | Reads raw data and exports derived datasets |
| `02_data_preparation.R` | Loads datasets, merges longitudinal data, handles missingness |
| `03_trajectory_clustering.R` | DTW trajectory clustering with stability/consensus analysis |
| `04_visualisations.R` | Generates all manuscript figures and supplementary material |

## Method

1. **Data preparation**: Longitudinal clinical and biomarker data are merged, filtered for completeness, and scaled
2. **Trajectory clustering**: DTW distance with subsampling-based stability analysis (100 iterations, 50% subsamples), followed by PAM consensus clustering
3. **Cluster characterisation**: Biomarker log-fold changes, paired statistical tests, and demographic comparisons across identified subgroups

## R Packages

```r
install.packages(c(
  "dplyr", "tidyverse", "ggplot2", "tidyr",
  "haven",         # Data file reading
  "dtwclust",      # DTW clustering
  "cluster",       # PAM
  "sharp",         # Consensus scoring
  "clValid",       # Cluster validation
  "dtw",           # Dynamic Time Warping distances
  "Rtsne",         # t-SNE visualisation
  "pheatmap",      # Heatmaps
  "openxlsx",      # Excel export
  "reshape2",      # Data reshaping
  "tableone",      # Demographic tables
  "lme4"           # Mixed-effects models
))
```

## Outputs

- **Figure 1**: Longitudinal trajectories by cluster and treatment group
- **Figure 2**: Biomarker fold-change heatmap with significance
- **Supplementary Figure 1**: Neutrophil-to-lymphocyte ratio changes
- **Supplementary Tables 1–3**: Demographics and statistical tests

## Note

Data files are not included due to data governance restrictions. Subject identifiers have been anonymised for public release.
