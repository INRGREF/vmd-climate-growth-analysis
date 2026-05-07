# VMD Climate–Growth Analysis Pipeline

[![DOI](https://zenodo.org/badge/1216983686.svg)](https://doi.org/10.5281/zenodo.20072373)

R scripts implementing Variational Mode Decomposition (VMD) and multi-scale statistical analysis of climate–growth relationships in *Pinus halepensis*. The workflow combines signal decomposition, cycle extraction, PCA reduction, and PCR bootstrap regression to quantify scale-dependent climate sensitivity of tree growth.

---

# 🌲 VMD Climate–Growth Analysis (Pinus halepensis)

## 📌 Overview

This repository provides a complete and reproducible R workflow for multi-scale dendroclimatological analysis based on Variational Mode Decomposition (VMD).

The workflow decomposes climate (temperature and precipitation) and tree-ring growth signals into intrinsic oscillatory modes, extracts shared climatic cycles, and quantifies their relationships with tree growth using Principal Component Analysis (PCA) and PCR bootstrap regression.

The approach allows the identification of climate–growth relationships across multiple temporal scales (interannual to decadal to multi-decadal).

---

## 🎯 Scientific objectives

- Decompose climate and growth signals using VMD
- Extract multi-scale climatic cycles
- Transform monthly climate into biological year structure
- Perform PCA on VMD-derived climate cycles
- Quantify climate–growth relationships using PCR bootstrap
- Assess statistical significance via confidence intervals
- Identify scale-dependent climatic drivers of tree growth

---

## 📂 Repository structure

raw_data/
scripts/
├── 01_preprocessing.R
├── 02_vmd_pipeline.R
├── 03_pca_analysis.R
├── 04_pcr_climate_growth.R
├── 05_figures_heatmap.R
├── 06_variance_analysis.R
results/
├── VMD_temperature_modes/
│   ├── modes/
│   ├── Cycle_analysis/
│       └── Grouped_cycles
        └── Merged_cycles
           └──Biological_year
│
├── VMD_climate_precipitation/
│    ├── modes/
│   ├── Cycle_analysis/
│       └── Grouped_cycles
        └── Merged_cycles
           └──Biological_year/
│              
├── VMD_growth/
│   ├── Cycle_analysis/
│   ├── Grouped_cycles/
        └──Mean_chronologies
├── figures/
├── data/

## Requirements
R version
R ≥ 4.0
Required packages
VMDecomp
dplyr
tidyr
FactoMineR
ggplot2
dplR

## Workflow execution

Run scripts sequentially:

VMD decomposition (climate + growth)
PCA on VMD climate cycles
PCR bootstrap climate–growth analysis
Visualization (heatmaps)
Variance contribution analysis

## Methodological workflow
VMD decomposition
Monthly climate and tree-ring series decomposed into intrinsic modes (IMFs)
Optimal number of modes estimated adaptively
Cycle extraction
Periodicity estimated from instantaneous frequencies
Standardized cycle structure stored in Cycle_analysis/
Biological year transformation
Climate variables reorganized (Oct → Sep)
PCA reduction
PCA applied on VMD climate cycles
Dimensionality reduction using eigenvalue threshold criterion
PCR bootstrap regression
Tree growth modeled as function of climate principal components
Bootstrap (1000 iterations) used for uncertainty estimation
Back-transformation to monthly climate space
Significance assessment
Confidence intervals (95%) used to evaluate stability of relationships

## Outputs
VMD decomposed climate and growth modes
Climate cycle periodicities
PCA scores and loadings per cycle
Climate–growth regression coefficients
Bootstrap confidence intervals
Heatmaps of scale-dependent relationships
Variance contribution of modes and cycles

## Data availability

Climate and tree-ring datasets are available from the corresponding author upon reasonable request.

## Citation

title: "Beyond Classical Response Functions: Variational Mode Decomposition Reveals Frequency-Specific Climate Controls on Pinus halepensis Growth in Semi-Arid Tunisia"
If you use this repository, please cite:

Khorchani, A., Bachtobji, B., Aouinti, H., & Khaldi, A. (2026).  
VMD Climate–Growth Analysis (Version 1.0.1) [Computer software].  
Zenodo. https://doi.org/10.5281/zenodo.20072373

status: "manuscript submitted to Dendrochronologia"
## License

This project is licensed under the MIT License.

## Reproducibility statement

This workflow ensures full reproducibility from raw climatic and dendrochronological data to final statistical inference. All intermediate outputs are explicitly saved to enable stepwise validation and independent replication.
