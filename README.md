# vmd-climate-growth-analysis
R scripts implementing Variational Mode Decomposition (VMD) for multi-scale analysis of climate and tree-growth relationships. Includes full workflow for signal decomposition, and climate–growth regression to ensure reproducible dendroclimatological analysis.

# VMD Climate–Growth Analysis (Pinus halepensis)

## 📌 Overview

This repository contains R scripts for performing a multi-scale climate–growth analysis using Variational Mode Decomposition (VMD) and Principal Component Analysis (PCA). The workflow is applied to *Pinus halepensis* tree-ring data and associated climate variables in Mediterranean environments.

The objective is to identify scale-dependent relationships between climate variability and radial tree growth.

---

## 🎯 Main features

- Variational Mode Decomposition (VMD) of climatic and growth signals  
- Principal Component Analysis (PCA) for dimensionality reduction  
- Multi-scale climate–growth regression analysis  
- Analysis of interannual to decadal climatic variability  
- Fully reproducible R workflow  

---

## 📂 Repository structure
├── scripts/
│ ├── 01_preprocessing.R
│ ├── 02_vmd_analysis.R
│ ├── 03_pca_analysis.R
│ ├── 04_regression_analysis.R
│ └── 05_figures.R
│
├── data/
│ ├── climate_data.csv
│ ├── tree_growth_data.csv
│
├── results/
├── figures/
├── docs/
├── README.md
├── LICENSE
└── .gitignore

---

## 💻 Requirements

R (≥ 4.0)

Required packages:
- stats  
- signal  
- ggplot2  
- dplyr  
- FactoMineR  

Optional reproducibility tool:
```r id="renv2"
install.packages("renv")
renv::init()
## How to run

Run scripts in order:
01_preprocessing.R
02_vmd_analysis.R
03_pca_analysis.R
04_regression_analysis.R
05_figures.R
## Data

Climate and tree-ring datasets are available from the corresponding author upon reasonable request.
## Citation

If you use this code, please cite:

Khorchani A. (Year). VMD-based multi-scale climate–growth analysis in Pinus halepensis. GitHub repository. DOI: [to be added]
---
## License

This project is licensed under the MIT License.

## Reproducibility

This repository follows open science principles to ensure reproducibility of all analyses from raw data processing to final outputs.
