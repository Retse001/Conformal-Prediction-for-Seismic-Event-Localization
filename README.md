# Conformal-Prediction-for-Seismic-Event-Localization

This repository contains the complete scientific pipeline, synthetic simulation frameworks, real-world regional validation datasets, and R/Python implementations developed for my thesis project. The overarching goal of this research is to bridge distribution-free predictive inference (**Conformal Prediction**) with clustering algorithms and machine learning-driven seismological architectures.

---

## 📁 Repository Structure & GitHub Directory Map

### 📂 DIR 1: Numerical Simulations and Clustering Foundations
*Located in the `DIR 1` folder, this section evaluates conformal inference methodologies applied to 2D Gaussian Mixture Models (GMM 2D with 3 components), reproducing and expanding upon the framework by Nath, Hur, and Allen (2026) [arXiv:2604.03488v1].*

* **`CSFCLestrapolation.R`**: Replicates the *Conformal Sets for Cluster Labels* (CSFCL) methodology across four data-overlapping and noise levels ($\sigma^2 \in \{0.3, 0.5, 0.8, 1.2\}$) corresponding to *very_small*, *small*, *medium*, and *large* variance scenarios.
* **`ConformalCentroid.R`**: **Proposed Contribution** — Implements our proposed extension to quantify and bound the spatial geometric uncertainty of estimated GMM cluster centroids using conformal prediction blocks.

### 📂 DIR 2: Controlled Seismic Simulations Framework
*Located in the `DIR 2` folder, this block handles synthetic seismic event catalogs, phase picking link simulations, and baseline conformal validations across discrete spatial separation environments.*

* **Simulation Engines**:
  * **`1.1Sim3clAS.R`**: Synthetic catalog generation under **High Separation** structural constraints (3 clusters).
  * **`1.2MedSep3cl.R`**: Synthetic catalog generation under **Medium Separation** constraints (3 clusters with mild boundaries).
  * **`1.3FSovrap3cl.R`**: Synthetic catalog generation under **Strong Overlap** dense constraints (highly intersecting clusters).
* **Phase Association Block**:
  * **`2.Prob3cl.R`**: Applies the state-of-the-art *Earthquake Phase Association Using a Bayesian Gaussian Mixture Model* (**GaMMA**) by Zhu et al. (Stanford University) [Zenodo Core](https://doi.org/10.5281/zenodo.6271310).
* **Conformal Localization Closures**:
  * **`3.1CP.R`**, **`3.2CP.R`**, **`3.3CP.R`**: Methodological contribution scripts applying conformal processing to calculate temporal quantiles and generate mapping regions for the High, Medium, and Strong Overlap simulation sets respectively.

### 📂 DIR 3: Real-Data Regional Validations 
*Located in the `DIR 3` folder, this section scales the conformal pipeline to real-data.*

* **`Ridgecrest.R`**: Conformal-Prediction-for-Seismic-Event-Localization for the Ridgecrest, California.
* **`Chile.R`**: Conformal-Prediction-for-Seismic-Event-Localization for Chile.
* **`NorCal.R`**: Conformal-Prediction-for-Seismic-Event-Localization for Northern California seismic network.


---

## 🚀 Getting Started

### Prerequisites
Ensure your local environment includes an active **R** installation equipped with the following core packages:
```R
install.packages(c("readr", "dplyr", "ggplot2", "ggforce", "patchwork", "MASS", "jsonlite", "purrr"))
