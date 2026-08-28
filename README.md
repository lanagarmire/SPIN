# SPIN: Spatial Predictive Integration Network

**SPIN (Spatial Predictive Integration Network)** is a modular framework for population-scale clinical prediction and biological discovery from spatial omics networks. SPIN converts each subject’s spatial transcriptomic or proteomic profile into a subject-specific molecular network, links these structured networks to clinical phenotypes using **BSNMani** (Bayesian scalar-on-network regression with manifold learning), and generates biologically interpretable outputs, including phenotype-associated subnetworks, hub genes or protein markers, enriched pathways, and subject-specific subnetwork loading factors.

This repository supports the manuscript:

> **Population-Scale Integration of Spatial Omics Networks for Clinical Prediction and Biological Discovery by SPIN**

---

## Overview

SPIN consists of three major modules.

### 1. Subject-specific spatial network construction

SPIN constructs subject-level spatial co-expression or co-abundance networks from spatial omics measurements, including:

- gene–gene co-expression networks from MERFISH spatial transcriptomics;
- marker–marker co-abundance networks from imaging mass cytometry (IMC) spatial proteomics.

The network-construction module is modular rather than hard-coded to one package. Any method that produces a compatible subject-specific symmetric molecular-network matrix can, in principle, be incorporated upstream of BSNMani.

### 2. Molecular-to-clinical prediction using BSNMani

BSNMani decomposes each subject-specific network into:

- population-shared latent subnetworks;
- subject-specific subnetwork loading factors;
- outcome and covariate effects linking the network representation to a clinical phenotype.

The framework supports binary, continuous, and time-to-event outcomes.

### 3. Visualization and biological interpretation

SPIN generates interpretable outputs, including:

- latent-subnetwork heatmaps;
- pathway-enrichment plots;
- hub-gene or hub-marker rankings;
- gene–gene or marker–marker network visualizations;
- integrated subnetwork–gene/marker–pathway networks;
- donor- or patient-level clinical and loading-factor heatmaps;
- network-construction and prediction-model benchmarking figures.

---

## Repository Structure

```text
SPIN/
├── README.md
├── LICENSE
├── data/
│   └── example_data/
├── scripts/
│   ├── preprocessing/
│   │   ├── SEA_AD/
│   │   └── IMC_BC/
│   ├── network_construction/
│   │   ├── SEA_AD/
│   │   │   ├── WGCNA/
│   │   │   ├── hdWGCNA/
│   │   │   ├── Smoothie/
│   │   │   ├── SpaceX/
│   │   │   └── celltype_oligodendrocyte/
│   │   └── IMC_BC/
│   ├── bsnmani_modeling/
│   │   ├── SEA_AD/
│   │   └── IMC_BC/
│   └── visualization/
│       ├── common/
│       ├── SEA_AD/
│       └── IMC_BC/
└── figures/
    ├── Figure1/
    ├── Figure2/
    ├── Figure3/
    ├── Figure4/
    └── Figure5/

```

The `results/`, `figures/`, `supplementary/`, and `environment/` directories are organized as destinations for analysis outputs, manuscript materials, and reproducibility files. Large raw datasets and generated intermediate objects are not stored in this Git repository.

---

## Data Sources

### SEA-AD MERFISH dataset

The Seattle Alzheimer’s Disease Brain Cell Atlas (SEA-AD) MERFISH dataset is publicly available through the Allen Institute and the Registry of Open Data on AWS:

- [SEA-AD AWS Open Data Registry](https://registry.opendata.aws/allen-sea-ad-atlas/)
- [Allen Brain Map SEA-AD portal](https://portal.brain-map.org/explore/seattle-alzheimers-disease)

The SEA-AD MERFISH panel targeted 140 genes and profiled more than 300,000 segmented cells at near-single-cell resolution. The raw MERFISH transcript tables contained 180 feature labels, comprising the 140 panel genes and 40 blank-control barcodes. Feature names were screened for mitochondrial-gene, blank-control, and negative-probe patterns. No mitochondrial genes or negative-control probes were detected, and the 40 blank-control barcodes were removed because they do not correspond to biologically meaningful gene-expression features. Consequently, all 140 panel genes were retained for spatial co-expression network construction.

The revised analysis used a corrected donor-matched cohort of 27 SEA-AD donors with matched MERFISH expression data, spatial coordinates, clinical metadata, and cell-type annotations.

For the main analysis, one representative MERFISH section was selected per donor based on the number of high-quality segmented cells, ensuring that each donor contributed one subject-specific network. A sensitivity analysis incorporated all available matched sections for each donor by summarizing the section-level information into a single donor-level network before BSNMani modeling.

### IMC breast cancer dataset

The imaging mass cytometry breast cancer dataset is available from the original publication and its associated Zenodo repository:

- [The Single-Cell Pathology Landscape of Breast Cancer](https://zenodo.org/records/3518284)

The source cohort contains 253 breast cancer patients with single-cell spatial proteomic and clinical information. Individual downstream analyses apply the quality-control and complete-case requirements appropriate to the corresponding model.

---

## Installation

The analysis uses both R and Python. Complete environment lockfiles are being standardized. Until those files are available, install the dependencies required by the individual scripts and network-construction methods.

Key R dependencies include:

```r
Giotto
data.table
dplyr
tidyr
DescTools
WGCNA
glmnet
survival
survminer
ggplot2
pheatmap
ComplexHeatmap
enrichR
igraph
ggraph
Rcpp
```

Additional R and Python dependencies are required for specific network-construction methods. hdWGCNA and SpaceX primarily use R, whereas the included Smoothie workflow contains both R preprocessing scripts and Python modules.

Some scripts were originally developed across multiple computing environments. Before execution, configure the input, output, and helper-file paths for the local system.

---

## Reproducibility Notes

- The SEA-AD analyses were rerun using the corrected donor-matched MERFISH cohort.
- Cross-validation splits were defined at the donor level for SEA-AD and at the patient level for IMC.
- For each fixed candidate configuration, held-out observations were not used to fit the model that generated their predictions.
- Candidate configurations were compared using performance calculated from aggregated out-of-fold predictions.
- Model selection was not performed using nested cross-validation.
- Uncertainty was quantified using 3,000 nonparametric bootstrap resamples of the aggregated held-out prediction–outcome pairs.
- The main SEA-AD analysis uses one representative section per donor.
- The all-section sensitivity analysis summarizes multiple matched sections into one donor-level network before modeling.
- Large raw datasets and generated intermediate objects are not redistributed through this Git repository.
- Some analysis scripts currently contain dataset- or computing-environment-specific paths. These paths must be configured for the user’s local environment before execution.

---

## Data and Code Availability

- SEA-AD source data: [SEA-AD AWS Open Data Registry](https://registry.opendata.aws/allen-sea-ad-atlas/)
- SEA-AD exploration portal: [Allen Brain Map SEA-AD](https://portal.brain-map.org/explore/seattle-alzheimers-disease)
- IMC source data: [Zenodo record 3518284](https://zenodo.org/records/3518284)
- Analysis code: this GitHub repository

Large raw datasets should be downloaded from their original public repositories. Processed files and intermediate results included in an associated archival record should be documented in the corresponding data and results directories after public release.

---

## Citation

When citing the statistical model underlying SPIN, please cite BSNMani:

```bibtex
@article{li2024bsnmani,
  title   = {BSNMani: Bayesian Scalar-on-Network Regression with Manifold Learning},
  author  = {Li, Yijun and Choi, Ki Sueng and Dunlop, Boadie W. and Craighead, Wade Edward and Mayberg, Helen S. and Garmire, Lana and Guo, Ying and Kang, Jian},
  journal = {arXiv preprint arXiv:2410.02965},
  year    = {2024}
}
```


## Contact

For questions, please contact:

- Lana X. Garmire: `lgarmire@uab.edu`
- Yijun Li: `liyijun@ds.dfci.harvard.edu`

For code-related questions, reproducibility issues, or bug reports, please open an issue in this GitHub repository.

---

## License

Please see the [`LICENSE`](LICENSE) file for licensing information.
