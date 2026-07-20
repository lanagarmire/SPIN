# SPIN: Spatial Predictive Integration Network

**SPIN (Spatial Predictive Integration Network)** is a modular framework for population-scale clinical prediction and biological discovery from spatial omics networks. SPIN converts each subject’s spatial transcriptomic or proteomic profile into a subject-specific molecular network, links these structured networks to clinical phenotypes using **BSNMani** (Bayesian scalar-on-network regression with manifold learning), and generates biologically interpretable outputs, including phenotype-associated subnetworks, hub genes or protein markers, enriched pathways, and subject-specific subnetwork loading factors.

This repository supports the manuscript:

> **Population-Scale Integration of Spatial Omics Networks for Clinical Prediction and Biological Discovery by SPIN**

---

## Overview

SPIN consists of three major modules.

### 1. Subject-specific spatial network construction

SPIN constructs subject-level spatial co-expression or co-abundance networks from spatial omics measurements, including:

* gene–gene co-expression networks from MERFISH spatial transcriptomics;
* marker–marker co-abundance networks from imaging mass cytometry (IMC) spatial proteomics.

The network-construction module is modular rather than hard-coded to one package. Any method that produces a compatible subject-specific symmetric molecular-network matrix can, in principle, be incorporated upstream of BSNMani.

### 2. Molecular-to-clinical prediction using BSNMani

BSNMani decomposes each subject-specific network into:

* population-shared latent subnetworks;
* subject-specific subnetwork loading factors;
* outcome and covariate effects linking the network representation to a clinical phenotype.

The framework supports binary, continuous, and time-to-event outcomes.

### 3. Visualization and biological interpretation

SPIN generates interpretable outputs, including:

* latent-subnetwork heatmaps;
* pathway-enrichment bubble plots;
* hub-gene or hub-marker ranking plots;
* gene–gene or marker–marker network visualizations;
* integrated subnetwork–gene/marker–pathway networks;
* donor- or patient-level clinical and loading-factor heatmaps.

---

## Repository Structure

```text
SPIN/
├── README.md
├── LICENSE
├── data/
│   ├── README.md
│   └── example_data/
├── scripts/
│   ├── preprocessing/
│   ├── network_construction/
│   ├── bsnmani_modeling/
│   ├── cross_validation/
│   ├── benchmarking/
│   └── visualization/
├── results/
│   ├── SEA_AD/
│   ├── IMC_BC/
│   └── celltype_SEAAD/
├── figures/
│   ├── Figure1/
│   ├── Figure2/
│   ├── Figure3/
│   ├── Figure4/
│   └── Figure5/
├── supplementary/
│   ├── figures/
│   ├── tables/
│   └── notes/
└── environment/
    ├── environment.yml
    └── renv.lock
```

---

## Data Sources

### SEA-AD MERFISH dataset

The Seattle Alzheimer’s Disease Brain Cell Atlas (SEA-AD) MERFISH dataset is publicly available through the Allen Institute and the Registry of Open Data on AWS:

* [SEA-AD AWS Open Data Registry](https://registry.opendata.aws/allen-sea-ad-atlas/)
* [Allen Brain Map SEA-AD portal](https://portal.brain-map.org/explore/seattle-alzheimers-disease)

The MERFISH panel contains 140 targets, including 100 biological marker genes and 40 blank control targets. The revised analysis uses a corrected donor-matched cohort of 27 SEA-AD donors with matched MERFISH expression data, spatial coordinates, clinical metadata, and cell-type annotations.

For the main analysis, one representative MERFISH section was selected per donor based on the number of high-quality segmented cells, ensuring that each donor contributed one subject-specific network. An additional sensitivity analysis incorporated all available matched sections for each donor by summarizing them into a single donor-level network before BSNMani modeling.

### IMC breast cancer dataset

The imaging mass cytometry breast cancer dataset is available from the original publication and its associated Zenodo repository:

* [The Single-Cell Pathology Landscape of Breast Cancer](https://zenodo.org/records/3518284)

The source cohort contains 253 breast cancer patients with single-cell spatial proteomic and clinical information. The analysis scripts apply the quality-control and complete-case requirements appropriate to each downstream model.

---

## Installation

### Option 1: Conda environment

```bash
conda env create -f environment/environment.yml
conda activate spin
```

### Option 2: R environment using `renv`

```r
install.packages("renv")
renv::restore()
```

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
```

Additional dependencies are required for specific network-construction methods, including hdWGCNA, Smoothie, and SpaceX.

---

## Analysis Workflow

### Step 1. Preprocess spatial omics data

#### SEA-AD MERFISH preprocessing

The raw MERFISH transcript table records detected transcripts. Preprocessing includes:

* discarding transcript records not assigned to segmented cells (`cell_id = -1`);
* counting detected transcripts for each gene–cell pair;
* converting the long-format transcript table into a gene-by-cell count matrix;
* encoding gene–cell pairs with no detected transcript as `0`;
* retaining these entries as non-detections rather than imputing nonzero values;
* screening feature names for mitochondrial-gene, blank-control, and negative-probe patterns;
* removing the 40 blank control targets;
* confirming that no mitochondrial genes or negative-control probes were present in the panel;
* filtering cells with fewer than 20 detected genes;
* normalizing expression values using a scale factor of 10,000;
* matching expression data with spatial coordinates, donor metadata, and cell-type annotations.

Example:

```bash
Rscript scripts/preprocessing/preprocess_SEAAD_MERFISH.R
```

#### IMC breast cancer preprocessing

IMC preprocessing includes:

* matching single-cell marker-intensity data with cell coordinates and patient metadata;
* applying sample- and cell-level quality control;
* retaining analysis-relevant protein markers;
* normalizing marker-intensity measurements;
* applying the complete-case requirements of the survival model;
* constructing patient-specific marker–marker co-abundance matrices.

Example:

```bash
Rscript scripts/preprocessing/preprocess_IMC_BC.R
```

---

### Step 2. Construct subject-specific molecular networks

For the SEA-AD MERFISH analysis, SPIN evaluated:

* WGCNA, included as a non-spatial reference baseline;
* hdWGCNA;
* Smoothie;
* SpaceX.

For the IMC breast cancer analysis, SPIN evaluated:

* WGCNA;
* hdWGCNA;
* Smoothie.

Each method produces a subject-specific symmetric gene–gene or marker–marker network matrix for downstream BSNMani modeling.

Example:

```bash
Rscript scripts/network_construction/run_network_construction_SEAAD.R
Rscript scripts/network_construction/run_network_construction_IMC.R
```

---

### Step 3. Fit BSNMani models

BSNMani decomposes each subject-specific network into shared latent subnetworks and subject-specific loading factors and then links these loading factors to the clinical outcome while adjusting for relevant covariates.

Example:

```bash
Rscript scripts/bsnmani_modeling/run_BSNMani_SEAAD.R
Rscript scripts/bsnmani_modeling/run_BSNMani_IMC.R
```

---

### Step 4. Evaluate candidate configurations using subject-level cross-validation

SEA-AD dementia classification was evaluated using donor-level leave-one-out cross-validation.

IMC breast cancer survival prediction was evaluated using 5-fold cross-validation.

For each fixed candidate configuration and each cross-validation split:

1. the BSNMani model was fitted using the training donors or patients;
2. the fitted model generated predictions for the held-out donor or fold;
3. the held-out prediction was stored;
4. predictions were aggregated across all folds to calculate performance.

Candidate network-construction methods and candidate values of the latent-subnetwork number (q) were compared using their aggregated out-of-fold performance.

Candidate (q) values were:

* (q \in {2,3,4,5}) for the SEA-AD all-cell analysis;
* (q \in {2,3,4,5}) for the oligodendrocyte-specific SEA-AD analysis;
* (q \in {2,\ldots,8}) for the IMC survival analysis.

Example:

```bash
Rscript scripts/cross_validation/run_LOOCV_SEAAD.R
Rscript scripts/cross_validation/run_5foldCV_IMC.R
```

---

### Step 5. Benchmark against conventional models

#### SEA-AD classification

BSNMani was compared with:

* Lasso logistic regression;
* Elastic Net logistic regression.

These baselines used vectorized network-edge features and the same clinical covariates used in the corresponding BSNMani analysis.

```bash
Rscript scripts/benchmarking/run_lasso_elasticnet_SEAAD.R
```

#### IMC survival prediction

SPIN was compared with a clinical-variable-only Cox proportional hazards model.

```bash
Rscript scripts/benchmarking/run_clinical_only_Cox_IMC.R
```

---

### Step 6. Generate visualization and interpretation outputs

The visualization scripts generate:

* latent-subnetwork heatmaps;
* pathway-enrichment bubble plots;
* hub-gene or hub-marker rankings;
* gene–gene or marker–marker network plots;
* integrated subnetwork–gene/marker–pathway networks;
* subject-level loading-factor heatmaps;
* Kaplan–Meier survival curves;
* network-method and prediction-model benchmarking figures.

Example:

```bash
Rscript scripts/visualization/plot_SEAAD_subnetworks.R
Rscript scripts/visualization/plot_IMC_subnetworks.R
Rscript scripts/visualization/plot_celltype_oligodendrocyte_SEAAD.R
```

---

## Main Analysis Summary

### SEA-AD MERFISH dementia classification

* **Cohort:** 27 SEA-AD donors
* **Modality:** MERFISH spatial transcriptomics
* **Outcome:** dementia status
* **Validation:** donor-level leave-one-out cross-validation
* **Selected network method:** hdWGCNA
* **Selected number of latent subnetworks:** (q = 3)
* **Selection basis:** best overall performance balance across candidate configurations
* **Accuracy:** 0.81
* **AUROC:** 0.72
* **AUPRC:** 0.73
* **MCC:** 0.62
* **Macro F1:** 0.76

### IMC breast cancer survival prediction

* **Source cohort:** 253 breast cancer patients
* **Modality:** imaging mass cytometry spatial proteomics
* **Outcome:** overall survival
* **Validation:** 5-fold cross-validation
* **Selected network method:** Smoothie
* **Selected number of latent subnetworks:** (q = 2)
* **C-index:** 0.776, reported as 0.78 after rounding
* **Risk stratification:** low-, intermediate-, and high-risk groups defined using risk-score tertiles
* **Log-rank test:** (p = 1.09 \times 10^{-9})

### Oligodendrocyte-specific SEA-AD analysis

* **Cell type:** oligodendrocytes
* **Outcome:** dementia status
* **Network method:** hdWGCNA, fixed to remain consistent with the all-cell SEA-AD analysis
* **Candidate latent-subnetwork numbers:** (q \in {2,3,4,5})
* **Selected number of latent subnetworks:** (q = 2)
* **Validation:** donor-level leave-one-out cross-validation
* **Accuracy:** 0.70
* **MCC:** 0.39
* **F1:** 0.62
* **Macro F1:** 0.68

---

## Sensitivity Analysis Using All Matched SEA-AD Sections

The main SEA-AD analysis uses one representative MERFISH section per donor. To assess sensitivity to representative-section selection, an additional analysis incorporates all matched MERFISH sections available for each donor.

Each section is preprocessed using the same pipeline. The section-level information is then summarized into one donor-level network before BSNMani modeling, ensuring that each donor contributes only one network to the prediction analysis and preventing section-level pseudoreplication.

The all-section analysis is reported as a sensitivity analysis and does not replace the main representative-section results.

---

## Evaluation Metrics

### SEA-AD classification

The following metrics are calculated from aggregated out-of-fold predictions:

* Accuracy
* AUROC
* AUPRC
* Matthews correlation coefficient
* F1 score
* Macro F1 score
* Precision
* Recall
* Specificity

### IMC survival prediction

The following measures are reported:

* concordance index;
* Kaplan–Meier survival curves;
* low-, intermediate-, and high-risk stratification;
* log-rank test results.

Uncertainty is quantified using nonparametric bootstrapping with (B = 3{,}000) resamples of the aggregated held-out predictions.

---

## Output Files

Typical output files include:

```text
results/
├── SEA_AD/
│   ├── representative_section_analysis/
│   ├── all_section_sensitivity_analysis/
│   ├── network_matrices/
│   ├── bsnmani_outputs/
│   ├── benchmark_results/
│   └── enrichment_results/
├── IMC_BC/
│   ├── network_matrices/
│   ├── bsnmani_outputs/
│   ├── survival_results/
│   └── enrichment_results/
└── celltype_SEAAD/
    ├── oligodendrocyte_networks/
    ├── bsnmani_outputs/
    └── enrichment_results/
```

---

## Reproducibility Notes

* The SEA-AD analyses were rerun using the corrected donor-matched MERFISH cohort.
* All cross-validation splits were defined at the donor or patient level.
* For each fixed candidate configuration, the held-out donor or patient was not used to fit the model that generated that subject’s prediction.
* Candidate configurations were compared using performance calculated from aggregated out-of-fold predictions.
* This repository does not claim nested or fold-internal hyperparameter selection.
* Uncertainty was quantified using 3,000 nonparametric bootstrap resamples.
* The main SEA-AD analysis uses one representative section per donor.
* The all-section sensitivity analysis summarizes multiple matched sections into one donor-level network before modeling.
* Large raw datasets are not redistributed through this repository.
* Local file paths in the analysis scripts may need to be updated before execution.

---

## Data and Code Availability

* SEA-AD source data: [SEA-AD AWS Open Data Registry](https://registry.opendata.aws/allen-sea-ad-atlas/)
* SEA-AD exploration portal: [Allen Brain Map SEA-AD](https://portal.brain-map.org/explore/seattle-alzheimers-disease)
* IMC source data: [Zenodo record 3518284](https://zenodo.org/records/3518284)
* Analysis code: this GitHub repository

Large raw datasets should be downloaded from their original public repositories. Processed files and intermediate results included in this repository or associated archival records are described in the corresponding data and results directories.

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

Please also cite the SPIN manuscript when its final citation becomes available:

```bibtex
@article{spin_spatial_omics,
  title   = {Population-Scale Integration of Spatial Omics Networks for Clinical Prediction and Biological Discovery by SPIN},
  author  = {Author list},
  journal = {Journal},
  year    = {Year}
}
```

---

## Contact

For questions, please contact:

* Lana X. Garmire: `lgarmire@uab.edu`
* Yijun Li: `liyijun@ds.dfci.harvard.edu`

For code-related questions, reproducibility issues, or bug reports, please open an issue in this GitHub repository.

---

## License

Please see the [`LICENSE`](LICENSE) file for licensing information.
