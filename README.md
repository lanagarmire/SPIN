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
│   │   ├── common/
│   │   ├── SEA_AD/
│   │   └── IMC_BC/
│   ├── cross_validation/
│   │   ├── SEA_AD/
│   │   │   ├── main_analysis/
│   │   │   └── all_section_sensitivity/
│   │   └── IMC_BC/
│   ├── benchmarking/
│   │   ├── SEA_AD/
│   │   └── IMC_BC/
│   └── visualization/
│       ├── common/
│       ├── SEA_AD/
│       ├── IMC_BC/
│       └── celltype_SEAAD/
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

## Analysis Workflow

### Step 1. Preprocess spatial omics data

#### SEA-AD MERFISH preprocessing

The raw MERFISH transcript tables record detected transcripts. Preprocessing includes:

- discarding transcript records not assigned to segmented cells (`cell_id = -1`);
- counting detected transcripts for each gene–cell pair;
- converting the long-format transcript table into a gene-by-cell count matrix;
- encoding gene–cell pairs with no detected transcript as `0`;
- retaining these zero entries as non-detections rather than imputing nonzero values;
- screening the 180 raw feature labels for mitochondrial-gene, blank-control, and negative-probe patterns;
- confirming that no mitochondrial genes or negative-control probes were detected;
- removing the 40 blank-control barcodes;
- retaining all 140 SEA-AD panel genes for downstream spatial co-expression analysis;
- filtering cells with fewer than 20 detected genes;
- normalizing expression values using a scale factor of 10,000;
- matching expression data with spatial coordinates, donor metadata, and cell-type annotations.

Representative entry points:

```bash
Rscript scripts/preprocessing/SEA_AD/preprocess_SEAAD_MERFISH.R
Rscript scripts/preprocessing/SEA_AD/annotate_SEAAD_cells.R
```

#### IMC breast cancer preprocessing

IMC preprocessing includes:

- matching single-cell marker-intensity data with cell coordinates and patient metadata;
- applying sample- and cell-level quality control;
- retaining analysis-relevant protein markers;
- normalizing marker-intensity measurements;
- applying the complete-case requirements of the survival model;
- constructing patient-specific marker–marker co-abundance matrices.

Representative entry point:

```bash
Rscript scripts/preprocessing/IMC_BC/preprocess_IMC_BC.R
```

---

### Step 2. Construct subject-specific molecular networks

For the SEA-AD MERFISH analysis, SPIN evaluated:

- WGCNA, included as a non-spatial reference method;
- hdWGCNA;
- Smoothie;
- SpaceX.

For the IMC breast cancer analysis, SPIN evaluated:

- WGCNA;
- hdWGCNA;
- Smoothie.

Each method produces a subject-specific symmetric gene–gene or marker–marker network matrix for downstream BSNMani modeling.

Representative SEA-AD entry points:

```bash
Rscript scripts/network_construction/SEA_AD/WGCNA/run_WGCNA_SEAAD.R
Rscript scripts/network_construction/SEA_AD/hdWGCNA/run_hdWGCNA_SEAAD.R
Rscript scripts/network_construction/SEA_AD/SpaceX/run_SpaceX_SEAAD.R
Rscript scripts/network_construction/SEA_AD/Smoothie/separate_sample_list.R
```

The Smoothie workflow additionally uses:

```text
scripts/network_construction/SEA_AD/Smoothie/smoothie_input_preparation.ipynb
scripts/network_construction/SEA_AD/Smoothie/src/
```

The oligodendrocyte-specific network-construction script is:

```bash
Rscript scripts/network_construction/SEA_AD/celltype_oligodendrocyte/run_hdWGCNA_oligodendrocyte_SEAAD.R
```

---

### Step 3. Fit BSNMani models

BSNMani decomposes each subject-specific network into shared latent subnetworks and subject-specific loading factors and links these loading factors to the clinical outcome while adjusting for relevant covariates.

Shared BSNMani helper functions and compiled routines are organized under:

```text
scripts/bsnmani_modeling/common/
```

The IMC BSNMani workflow is organized under:

```text
scripts/bsnmani_modeling/IMC_BC/
```

Representative IMC entry point:

```bash
Rscript scripts/bsnmani_modeling/IMC_BC/BSNMani.R
```

For Slurm-based execution:

```bash
sbatch scripts/bsnmani_modeling/IMC_BC/run_BSNMani_IMC.sh
```

The SEA-AD fitting and evaluation workflow is integrated with the donor-level cross-validation scripts described below.

---

### Step 4. Evaluate candidate configurations using subject-level cross-validation

SEA-AD dementia classification was evaluated using donor-level leave-one-out cross-validation.

IMC breast cancer survival prediction was evaluated using 5-fold cross-validation.

For each fixed candidate configuration and each cross-validation split:

1. the model was fitted using the training donors or patients;
2. the fitted model generated predictions for the held-out SEA-AD donor or held-out IMC patient fold;
3. the held-out predictions were stored;
4. predictions were aggregated across all cross-validation splits to calculate performance.

Candidate network-construction methods and candidate values of the latent-subnetwork number, `q`, were compared using aggregated out-of-fold performance. Model selection was not performed using nested cross-validation.

Candidate values were:

- `q = 2, 3, 4, or 5` for the SEA-AD all-cell analysis;
- `q = 2, 3, 4, or 5` for the oligodendrocyte-specific SEA-AD analysis;
- `q = 2–8` for the IMC survival analysis.

Representative SEA-AD scripts:

```bash
Rscript scripts/cross_validation/SEA_AD/main_analysis/Part1_LOOCV_Evaluation.R
Rscript scripts/cross_validation/SEA_AD/main_analysis/Part2_LOOCV_Evaluation.R
```

For Slurm-based execution:

```bash
sbatch scripts/cross_validation/SEA_AD/main_analysis/submit_loocv.sh
```

Representative IMC evaluation scripts:

```bash
Rscript scripts/cross_validation/IMC_BC/Yiwen_Results_withBootstrap.R
Rscript scripts/cross_validation/IMC_BC/Check_Yiwen_Smoothie_Results.R
```

---

### Step 5. Benchmark against conventional models

#### SEA-AD classification

BSNMani was compared with:

- Lasso logistic regression;
- Elastic Net logistic regression.

These baselines used vectorized network-edge features and the same clinical covariates used in the corresponding BSNMani analysis.

Representative scripts:

```bash
Rscript scripts/benchmarking/SEA_AD/benchmark_lasso_elasticnet_SEAAD.R
Rscript scripts/benchmarking/SEA_AD/clinical_model_LOOCV_SEAAD.R
```

#### IMC survival prediction

SPIN was compared with a clinical-variable-only Cox proportional hazards model. Additional scripts compare the candidate network-construction approaches.

Representative scripts:

```bash
Rscript scripts/benchmarking/IMC_BC/Baseline.R
Rscript scripts/benchmarking/IMC_BC/BenchmarkSpatialCoexpr.R
```

---

### Step 6. Generate visualization and interpretation outputs

The visualization scripts generate:

- latent-subnetwork heatmaps;
- pathway-enrichment plots;
- hub-gene or hub-marker rankings;
- gene–gene or marker–marker network plots;
- integrated subnetwork–gene/marker–pathway networks;
- subject-level loading-factor heatmaps;
- Kaplan–Meier survival curves;
- network-method and prediction-model benchmarking figures.

Representative SEA-AD scripts:

```bash
Rscript scripts/visualization/SEA_AD/network_method_comparison/plot_four_method_heatmaps_SEAAD.R
Rscript scripts/visualization/SEA_AD/subnetwork_interpretation/cluster_SEAAD_subnetworks.R
Rscript scripts/visualization/SEA_AD/lambda_exploration.R
```

The four-method heatmap script compares WGCNA, hdWGCNA, SpaceX, and Smoothie for the same donor using a common row/column ordering and color scale.

Representative IMC visualization scripts are located under:

```text
scripts/visualization/IMC_BC/
```

Shared pathway and enrichment visualization utilities are located under:

```text
scripts/visualization/common/
```

---

## Main Analysis Summary

### SEA-AD MERFISH dementia classification

- **Cohort:** 27 SEA-AD donors
- **Modality:** MERFISH spatial transcriptomics
- **Outcome:** dementia status
- **Validation:** donor-level leave-one-out cross-validation
- **Selected network method:** hdWGCNA
- **Selected number of latent subnetworks:** `q = 3`
- **Selection basis:** best overall performance balance across candidate configurations
- **Accuracy:** 0.81
- **AUROC:** 0.72
- **AUPRC:** 0.73
- **MCC:** 0.62
- **Macro F1:** 0.76

### IMC breast cancer survival prediction

- **Source cohort:** 253 breast cancer patients
- **Modality:** imaging mass cytometry spatial proteomics
- **Outcome:** overall survival
- **Validation:** 5-fold cross-validation
- **Selected network method:** Smoothie
- **Selected number of latent subnetworks:** `q = 2`
- **C-index:** 0.776, reported as 0.78 after rounding
- **Risk stratification:** low-, intermediate-, and high-risk groups defined using risk-score tertiles
- **Log-rank test:** `p = 1.09 × 10⁻⁹`

### Oligodendrocyte-specific SEA-AD analysis

- **Cell type:** oligodendrocytes
- **Outcome:** dementia status
- **Network method:** hdWGCNA, fixed to remain consistent with the all-cell SEA-AD analysis
- **Candidate latent-subnetwork numbers:** `q = 2, 3, 4, or 5`
- **Selected number of latent subnetworks:** `q = 2`
- **Validation:** donor-level leave-one-out cross-validation
- **Accuracy:** 0.70
- **MCC:** 0.39
- **F1:** 0.62
- **Macro F1:** 0.68

---

## Sensitivity Analysis Using All Matched SEA-AD Sections

The main SEA-AD analysis uses one representative MERFISH section per donor. To assess sensitivity to representative-section selection, an additional analysis incorporates all matched MERFISH sections available for each donor.

Each section is preprocessed using the same pipeline. The section-level information is then summarized into one donor-level network before BSNMani modeling, ensuring that each donor contributes only one network to the prediction analysis and preventing section-level pseudoreplication.

Representative scripts:

```bash
Rscript scripts/cross_validation/SEA_AD/all_section_sensitivity/BSNMani_allslides_fit_donorLOOCV.R
Rscript scripts/cross_validation/SEA_AD/all_section_sensitivity/BSNMani_allslides_diag_predict_donorLOOCV.R
Rscript scripts/cross_validation/SEA_AD/all_section_sensitivity/collect_donorLOOCV_predictions.R
```

The all-section analysis is reported as a sensitivity analysis and does not replace the main representative-section results.

---

## Evaluation Metrics

### SEA-AD classification

The following metrics are calculated from aggregated out-of-fold predictions:

- Accuracy
- AUROC
- AUPRC
- Matthews correlation coefficient
- F1 score
- Macro F1 score
- Precision
- Recall
- Specificity

### IMC survival prediction

The following measures are reported:

- concordance index;
- Kaplan–Meier survival curves;
- low-, intermediate-, and high-risk stratification;
- log-rank test results.

Uncertainty was quantified using nonparametric bootstrapping with `B = 3,000` resamples of the aggregated held-out prediction–outcome pairs.

---

## Output Files

Typical output categories include:

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

These subdirectories describe the intended organization of generated outputs. Large intermediate objects and raw datasets should remain outside Git version control.

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

Please also cite the SPIN manuscript when its final citation becomes available:

```bibtex
@unpublished{liu_spin,
  title  = {Population-Scale Integration of Spatial Omics Networks for Clinical Prediction and Biological Discovery by SPIN},
  author = {Liu, Tong and Ko, Euiseong and Yang, Yiwen and Wu, Haowen and Unjitwattana, Thatchayut and Wang, Suyuan and Kang, Jian and Li, Yijun and Garmire, Lana X.},
  note   = {Manuscript under review}
}
```

---

## Contact

For questions, please contact:

- Lana X. Garmire: `lgarmire@uab.edu`
- Yijun Li: `liyijun@ds.dfci.harvard.edu`

For code-related questions, reproducibility issues, or bug reports, please open an issue in this GitHub repository.

---

## License

Please see the [`LICENSE`](LICENSE) file for licensing information.
