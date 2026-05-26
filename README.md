# SPIN: Spatial Predictive Integration Network

**SPIN (Spatial Predictive Integration Network)** is a modular framework for population-scale clinical prediction and biological discovery from spatial omics networks. SPIN converts each subject's spatial transcriptomic or proteomic profile into a subject-specific molecular network, links these networks to clinical phenotypes using **BSNMani** (Bayesian scalar-on-network regression with manifold learning), and provides interpretable visualizations of subnetworks, hub genes/protein markers, enriched pathways, and patient-level loading factors.

This repository supports the manuscript:

> **Population-Scale Integration of Spatial Omics Networks for Clinical Prediction and Biological Discovery by SPIN**

---

## Overview

SPIN contains three major modules:

1. **Subject-specific spatial network construction**  
   Constructs subject-level spatial co-expression or co-abundance networks from spatial omics data, including:
   - gene-gene co-expression networks from MERFISH spatial transcriptomics data
   - marker-marker co-abundance networks from IMC spatial proteomics data

2. **Molecular-to-clinical prediction using BSNMani**  
   Uses Bayesian scalar-on-network regression with manifold learning to decompose each subject-specific network into:
   - population-shared latent subnetworks
   - subject-specific subnetwork loading factors
   - clinical regression coefficients associated with patient phenotypes

3. **Visualization and biological interpretation**  
   Generates interpretable outputs, including:
   - subnetwork heatmaps
   - pathway enrichment bubble plots
   - hub gene/protein marker ranking plots
   - gene-gene or marker-marker correlation network visualizations
   - integrated subnetwork-gene/marker-pathway networks
   - patient-level clinical/loading heatmaps

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
├── supplementary_tables/
└── environment/
    ├── environment.yml
    └── renv.lock
```

Please update the folder names above if your final repository uses a different organization.

---

## Data Sources

### SEA-AD MERFISH dataset

The SEA-AD MERFISH middle temporal gyrus dataset is publicly available through the Allen Brain Map / SEA-AD data portal:

- SEA-AD data portal: <https://registry.opendata.aws/allen-sea-ad-atlas/>
- Allen Brain Map portal: <https://portal.brain-map.org/explore/seattle-alzheimers-disease>

In the revised analysis, we used the corrected donor-matched SEA-AD MERFISH cohort. Only donors with matched MERFISH expression data, spatial coordinates, clinical metadata, and cell-type annotations were retained.

### IMC breast cancer dataset

The Imaging Mass Cytometry breast cancer dataset is available from the original publication and associated Zenodo repository:

- Zenodo record: <https://zenodo.org/records/3518284>

After quality control, the analysis used 253 patients with matched spatial proteomic measurements and clinical survival annotations.

---

## Installation

### Option 1: Conda environment

```bash
conda env create -f environment/environment.yml
conda activate spin
```

### Option 2: R package environment

If using `renv`:

```r
install.packages("renv")
renv::restore()
```

Key R packages used in the analysis include:

```r
Giotto
data.table
dplyr
tidyr
ggplot2
pheatmap
ComplexHeatmap
survival
survminer
glmnet
enrichR
igraph
ggraph
```

Additional dependencies may be required depending on the network-construction method.

---

## Analysis Workflow

### Step 1. Preprocess spatial omics data

SEA-AD MERFISH preprocessing includes:

- parsing specimen-level transcript tables
- aggregating transcript counts at the cell level
- matching donor metadata and cell-type annotations
- filtering low-quality cells
- removing blank control targets
- generating expression matrices and spatial coordinate files

Example:

```bash
Rscript scripts/preprocessing/preprocess_SEAAD_MERFISH.R
```

IMC breast cancer preprocessing includes:

- filtering samples with low cell counts
- removing patients with missing survival information
- selecting informative protein markers
- z-score normalizing protein marker intensities
- matching expression matrices with spatial coordinates

Example:

```bash
Rscript scripts/preprocessing/preprocess_IMC_BC.R
```

---

### Step 2. Construct subject-specific spatial networks

For the SEA-AD MERFISH analysis, SPIN evaluated:

- WGCNA
- SpaceX
- Smoothie
- hdWGCNA

For the IMC breast cancer analysis, SPIN evaluated:

- WGCNA
- Smoothie
- hdWGCNA

Example:

```bash
Rscript scripts/network_construction/run_network_construction_SEAAD.R
Rscript scripts/network_construction/run_network_construction_IMC.R
```

---

### Step 3. Run BSNMani modeling

BSNMani decomposes each subject-specific network into latent subnetworks and subject-specific loading factors.

Example:

```bash
Rscript scripts/bsnmani_modeling/run_BSNMani_SEAAD.R
Rscript scripts/bsnmani_modeling/run_BSNMani_IMC.R
```

---

### Step 4. Perform cross-validation and model selection

For SEA-AD dementia classification, subject-level leave-one-out cross-validation was used.

For IMC breast cancer survival prediction, 5-fold cross-validation was used.

Example:

```bash
Rscript scripts/cross_validation/run_LOOCV_SEAAD.R
Rscript scripts/cross_validation/run_5foldCV_IMC.R
```

Model selection was performed over candidate network-construction methods and candidate numbers of latent subnetworks `q`.

---

### Step 5. Benchmark against baseline models

For the SEA-AD classification task, BSNMani was compared against:

- Lasso logistic regression
- Elastic Net logistic regression

These baselines used vectorized lower-triangular network edges together with the same clinical covariates.

Example:

```bash
Rscript scripts/benchmarking/run_lasso_elasticnet_SEAAD.R
```

For the IMC survival task, SPIN was compared against a clinical-variable-only Cox proportional hazards model.

Example:

```bash
Rscript scripts/benchmarking/run_clinical_only_Cox_IMC.R
```

---

### Step 6. Generate visualization and interpretation outputs

Visualization scripts generate:

- subnetwork heatmaps
- pathway enrichment bubble plots
- hub gene/protein marker ranking plots
- gene-gene or marker-marker correlation network plots
- integrated subnetwork-gene/marker-pathway networks
- patient-level loading heatmaps

Example:

```bash
Rscript scripts/visualization/plot_SEAAD_subnetworks.R
Rscript scripts/visualization/plot_IMC_subnetworks.R
Rscript scripts/visualization/plot_celltype_oligodendrocyte_SEAAD.R
```

---


## Main Analysis Summary

### SEA-AD MERFISH dementia prediction

- Cohort: 27 SEA-AD donors
- Modality: MERFISH spatial transcriptomics
- Outcome: dementia status
- Final selected network method: hdWGCNA
- Final selected number of latent subnetworks: `q = 3`
- Evaluation: subject-level leave-one-out cross-validation
- Main performance: accuracy = 0.81

### IMC breast cancer survival prediction

- Cohort: 253 breast cancer patients
- Modality: Imaging Mass Cytometry spatial proteomics
- Outcome: overall survival
- Final selected network method: Smoothie
- Final selected number of latent subnetworks: `q = 2`
- Evaluation: 5-fold cross-validation
- Main performance: C-index = 0.78

### Oligodendrocyte-specific SEA-AD analysis

- Cell type: oligodendrocytes
- Outcome: dementia status
- Network method: hdWGCNA
- Final selected number of latent subnetworks: `q = 2`
- Evaluation: subject-level leave-one-out cross-validation

---

## Output Files

Typical output files include:

```text
results/
├── SEA_AD/
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

## Notes on Reproducibility

- The SEA-AD analyses were rerun using the corrected donor-matched MERFISH files released by the SEA-AD team.
- All model selection and evaluation were performed at the subject level.
- Held-out samples were not used for model fitting or parameter selection.
- Uncertainty was quantified using nonparametric bootstrapping with `B = 3,000` resamples over held-out prediction pairs.
- Exact file paths may need to be updated depending on where raw and processed datasets are stored locally.
- Large raw datasets are not stored in this repository. Please download raw data from the public data portals listed above.

---

## Citation

Please cite BSNMani:

```bibtex
@article{BSNMani2024,
  title   = {BSNMani: Bayesian Scalar-on-network Regression with Manifold Learning},
  author  = {Li, Yijun and others},
  year    = {2024},
  archivePrefix = {arXiv}
}
```

---

## Contact

For questions, please contact:

- Lana X. Garmire: lgarmire@uab.edu
- Yijun Li: liyijun@ds.dfci.harvard.edu

For code-related issues, please open an issue on this GitHub repository.

---

## License

Please see the `LICENSE` file for details.
