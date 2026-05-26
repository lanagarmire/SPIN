# SPIN: Spatial Predictive Integration Network

**SPIN (Spatial Predictive Integration Network)** is a modular framework for population-scale clinical prediction and biological discovery from spatial omics networks. SPIN converts each subject’s spatial transcriptomic or proteomic profile into a subject-specific molecular network, links these networks to clinical phenotypes using Bayesian scalar-on-network regression with manifold learning (**BSNMani**), and provides interpretable downstream visualizations including subnetworks, hub genes/protein markers, pathway enrichment, and patient-level loading heatmaps.

This repository contains the code used for the manuscript:

> **Population-Scale Integration of Spatial Omics Networks for Clinical Prediction and Biological Discovery by SPIN**

---

## Overview

SPIN consists of three major modules:

1. **Subject-specific spatial network construction**  
   Constructs subject-level spatial co-expression or co-abundance networks from spatial omics data, such as:
   - gene–gene co-expression networks from MERFISH spatial transcriptomics data
   - marker–marker co-abundance networks from IMC spatial proteomics data

2. **Molecular-to-clinical prediction using BSNMani**  
   Uses Bayesian scalar-on-network regression with manifold learning to decompose each subject-level network into:
   - population-shared latent subnetworks
   - subject-specific subnetwork loading factors
   - clinical regression coefficients associated with patient phenotypes

3. **Visualization and biological interpretation**  
   Generates interpretable outputs, including:
   - subnetwork heatmaps
   - hub gene/protein marker ranking plots
   - gene–gene or marker–marker correlation network visualizations
   - pathway enrichment plots
   - integrated subnetwork–gene/marker–pathway networks
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
│   ├── benchmarking/
│   └── visualization/
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
│   ├── oligodendrocyte_networks/
│   ├── bsnmani_outputs/
│   └── enrichment_results/
├── figures/
│   ├── Figure1/
│   ├── Figure2/
│   ├── Figure3/
│   ├── Figure4/
│   └── Figure5/
└── supplementary_tables/
