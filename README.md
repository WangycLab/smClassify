# Microbial Single-Cell Transcriptomics Links Gut Microbiota Functional States to Metabolic Remodeling in Mice

## Overview
This repository contains the analysis code and data processing pipelines associated with the paper:  
**“Microbial Single-Cell Transcriptomics Links Gut Microbiota Functional States to Metabolic Remodeling in Mice”**
## Repository Structure

### `scripts/`
Original scripts used in the study are organized within this directory.

---

### 1. `mapping_and_species_identification/`
This folder contains all scripts used for processing FASTQ files and converting them into gene expression matrix. The workflow includes:

- Genome alignment using **STAR**  
- Filtering, UMI-aware quantification, and matrix generation  
- Species identification using Kraken2

Parameters and execution details are documented within each script.

---

### 2. `downstream_analysis_and_visualization/`
This folder contains all downstream R scripts used to generate the figures in the manuscript. Each script corresponds directly to a figure or figure panel.

#### Included scripts:

**Figure 2**  
Workflow of smClassify, community composition, comparison of classification performance by different databases.

**Figure 3**  
Clustering, differential gene expression, KEGG enrichment, functional cluster distribution, species proportion analysis, subcluster exploration.

**Figure 4**  
PLS-DA, differential metabolites, network analysis, KEGG pathway visualization, metabolite-related boxplots and dotplots (amino acids, bile acids, LCFA, PULs).

**Figure 5**  
Sample–species correlations, subcluster analysis of *Duncaniella* and *Parabacteroides*, multi-species correlation.

**Figure 6**  
Pseudotime analysis, subcluster analysis of *Muribaculum*, nitrogen/amino-acid metabolism analysis.

