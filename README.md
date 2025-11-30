# Microbial Single-Cell Transcriptomics Links Gut Microbiota Functional States to Metabolic Remodeling in Mice

## Overview
This repository contains the analysis code and data processing pipelines associated with the paper:  
**“Microbial Single-Cell Transcriptomics Links Gut Microbiota Functional States to Metabolic Remodeling in Mice”**
## Repository Structure

### `scripts/`
Original scripts used in the study are organized within this directory.  

---

### 1. `Reference_genome_comparison_and_species_classification/`
This directory contains all scripts used for benchmarking the three reference genome resources and performing species identification using a Kraken2-based classifier. The included scripts cover:  

**FASTQ preprocessing for taxonomic classification**  
Read filtering and quality control, Barcode/UMI extraction  
**Species classification scripts**  
Batch classification of FASTQ files with the previously described Kraken2-based workflow  
Parsing, summarizing, and exporting genus-level and species-level reports  
**Reference genome comparison utilities**  
Scripts to evaluate mapping / assignment rates across the three reference sets  
Scripts to compare database content (genus/species overlap, taxonomic breadth, host specificity)  
FastANI-based genome similarity comparison between reference catalogs  
  
Each script contains detailed parameters, tool versions, and instructions for execution.  

---

### 2. `smClassify_and_downstream_analysis/`
This directory contains all downstream R scripts used to generate the figures in the manuscript. Each script corresponds directly to a specific figure or figure panel.  

#### Included scripts:

**Figure 2**  
smClassify pipeline scripts  
Community composition summaries  

**Figure 3**  
 Cell clustering and dimensionality reduction  
 Differential gene expression analysis  
 KEGG pathway enrichment  
 Functional-cluster distribution  
 Species proportion analysis  
 Subcluster characterization  

**Figure 4**  
 PLS-DA analysis  
 Differential metabolite analysis  
 Network construction and visualization  
 KEGG pathway mapping  
 Metabolite boxplots and dotplots (amino acids, bile acids, LCFA, PUL-related metabolites)  

**Figure 5**  
 Analysis of community-embedded functional states  
 Cross-species functional complementation within the gut microbiota  
  
**Figure 6**  
 Subpopulation analysis of *M. gordoncarteri*  
 Disrupted nitrogen- and stress-response programs linked to host nitrogen imbalance
