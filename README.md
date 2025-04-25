# VertMolRate: Datasets and scripts for analyzing latitudinal gradient in molecular rates and testing evolutionary speed hypothesis

This repository contains scripts and datasets to replicate the analyses presented in:

T. Cai, Z. Wen, Z. Jiang and Y. Zhen, Distinct latitudinal patterns of molecular rates across vertebrates. Proc. Natl. Acad. Sci. U.S.A. (2025). DOI: 10.1073/pnas.2423386122.

## Contents

This repository is structured as follows:

### 1. DataFiles

   A file including input data used for analyses.

- trees: Phylogenetic trees used for PGLMMs and Fig.1.
- MolEvolRate:
  - MolRate_mtDNA.csv: Molecular rate and predictors for mitochondrial DNA (mtDNA).
  - MolRate_nuDNA.csv: Molecular rate and predictors for nuclear DNA (nuDNA).
- GISLayers: GIS layers used for used for visualizing global spatial patterns of molecular rates.
- sisters_family: Data on sister taxa pairs used to examine the relationship between diversification rates and molecular rates.
- SpatialJoinFiles: Contains data on latitude, longitude, species names, molecular rates in ecoregions or grid cells.
- BMR: Basal metabolic rates data of birds.
- calibrations.csv: Calibration pionts used for ploting Fig.1.

### 2. Scripts

   The R scripts used to conduct the analyses.

- source_functions.R: R script includes some functions used for analyses.
- test_ESH_mtDNA.R: R script testing the evolutionary speed hypothesis (ESH) using mtDNA molecular rates. The script is designed to replicate the analyses from both the main text and the supplementary materials.
- test_ESH_nuDNA.R: R script testing the evolutionary speed hypothesis (ESH) using nuDNA molecular rates. Similar to the mtDNA analysis script, it includes procedures for both the primary and supplementary analyses.

### 3. Outputs

   The results of the analyses, including figures and tables present in main text and appendix.
   
- MainFigures: Results of main figures
- Supplementary: Results of supplementary figures.

### 4. SI Datasets

   Supporting Information for Datasets S1–S28 in main text and appendix.

