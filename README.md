# VertMolRate

Included in this data package are scripts and files to replicate analyses in:

Cai, T., Wen, Z., Jiang, Z. and Zhen, Y. 2024. Distinct latitudinal gradients of molecular rates and global diversity of vertebrates. 

Included are:

## 1. DataFiles

  A folder containing the input data used for the analyses.

- trees: A time tree used for all analyses.
- Molecular rate and predictors: VertMolRate.csv.
- GISLayers: Some GIS layers are used to present latitudinal gradient patterns of molecular rates across grids.
- sister_family: Sister pairs data used for examining relationships between diversification rates and molecular rates at species level.
- SpatialjoinFiles: Some GIS layers include species names in each grid and ecoregion.
- BMR: Basal metabolic rates and coresponding phylogeny of birds.
- F34F61: Molecular rates are estimated by F34 and F61 models in codeml.

## 2. Scripts

  Scripts to perform the analyses used in the main text and supplementary text, respectively.

- Rscript1.R: An R script for analyzing the variation in molecular rates across different species and taxonomic groups.
- Rscript2.R: An R script for analyzing latitudinal gradient in molecular rates at species and assemblage levels.
- Rscript3.R: An R script to predict molecular rates using multiple PGLMMs.
- Rscript4.R: An R script to examine relationships between diversification rates and molecular rates.



## 3. Outputs

- MainFigures: Results of the main figures produced by the R scripts.
- Supplementary: Results of the supplementary figures produced by the R scripts.
  


## 4. source_functions.R

- source_functions.R: An R script containing functions used for the analyses.
