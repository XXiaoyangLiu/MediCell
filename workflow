### Analysis Workflow Pipeline
## The entire analysis pipeline consists of 7 steps.

### Step 1 Cell proportion deconvolution estimation
## All input files for CIBERSORTx were prepared following the official guidelines hosted on the CIBERSORTx website. Instead of utilizing the full transcriptomic repertoire of single-cell RNA-seq profiles to build the signature reference panel, we suguest identified cell type marker genes using FindMarkers function implemented in Seurat.
## For ENIGMA deconvolution, both the single cell reference atlas and bulk RNA-seq expression matrices were normalized to TPM values.
## Script：M01_ENIGMA.R
##After deconvolution from ENIGMA, we recommend retaining exclusively protein coding genes and removing lowly expressed genes to mitigate technical noise.

### Step 2 Genotype preprocessing
## Purpose： Perform quality control and filtering on genotype data (SNP array as example) matched to bulk RNA-seq samples, and compute genotype dosages

### Step 3 Covariate selection
## Script：M02_DataProcess_covariance_fordb.R

### Step 4 TWAS model training
## Purpose：Train prediction models to generate input files for S-PrediXcan
## Script：M03_ModelTraining.R

### Step 5 LD score preparation from matched 1000 Genomes population as covariates for S-PrediXcan
## Script：M05_DataProcess_covariance_for_TWAS

### Step 6 spredixcan calculation
## Script：M06_TWAS.sh

### Step 7 Downstream statistical analysis
## Script：M07_ResultStatistic.R
