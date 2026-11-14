### Analysis Workflow Pipeline
## The entire analysis pipeline consists of 7 steps.

### Step 1 Cell proportion deconvolution estimation
## All input files for CIBERSORTx were prepared following the official guidelines hosted on the CIBERSORTx website. Instead of utilizing the full transcriptomic repertoire of single-cell RNA-seq profiles to build the signature reference panel, we suguest identified cell type marker genes using FindMarkers function implemented in Seurat.
## For ENIGMA deconvolution, both the single cell reference atlas and bulk RNA-seq expression matrices were normalized to TPM values.
## Script：M01_ENIGMA.R
##After deconvolution from ENIGMA, we recommend retaining exclusively protein coding genes and removing lowly expressed genes to mitigate technical noise.
## The retained-gene filtering and inverse-normal transform are script M02_DataProcess_exp_fordb.R (driver step M02exp).

### Step 2 Genotype preprocessing
## Purpose： Perform quality control and filtering on genotype data (SNP array as example) matched to bulk RNA-seq samples, and compute genotype dosages

### Step 3 Covariate selection
## Script：M02_DataProcess_covariance_for_db.R (driver step M02cov; PEER factors + covariates)

### Step 4 TWAS model training
## Purpose：Train prediction models to generate input files for S-PrediXcan
## Script：M03_ModelTraining.R (helper: script/helpers/gtex_v7_nested_cv_elnet.R), then M04_Makedb.R (driver step M04db) builds the SQLite weight databases

### Step 5 LD score preparation from matched 1000 Genomes population as covariates for S-PrediXcan
## Script：M06_LD.sh (driver step M06ld, dosage extraction/alignment) and M06_Makeld.sh (driver step M06cov, LD covariance matrices)

### Step 6 spredixcan calculation
## Purpose：Meta-analyse GWAS summary statistics (M05_DataProcess_GWAS_align4METAL.R / M05_DataProcess_GWAS_metal.sh / M05_DataProcess_GWAS_metal_filter.R, driver steps M05align/M05metal/M05filter), then run S-PrediXcan per cell type
## Script：M07_TWAS.sh (driver step M07)

### Step 7 Downstream statistical analysis
## Script：M08_ResultStatistic.R (driver step M08)
