### Analysis Workflow Pipeline
## The entire analysis pipeline consists of 7 steps.

### Step 1 Cell proportion deconvolution estimation
## Purpose: Process scRNA-seq data to generate reference matrices for Cibersortx and ENIGMA estimation, and infer cell-type expression profiles from bulk RNA-seq data
## Script：
    pre_sc_RNA-seq.R
    Cibersort.R
    ENIGMA.R
    ENIGMA_qc.R
## Input：
## Output：

### Step 2 Genotype preprocessing
## Purpose： Perform quality control and filtering on genotype data (SNP array as example) matched to bulk RNA-seq samples, and compute genotype dosages
## Script：
    genotype.R
## Input:
## Output：

### Step 3 Covariate selection
## Purpose：Compute PEER factors to correct for latent confounding variables
## Script：
    peer.R

## Input：
## Output：

### Step 4 TWAS model training
## Purpose：Train prediction models to generate input files for S-PrediXcan
## Script：
    
## Input：
## Output：

### Step 5 LD score preparation from matched 1000 Genomes population as covariates for S-PrediXcan
## Purpose：
## Script：

## Input：
## Output：

### Step 6 spredixcan calculation
## Purpose：
## Script：

## Input：
## Output：

### Step 7 Downstream statistical analysis
## Purpose：
## Script：

## Input：
## Output：
