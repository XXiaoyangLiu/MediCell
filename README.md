# MediCell
1.Overview
MediCell is an integrative computational framework that combines single-cell transcriptomic reference atlases and transcriptome-wide association study (TWAS) methodologies to dissect cell-type-specific genetic effects underlying complex human traits. Four categories of input data are mandatory for end-to-end pipeline execution:
Paired bulk RNA-seq and genomic genotype data from the same cohort
Healthy tissue single-cell RNA-seq reference atlas
Summary statistics of genome-wide association study (GWAS) for target complex traits
Pre-trained SPrediXcan/TWAS expression weight panels

2. Data Requirements
MediCell requires 6 categories of input datasets with 3 cohort stratification:
  Discovery cohort (healthy donors): Paired Bulk RNA-seq transcriptomes, individual genotypes and phenotype covariates used for genetically regulated expression imputation and cell-type-specific TWAS calculation
  scRNA-seq reference cohort (healthy donors): Single-cell atlas in matched tissue exclusively serving as deconvolution reference
  GWAS cohort: Independent case/control or population cohort with trait GWAS summary statistics for interested phenotypes
Detailed specifications for each dataset are listed below:
2.1 Bulk transcriptomic profiles
  Gene expression matrix generated from bulk RNA-seq of the healthy discovery cohort.
  TPM-normalized expression values are recommended as the optimal input for cellular fraction deconvolution across this strategy.
2.2 Individual level genotype data
  Genetic variant profiles perfectly matched to samples in the healthy discovery cohort (identical individuals).
  Data can be derived from either whole-genome sequencing (WGS) or genome-wide SNP arrays; genotype dosage values are required.
2.3 Tissue matched scRNA-seq reference atlas
  scRNA-seq profiles acquired from the same tissue origin as bulk RNA-seq, obtained exclusively from healthy donors (independent reference cohort, no overlap with discovery cohort).
2.4 Covariate matrix for model training
  Covariates derived entirely from the healthy discovery cohort, consisting of three components:
  Individual covariates: Individual phenotypic information including age, sex and other required confounders.
  Genetic principal components (PCs): Ancestry PCs calculated from cohort genotypes.
  PEER factors: Residual confounding factors estimated via the PEER algorithm. PEER adjustment is performed on expression profiles deconvoluted by ENIGMA cellular and filtered low expressed gene for each cell type.
  All covariates are jointly incorporated to adjust for technical and biological confounding during model training.
2.5 LD reference panel
  Population matched linkage disequilibrium reference panel retrieved from the 1000 Genomes Project Phase 3(ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/).
2.6 GWAS summary statistics
  Genome-wide association study summary results for interested traits, generated from an independent GWAS cohort (distinct from the healthy discovery cohort).

3.Software & Computational Tool Requirements
3.1.shell & Python
MetaXcan(https://github.com/hakyimlab/MetaXcan)

3.2.R
ENIGMA and other required R packages on r4.3.1
PEER（need a independent enviroments on r3.6.3）

3.3.Network Tools
CIBERSORTx



