# MediCell
1.Overview
MediCell is an integrative computational framework that combines single-cell transcriptomic reference atlases and   transcriptome-wide association study (TWAS) methodologies to dissect cell-type-specific genetic effects underlying complex human traits. Four categories of input data are mandatory for end-to-end pipeline execution:
Paired bulk RNA-seq and genomic genotype data from the same cohort
Healthy tissue single-cell RNA-seq reference atlas
Summary statistics of genome-wide association study (GWAS) for target complex traits
Pre-trained SPrediXcan/TWAS expression weight panels

2.Data Requirements
  2.1.Bulk transcriptomic profiling: BulkRNA-seq expression matrix from 大队列 study cohort individuals（bulk链接）
  2.2.Genotype data: Genetic variants matched with 2.1, either whole-genome sequencing (WGS) or SNP array.参考（genotype的链接）
  2.3.scRNA-seq reference atlas: scRNA-seq profiles derived from the same human tissue as bulk RNA-seq samples (healthy individuals only)Annotated cell type metadata is required (cell type labels assigned to each cell); serves to decompose bulk transcriptomic signals and map global genetic associations to distinct cell populations
  2.4.来自1000G的ld reference:下载自（1000G），匹配种族
  2.5.GWAS summary statistics: Standardized GWAS summary results for complex traits of interest (e.g., kidney function indices, autoimmune diseases)Must contain mandatory columns: variant ID, chromosome, position, effect allele, reference allele, regression beta, standard error, P-value; compatible with SPrediXcan input format
Precomputed TWAS expression weights.Trained cis-expression quantitative trait loci (cis-eQTL) weights compatible with SPrediXcan/MetaXcan, used to impute  genetically regulated gene expression (GReX) from individual genotypes

3.Software Requirements

4.Result
