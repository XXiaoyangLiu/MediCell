setwd('/node1/liuxy/MediCell/kidney_3/tub/predictdb/output')
library(stringr)
library(dplyr)
library(parallel)
"%&%" <- function(a,b) paste(a,b, sep='')
argv <- commandArgs(trailingOnly = TRUE)
CELLTYPE <- argv[1]
CHROM <- argv[2]

snp_annot_file <- "/node1/liuxy/MediCell/kidney/tub/predictdb_2/usefile/snp_annot/split_chr/snp_annot.chr" %&% CHROM %&% ".txt"
gene_annot_file <- "/node1/liuxy/MediCell/kidney/glom/predictdb/usefile/gene_annot.parsed.txt"
genotype_file <- "/node1/liuxy/MediCell/kidney/tub/predictdb_2/usefile/genotype/split_chr/genotype.chr" %&% CHROM %&% ".txt"
covariates_file <- "/node1/liuxy/MediCell/kidney_3/tub/predictdb/usefile/cov/" %&% CELLTYPE %&% "_peer15.txt"
expression_file <- "/node1/liuxy/chplot/kidney_new/output/tub/exp/" %&% CELLTYPE %&% "_transformed_expression.txt"

prmtr <- data.frame(chrom = CHROM, snp_annot_file = snp_annot_file, gene_annot_file = gene_annot_file, genotype_file = genotype_file, covariates_file = covariates_file, expression_file = expression_file)
prmtr$celltype <- str_split_i(str_split_i(prmtr$covariates_file, "/", -1), "_", 1)
prmtr$prefix <- str_split_i(str_split_i(prmtr$covariates_file, "/", -1), pattern = "_peer15.txt", 1)
prmtr$prefix <- paste0("Model_training_", prmtr$prefix)
source("/node1/liuxy/add_plot/pbmc_compare/script/gtex_v7_nested_cv_elnet.R")
mapply(main,prmtr$snp_annot_file,prmtr$gene_annot_file,prmtr$genotype_file,prmtr$expression_file,as.character(prmtr$covariates_file),as.numeric(prmtr$chrom),prmtr$prefix,MoreArgs = list(null_testing = FALSE))
