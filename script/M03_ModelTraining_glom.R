library(stringr)
library(dplyr)
library(parallel)
"%&%" <- function(a,b) paste(a,b, sep='')
argv <- commandArgs(trailingOnly = TRUE)
CELLTYPE <- argv[1]
CHROM <- argv[2]
PEER <- argv[3]
#snp_annot_file <- "/node1/liuxy/MediCell/kidney/glom/predictdb/usefile/split_chr/snp_annot/snp_annot.chr" %&% CHROM %&% ".txt"
snp_annot_file <- "/node1/liuxy/test_peer_glom/exp0.2/predictDB/usefile/new_gt/snp_annot/snp_annot.chr" %&% CHROM %&% ".txt"
gene_annot_file <- "/node1/liuxy/MediCell/kidney/glom/predictdb/usefile/gene_annot.parsed.txt"
genotype_file <- "/node1/liuxy/test_peer_glom/exp0.2/predictDB/usefile/new_gt/geno/genotype.chr" %&% CHROM %&% ".txt"
covariates_file <- "/node1/liuxy/MediCell/kidney_3/glom/predictdb/usefile/cov/" %&% CELLTYPE %&% "_peer15.txt"
expression_file <- "/node1/liuxy/chplot/kidney_new/output/glom/exp/" %&% CELLTYPE %&% "_transformed_expression.txt"

out_dir <- "/node1/liuxy/MediCell/kidney_3/glom/predictdb/output/"
dir.create(path = out_dir, recursive = TRUE, showWarnings = FALSE)
setwd("/node1/liuxy/MediCell/kidney_3/glom/predictdb/output/")
prmtr <- data.frame(chrom = CHROM, snp_annot_file = snp_annot_file, gene_annot_file = gene_annot_file, genotype_file = genotype_file, covariates_file = covariates_file, expression_file = expression_file)
prmtr$celltype <- str_split_i(str_split_i(prmtr$covariates_file, "/", -1), "_", 1)
prmtr$prefix <- str_split_i(str_split_i(prmtr$covariates_file, "/", -1), pattern = "_peer15.txt", 1)
prmtr$prefix <- paste0("Model_training_", prmtr$prefix)
source("/node1/liuxy/add_plot/pbmc_compare/script/gtex_v7_nested_cv_elnet.R")
mapply(main,prmtr$snp_annot_file,prmtr$gene_annot_file,prmtr$genotype_file,prmtr$expression_file,as.character(prmtr$covariates_file),as.numeric(prmtr$chrom),prmtr$prefix,MoreArgs = list(null_testing = FALSE))
