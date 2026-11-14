library(stringr)
library(dplyr)
library(parallel)
"%&%" <- function(a,b) paste(a,b, sep='')
argv <- commandArgs(trailingOnly = TRUE)
CELLTYPE <- argv[1]
CHROM <- argv[2]

sub_chr <- function(pattern) sub("{chr}", CHROM, pattern, fixed = TRUE)
snp_annot_file <- sub_chr(Sys.getenv("MC_SNP_ANNOT_PATTERN"))
gene_annot_file <- Sys.getenv("MC_GENE_ANNOT")
genotype_file <- sub_chr(Sys.getenv("MC_GENOTYPE_PATTERN"))
covariates_file <- Sys.getenv("MC_COV_FILE")
expression_file <- Sys.getenv("MC_EXP_FILE")

out_dir <- Sys.getenv("MC_M03_OUT")
dir.create(path = out_dir, recursive = TRUE, showWarnings = FALSE)
# the elnet core writes to relative subdirectories of the working directory
dir.create(file.path(out_dir, "summary"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "weights"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "covariances"), recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)
prmtr <- data.frame(chrom = CHROM, snp_annot_file = snp_annot_file, gene_annot_file = gene_annot_file, genotype_file = genotype_file, covariates_file = covariates_file, expression_file = expression_file)
prmtr$celltype <- CELLTYPE
prmtr$prefix <- paste0("Model_training_", CELLTYPE)
source(Sys.getenv("MC_ELNET_SCRIPT"))
seed_raw <- Sys.getenv("MC_SEED", "")
seed <- if (!nzchar(seed_raw)) NA else mc_number("MC_SEED", "", lower = 0, integer = TRUE)
mapply(main,prmtr$snp_annot_file,prmtr$gene_annot_file,prmtr$genotype_file,prmtr$expression_file,as.character(prmtr$covariates_file),as.numeric(prmtr$chrom),prmtr$prefix,MoreArgs = list(null_testing = FALSE, seed = seed))
