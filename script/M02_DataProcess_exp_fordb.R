## R = 3.6.3
# Expression matrices for prediction DB: filter low-expression genes, INT-transform,
# output per-cell-type matrices for predictDB / S-PrediXcan.
library(stringr)
library(parallel)
library(dplyr)
script_file <- sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(script_file)), "helpers", "input_contracts.R"))

load(Sys.getenv("MC_ENIGMA_RDATA"))
exp_dir <- Sys.getenv("MC_EXP_DIR")

# Gene filter thresholds: keep genes with inferred TPM > MC_TPM_MIN in at least
# max(ceil(MC_MIN_PCT * n_samples), MC_MIN_SAMPLES) samples.
# Defaults preserve the published filtering (0.1 / 1% / 5).
tpm_min <- mc_number("MC_TPM_MIN", "0.1", lower = 0)
min_pct <- mc_number("MC_MIN_PCT", "0.01", lower = 0, upper = 1)
min_samples_floor <- mc_number("MC_MIN_SAMPLES", "5", lower = 0, integer = TRUE)

if (!exists("enigma.cse.test") || length(dim(enigma.cse.test)) != 3L || any(dim(enigma.cse.test) == 0))
  stop("MC_ENIGMA_RDATA must contain a nonempty genes x samples x cell-types array enigma.cse.test")

for (i in 1:3) mc_ids(dimnames(enigma.cse.test)[[i]], paste("ENIGMA dimension", i))
cell_types <- mc_celltypes(dimnames(enigma.cse.test)[[3]])
filtered_genes <- list()
for (cell_type in cell_types) {
  expr_matrix <- matrix(enigma.cse.test[, , cell_type, drop = FALSE],
                        nrow = dim(enigma.cse.test)[1], ncol = dim(enigma.cse.test)[2],
                        dimnames = dimnames(enigma.cse.test)[1:2])
  mc_numeric(expr_matrix, paste("Expression for", cell_type), TRUE)
  # number of samples in which each gene is expressed above the TPM threshold
  n_samples_above_threshold <- rowSums(expr_matrix > tpm_min)
  # minimum sample count (min_pct fraction), and at least min_samples_floor samples
  min_samples_1 <- ceiling(min_pct * ncol(expr_matrix))
  min_samples <- max(min_samples_1, min_samples_floor)
  keep_genes <- n_samples_above_threshold >= min_samples

  if (!any(keep_genes)) stop("No genes pass expression thresholds for ", cell_type)
  filtered_genes[[cell_type]] <- rownames(expr_matrix)[keep_genes]
  message(sprintf("Cell type: %s, Genes kept: %d/%d (%.1f%%) [TPM > %g in >= %d samples]",
                  cell_type,
                  sum(keep_genes),
                  nrow(expr_matrix),
                  100 * sum(keep_genes)/nrow(expr_matrix),
                  tpm_min,
                  min_samples))
}

sapply(filtered_genes, length)

# keep only the filtered genes in enigma.cse.test
filtered_matrices <- list()
for (cell_type in cell_types) {
  all_genes <- dimnames(enigma.cse.test)[[1]]
  keep_genes <- filtered_genes[[cell_type]]
  keep_indices <- which(all_genes %in% keep_genes)
  filtered_mat <- matrix(enigma.cse.test[keep_indices, , cell_type, drop = FALSE],
                         nrow = length(keep_indices), ncol = dim(enigma.cse.test)[2],
                         dimnames = list(all_genes[keep_indices], dimnames(enigma.cse.test)[[2]]))
  filtered_matrices[[cell_type]] <- filtered_mat
}

# INT-transform the filtered matrices and write one samples x genes TXT per cell type.
cell_types <- names(filtered_matrices)
for (ctname in cell_types) {
  ct_matrix <- filtered_matrices[[ctname]]
  # transpose (samples x genes)
  tmp <- as.data.frame(t(ct_matrix))
  # INT transform
  tmp.INT <- matrix(0, nrow(tmp), ncol(tmp))
  for (i in seq_len(ncol(tmp))) {
    tmp.INT[, i] <- qqnorm(tmp[, i], plot.it = FALSE)$x
  }
  rownames(tmp.INT) <- rownames(tmp)
  colnames(tmp.INT) <- colnames(tmp)
  tmp.INT <- as.data.frame(tmp.INT)

  write.table(
    tmp.INT,
    file = file.path(exp_dir, paste0(ctname, "_transformed_expression.txt")),
    sep = "\t",
    row.names = TRUE
  )
}
