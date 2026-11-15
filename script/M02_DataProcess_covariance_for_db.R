# M02: PEER residual factors as covariates for model training (generic, all tissues)
# conda env: r3.6.3 (PEER)
#
# Covariates input (MC_COVARIATES): tab- or space-separated txt with header.
# Two layouts are auto-detected against the expression sample IDs:
#   - recommended: rows = samples, first column = sample ID, then covariate columns
#   - legacy: rows = covariates, columns = samples, first column = covariate names
# Numeric columns pass through in file order; columns listed in MC_CAT_COVS
# (comma-separated) are dummy-encoded, and dummy levels with fewer than
# MC_MIN_LEVEL_COUNT samples are dropped.
#
# For each cell type in MC_CELLTYPES, PEER estimates MC_PEER_NK hidden factors
# from the M02exp expression matrix while controlling for the known covariates.
# The first MC_PEER_KEEP factors are appended and written as {ct}_covariates.txt
# (rows = covariates, columns = samples, aligned to expression sample order).
library(peer)
library(parallel)
library(data.table)

script_file <- sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(script_file)), "helpers", "input_contracts.R"))
exp_dir <- Sys.getenv("MC_EXP_DIR")
cov_dir <- Sys.getenv("MC_COV_DIR")
peer_nk <- mc_number("MC_PEER_NK", "2", lower = 0, integer = TRUE)
peer_keep <- mc_number("MC_PEER_KEEP", "2", lower = 0, upper = peer_nk, integer = TRUE)
min_level_count <- mc_number("MC_MIN_LEVEL_COUNT", "5", lower = 1, integer = TRUE)
cat_covs <- strsplit(Sys.getenv("MC_CAT_COVS", ""), ",")[[1]]
cat_covs <- cat_covs[nchar(cat_covs) > 0]
celltypes <- mc_celltypes()
layout <- mc_choice("MC_COV_LAYOUT", c("auto", "samples", "covariates"))

covs_raw <- read.table(Sys.getenv("MC_COVARIATES"), header = TRUE,
                       stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
mc_ids(colnames(covs_raw), "Known covariate header")
samples_ref <- rownames(mc_matrix(file.path(exp_dir, paste0(celltypes[1], "_transformed_expression.txt")), "Expression"))
if (layout == "auto") {
  target <- mc_sample_key(samples_ref, "Expression samples")
  row_match <- all(target %in% mc_sample_key(covs_raw[[1]], "Known covariate first column"))
  col_match <- all(target %in% mc_sample_key(colnames(covs_raw)[-1], "Known covariate sample header"))
  if (row_match == col_match) stop("MC_COV_LAYOUT=auto requires exactly one complete layout match; specify samples | covariates")
  layout <- if (row_match) "samples" else "covariates"
}
if (layout == "covariates") {
  sample_ids <- colnames(covs_raw)[-1]
  covs <- as.data.frame(t(covs_raw[, -1, drop = FALSE]), stringsAsFactors = FALSE)
  colnames(covs) <- as.character(covs_raw[[1]])
} else {
  sample_ids <- as.character(covs_raw[[1]])
  covs <- covs_raw[, -1, drop = FALSE]
}
invisible(mc_sample_key(sample_ids, "Known covariate sample IDs"))
rownames(covs) <- sample_ids
mc_ids(colnames(covs), "Known covariate names")
if (!all(cat_covs %in% colnames(covs))) stop("MC_CAT_COVS must name existing covariates: ", paste(colnames(covs), collapse = ", "))

is_cat <- colnames(covs) %in% cat_covs
# Preserve the existing categorical model.matrix treatment and rare-level filter.
known <- matrix(numeric(0), nrow = nrow(covs), ncol = 0, dimnames = list(sample_ids, NULL))
if (any(!is_cat)) {
  num_covs <- as.data.frame(lapply(covs[, !is_cat, drop = FALSE], function(x) suppressWarnings(as.numeric(x))))
  colnames(num_covs) <- colnames(covs)[!is_cat]
  known <- as.matrix(num_covs)
  rownames(known) <- sample_ids
  mc_numeric(known, "Known numeric covariates")
}
if (any(is_cat)) {
  cat_raw <- covs[, is_cat, drop = FALSE]
  if (anyNA(cat_raw) || any(as.matrix(cat_raw) == "")) stop("Categorical covariates must not contain missing values")
  cat_raw[] <- lapply(cat_raw, factor)
  dummies <- model.matrix(~ ., data = cat_raw)[, -1, drop = FALSE]
  dummies <- dummies[, colSums(dummies) >= min_level_count, drop = FALSE]
  known <- cbind(known, dummies)
}

run_celltype <- function(ctname) {
  expr <- mc_matrix(file.path(exp_dir, paste0(ctname, "_transformed_expression.txt")), "Expression")
  storage.mode(expr) <- "double"
  samples_ordr <- rownames(expr)
  idx <- mc_sample_match(samples_ordr, sample_ids, "Known covariate sample IDs")
  thiscovs <- known[idx, , drop = FALSE]
  if (peer_nk > 0 && peer_nk + ncol(thiscovs) + 1L >= min(dim(expr)))
    stop("PEER requires MC_PEER_NK + known covariates + mean < both sample and gene counts; use MC_PEER_NK=0/MC_PEER_KEEP=0 for known covariates only")
  # Exact mode never removes a real X prefix. Legacy output retains original known IDs.
  output_ids <- if (mc_choice("MC_SAMPLE_ID_POLICY", c("exact", "legacy-rnames")) == "legacy-rnames") sample_ids[idx] else samples_ordr
  if (peer_nk == 0L) {
    out <- t(thiscovs) # explicit known-covariates-only mode; no PEER fit
  } else {
    model <- PEER()
    PEER_setPhenoMean(model, expr)
    PEER_setNk(model, peer_nk)
    PEER_setCovariates(model, thiscovs)
    PEER_setAdd_mean(model, TRUE)
    PEER_update(model)
    x <- PEER_getX(model)
    res <- list(factors = x, weights = PEER_getW(model), precision = PEER_getAlpha(model), residuals = PEER_getResiduals(model))
    save(res, file = file.path(cov_dir, paste0(ctname, "_peer_res.RData")))
    mean_col <- ncol(thiscovs) + 1L
    if (ncol(x) != mean_col + peer_nk || !all(x[, mean_col] == 1)) stop("Unexpected PEER factor layout")
    out <- t(x[, -mean_col, drop = FALSE])
    rownames(out) <- c(colnames(thiscovs), paste0("peer", seq_len(peer_nk)))
    keep_rows <- c(seq_len(ncol(thiscovs)), ncol(thiscovs) + seq_len(peer_keep))
    out <- out[keep_rows, , drop = FALSE]
  }
  colnames(out) <- output_ids
  mc_numeric(out, "Output covariates")
  file <- file.path(cov_dir, paste0(ctname, "_covariates.txt"))
  write.table(as.data.frame(out), file, sep = "\t", row.names = TRUE)
  file
}
results <- mclapply(celltypes, function(ct) try(run_celltype(ct), silent = TRUE), mc.cores = 6)
failed <- vapply(results, inherits, logical(1), "try-error")
if (any(failed)) stop("M02cov failed for ", paste(celltypes[failed], collapse = ", "), ": ",
                      paste(unlist(results[failed]), collapse = "\n"))
print(unlist(results))
