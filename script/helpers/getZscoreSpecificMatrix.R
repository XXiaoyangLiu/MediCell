# Compute drug x cell-type prioritization scores by combining |TWAS Z-score|
# with cross-cell-type expression specificity (expZscore) across leading-edge genes.
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript getZscoreSpecificMatrix.R <zscore_txt> <exp_csv> <drug_pairs_txt> <out_rdata> <out_csv>")
}

zscore_file <- args[1]
exp_file <- args[2]
drugpair_file <- args[3]
out_rdata <- args[4]
out_csv <- args[5]

for (f in c(zscore_file, exp_file, drugpair_file)) {
  if (!file.exists(f)) stop("Required input file not found: ", f)
}

library(tidyr)

df <- read.table(zscore_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
exp <- read.csv(exp_file, header = TRUE, na.strings = c("", "NA"), row.names = 1, check.names = FALSE)
if (ncol(exp) < 2L) {
  stop("Cross-cell-type LINCS scoring requires at least 2 cell types in expression matrix; found ",
       ncol(exp), " (", paste(colnames(exp), collapse = ", "), ").")
}
exp[is.na(exp)] <- 0

dir.create(dirname(out_rdata), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

if (nrow(df) == 0 || !file.info(drugpair_file)$size) {
  result <- matrix(numeric(0), nrow = 0, ncol = ncol(exp), dimnames = list(character(0), colnames(exp)))
  save(result, file = out_rdata)
  write.csv(data.frame(pert_iname = character(0), check.names = FALSE), file = out_csv, row.names = FALSE)
  message("No overlapping TWAS z-scores or drug-gene pairs; wrote empty score matrix.")
  quit(save = "no", status = 0)
}

mat <- tidyr::pivot_wider(
  df,
  names_from = CellType,
  values_from = Zscore,
  values_fill = 0
)
twas_df <- as.data.frame(mat, stringsAsFactors = FALSE)
rownames(twas_df) <- twas_df$Gene
twas_mat <- as.matrix(twas_df[, -1, drop = FALSE])
storage.mode(twas_mat) <- "numeric"
twas_mat <- abs(twas_mat)

common_genes <- intersect(rownames(twas_mat), rownames(exp))
common_cts <- intersect(colnames(twas_mat), colnames(exp))
if (length(common_genes) == 0 || length(common_cts) == 0) {
  stop("No overlap between TWAS genes/cell types and expression matrix (",
       length(common_genes), " genes, ", length(common_cts), " cell types).")
}
if (length(common_cts) < 2L) {
  stop("Cross-cell-type LINCS scoring requires at least 2 overlapping cell types between TWAS and expression matrix; found ",
       length(common_cts), " (", paste(common_cts, collapse = ", "), ").")
}

twas_sub <- twas_mat[common_genes, common_cts, drop = FALSE]
subexp <- as.matrix(exp[common_genes, common_cts, drop = FALSE])
storage.mode(subexp) <- "numeric"

exp_mean <- rowMeans(subexp, na.rm = TRUE)
exp_sd <- apply(subexp, 1, sd, na.rm = TRUE)
exp_zscore <- (subexp - exp_mean) / exp_sd
zero_sd <- is.na(exp_sd) | exp_sd == 0 | !is.finite(exp_sd)
if (any(zero_sd)) {
  exp_zscore[zero_sd, ] <- 0
}
exp_zscore[!is.finite(exp_zscore)] <- 0

rank_score <- twas_sub * exp_zscore

lededge <- read.table(drugpair_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "")
drugs <- unique(lededge[, 1])
rows <- lapply(drugs, function(drug) {
  genes <- unique(lededge[lededge[, 1] == drug, 2])
  genes <- intersect(genes, rownames(rank_score))
  if (length(genes) == 0) {
    return(setNames(numeric(ncol(rank_score)), colnames(rank_score)))
  }
  colSums(rank_score[genes, , drop = FALSE], na.rm = TRUE)
})
result <- do.call(rbind, rows)
rownames(result) <- drugs
colnames(result) <- colnames(rank_score)

save(result, file = out_rdata)
out_df <- data.frame(pert_iname = rownames(result), as.data.frame(result, check.names = FALSE),
                     row.names = NULL, check.names = FALSE)
write.csv(out_df, file = out_csv, row.names = FALSE)
message("Wrote drug x cell-type matrix (", nrow(result), " drugs x ", ncol(result), " cell types) -> ", out_csv)
