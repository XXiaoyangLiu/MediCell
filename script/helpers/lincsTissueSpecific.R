# Filter LINCS enrichment hits, exclude unwanted perturbation classes,
# rank by p-value/padj, and export top-N perturbations with their leading-edge genes.
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript lincsTissueSpecific.R <lincs_csv> <twas_label> <out_txt> <out_rdata> [top_n] [exclude_pattern] [rank_col]")
}

lincs_csv <- args[1]
twas_label <- args[2]
out_txt <- args[3]
out_rdata <- args[4]
top_n <- if (length(args) >= 5 && nzchar(args[5])) as.integer(args[5]) else 50L
exclude_arg <- if (length(args) >= 6) args[6] else "^BRD-|trt_sh\\.css"
exclude_pattern <- if (exclude_arg == "NONE") "" else exclude_arg
rank_col <- if (length(args) >= 7 && nzchar(args[7])) args[7] else "pval"

if (!file.exists(lincs_csv)) {
  stop("LINCS results CSV not found: ", lincs_csv)
}
if (is.na(top_n) || top_n < 1L) {
  stop("Invalid top_n value: ", args[5])
}

dat <- read.csv(lincs_csv, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("pert_iname", rank_col, "leadingEdge")
missing_cols <- setdiff(required_cols, colnames(dat))
if (length(missing_cols) > 0) {
  stop("Missing required columns in LINCS CSV: ", paste(missing_cols, collapse = ", "))
}
if (!"pert_type" %in% colnames(dat)) {
  dat$pert_type <- ""
}

if ("TWAS" %in% colnames(dat) && nzchar(twas_label)) {
  matched <- dat[dat$TWAS %in% twas_label, , drop = FALSE]
  if (nrow(matched) > 0) {
    dat <- matched
  } else {
    warning("No rows matched TWAS == ", twas_label, "; using all rows in ", lincs_csv)
  }
}

dat <- dat[!is.na(dat$pert_iname) & nzchar(trimws(dat$pert_iname)), , drop = FALSE]
dat <- dat[!is.na(dat$leadingEdge) & nzchar(trimws(dat$leadingEdge)), , drop = FALSE]

if (nzchar(exclude_pattern)) {
  combined_tag <- trimws(paste(dat$pert_iname, dat$pert_type))
  dat <- dat[!grepl(exclude_pattern, combined_tag), , drop = FALSE]
}

dat[[rank_col]] <- suppressWarnings(as.numeric(dat[[rank_col]]))
dat <- dat[is.finite(dat[[rank_col]]), , drop = FALSE]
dat <- dat[order(dat[[rank_col]]), , drop = FALSE]
dat <- dat[!duplicated(dat$pert_iname), , drop = FALSE]

n_keep <- min(top_n, nrow(dat))
top_df <- if (n_keep > 0) dat[seq_len(n_keep), , drop = FALSE] else dat[0, , drop = FALSE]
genes_df <- top_df[, c("pert_iname", "leadingEdge"), drop = FALSE]

dir.create(dirname(out_txt), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_rdata), recursive = TRUE, showWarnings = FALSE)

write.table(
  genes_df,
  file = out_txt,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
save(top_df, file = out_rdata)
message("Selected ", nrow(top_df), " perturbations (requested top ", top_n, ", ranked by ", rank_col, ") -> ", out_txt)
