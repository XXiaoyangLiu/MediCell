# L02: Run LINCS fgsea enrichment on directional TWAS genes and annotate hits.
options(stringsAsFactors = FALSE)

req_env <- function(key) {
  val <- Sys.getenv(key)
  if (!nzchar(val)) stop("Required environment variable is missing: ", key)
  val
}

opt_env <- function(key, default = "") {
  val <- Sys.getenv(key)
  if (!nzchar(val)) default else val
}

pick_col <- function(cols, explicit, candidates) {
  if (nzchar(explicit) && explicit %in% cols) return(explicit)
  hit <- candidates[candidates %in% cols]
  if (length(hit) > 0) hit[1] else ""
}

script_file <- sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(script_file)), "helpers", "lincs_fgsea_lib.R"))

l02_dir <- req_env("MC_LINCS_L02_DIR")
twas_dir_csv <- req_env("MC_LINCS_DIR_CSV")
trait <- req_env("MC_LINCS_TRAIT")
label <- req_env("MC_LINCS_LABEL")
direction <- req_env("MC_LINCS_DIRECTION")
gmt_file <- req_env("MC_LINCS_GMT")
pert1_file <- opt_env("MC_LINCS_PERT_PHASE1")
pert2_file <- opt_env("MC_LINCS_PERT_PHASE2")

gene_col_in <- opt_env("MC_LINCS_GENE_COL")
z_col_in <- opt_env("MC_LINCS_Z_COL")
padj_col_in <- opt_env("MC_LINCS_PADJ_COL")

zscore_cutoff <- as.numeric(opt_env("MC_LINCS_ZSCORE_CUTOFF", "0.25"))
padj_cutoff <- as.numeric(opt_env("MC_LINCS_TWAS_PADJ_CUTOFF", "0.05"))
min_size <- as.integer(opt_env("MC_LINCS_MIN_SIZE", "10"))
max_size <- as.integer(opt_env("MC_LINCS_MAX_SIZE", "500"))
cores <- as.integer(opt_env("MC_LINCS_CORES", "4"))
sig_col <- opt_env("MC_LINCS_SIG_COL", "pval")
sig_cutoff <- as.numeric(opt_env("MC_LINCS_SIG_CUTOFF", "0.05"))
nes_filter <- opt_env("MC_LINCS_NES_FILTER", "neg")
remove_ribo <- opt_env("MC_LINCS_REMOVE_RIBO", "0") == "1"

if (!file.exists(twas_dir_csv)) stop("Directional TWAS CSV not found: ", twas_dir_csv)
dir.create(l02_dir, recursive = TRUE, showWarnings = FALSE)

deg_table <- read.csv(twas_dir_csv, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
cols <- colnames(deg_table)
gene_col <- pick_col(cols, gene_col_in, c("Symbol", "gene_name", "Gene", "symbol"))
z_col <- pick_col(cols, z_col_in, c("Zscore", "zscore", "Z"))
padj_col <- pick_col(cols, padj_col_in, c("FDR_P", "fdr_p", "FDR P.value", "padj", "p_val_adj"))
if (!nzchar(gene_col) || !nzchar(z_col)) {
  stop("Missing gene or z-score column in directional TWAS CSV: ", twas_dir_csv)
}
padj_arg <- if (nzchar(padj_col)) padj_col else NULL

fgsea_dir <- if (direction == "pos") "up" else if (direction == "neg") "down" else infer_direction_by_zscore(deg_table, z_col)

message("Loading GMT: ", gmt_file)
pathways <- parse_lincs_gmt(gmt_file)
gmt_desc <- attr(pathways, "gmt_desc")
message("Running fgsea (", fgsea_dir, ", minSize=", min_size, ", maxSize=", max_size, ", cores=", cores, ") ...")

res <- run_lincs_fgsea(
  deg_table = deg_table,
  direction = fgsea_dir,
  pathways = pathways,
  gene_symbol_col = gene_col,
  zscore_col = z_col,
  padj_col = padj_arg,
  zscore_cutoff = zscore_cutoff,
  padj_cutoff = padj_cutoff,
  min_size = min_size,
  max_size = max_size,
  cores = cores,
  remove_ribo = remove_ribo,
  sig_col = sig_col,
  sig_cutoff = sig_cutoff
)

pert_anno <- load_pert_annotation(c(pert1_file, pert2_file))
twas_tag <- paste0(trait, "_", label, "_", direction)

all_anno <- annotate_lincs_fgsea(res$fgsea_all, pert_anno = pert_anno, gmt_desc = gmt_desc)
sig_anno <- annotate_lincs_fgsea(res$fgsea_sig, pert_anno = pert_anno, gmt_desc = gmt_desc)

all_anno$TWAS <- if (nrow(all_anno) > 0) rep(twas_tag, nrow(all_anno)) else character(0)
sig_anno$TWAS <- if (nrow(sig_anno) > 0) rep(twas_tag, nrow(sig_anno)) else character(0)

if (nes_filter == "neg" && nrow(sig_anno) > 0) {
  sig_anno <- sig_anno[!is.na(sig_anno$NES) & sig_anno$NES < 0, , drop = FALSE]
} else if (nes_filter == "pos" && nrow(sig_anno) > 0) {
  sig_anno <- sig_anno[!is.na(sig_anno$NES) & sig_anno$NES > 0, , drop = FALSE]
}

all_out <- file.path(l02_dir, paste0("LINCS_all_", twas_tag, ".csv"))
sig_out <- file.path(l02_dir, paste0("LINCS_sig_", twas_tag, ".csv"))
write.csv(all_anno, file = all_out, row.names = FALSE)
write.csv(sig_anno, file = sig_out, row.names = FALSE)

message("LINCS enrichment complete: ", nrow(all_anno), " tested, ",
        nrow(sig_anno), " significant (", sig_col, "<", sig_cutoff,
        ", NES=", nes_filter, ") -> ", sig_out)
