# M08: selected CTs, named columns, per-CT corrections (historical output layout).
required <- function(key) {
  value <- Sys.getenv(key)
  if (!nzchar(value)) stop("Required environment: ", key)
  value
}
TWAS_dir <- required("MC_TWAS_DIR")
gwas_id <- required("MC_GWAS_ID")
output_dir <- required("MC_SUMMARY_OUT")
celltypes <- strsplit(required("MC_CELLTYPES"), "[[:space:]]+")[[1]]
if (anyDuplicated(celltypes) || any(!grepl("^[A-Za-z0-9_.-]+$", celltypes))) stop("Invalid CT selection")
genes <- read.csv(required("MC_GENE_POS"), stringsAsFactors=FALSE, check.names=FALSE)
if (!all(c("gene_id", "gene_name") %in% names(genes)) || anyNA(genes$gene_id) ||
    any(!nzchar(genes$gene_id)) || anyDuplicated(genes$gene_id)) stop("Invalid/duplicate gene_pos gene_id mapping")
columns <- c("gene", "gene_name", "zscore", "effect_size", "pvalue", "var_g", "pred_perf_r2",
             "pred_perf_pval", "pred_perf_qval", "n_snps_used", "n_snps_in_cov", "n_snps_in_model")
cores <- suppressWarnings(as.integer(Sys.getenv("MC_CORES", "4")))
if (is.na(cores) || cores < 1L) cores <- 1L
df_list <- parallel::mclapply(celltypes, function(ct) try({
  file <- file.path(TWAS_dir, paste0(ct, ".csv"))
  dat <- read.csv(file, stringsAsFactors=FALSE, check.names=FALSE)
  if (!nrow(dat) || anyDuplicated(names(dat)) || !setequal(columns, names(dat))) stop("Empty/unknown TWAS schema: ", file)
  dat <- dat[columns]
  if (!is.numeric(dat$pvalue) || any(!is.finite(dat$pvalue)) || any(dat$pvalue < 0 | dat$pvalue > 1) ||
      anyNA(dat$gene) || any(!nzchar(dat$gene)) || anyDuplicated(dat$gene)) stop("Invalid TWAS gene/pvalue: ", file)
  mapping <- match(dat$gene, genes$gene_id)
  if (anyNA(mapping) || anyNA(genes$gene_name[mapping]) || any(!nzchar(genes$gene_name[mapping]))) stop("Unmatched/missing gene_pos mapping: ", file)
  dat$fdr_p <- p.adjust(dat$pvalue, method="BH")
  dat$cell_type <- ct
  dat$ngeneindb <- nrow(dat)
  dat$bonferroni_p <- p.adjust(dat$pvalue, method="bonferroni")
  dat$gwas_id <- gwas_id
  dat$gene_name <- genes$gene_name[mapping]
  names(dat)[names(dat)=="gene"] <- "gene_id"
  dat
}, silent = TRUE), mc.cores = min(cores, length(celltypes)))
failed <- vapply(df_list, inherits, logical(1), "try-error")
if (any(failed)) stop("M08 failed for ", paste(celltypes[failed], collapse = ", "), ": ",
                      paste(unlist(df_list[failed]), collapse = "
"))
final_df <- do.call(rbind, df_list)
bonf_sig_df <- subset(final_df, bonferroni_p < 0.05)
fdr_sig_df <- subset(final_df, fdr_p < 0.05)
all_sig_df <- unique(rbind(bonf_sig_df, fdr_sig_df))
final_df$FDR.sig <- ifelse(final_df$fdr_p < 0.05, "TRUE", "FALSE")
final_df$bonferroni_p.sig <- ifelse(final_df$bonferroni_p < 0.05, "TRUE", "FALSE")
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
write.csv(final_df, file.path(output_dir, "allcelltype_result.csv"), row.names=FALSE)
write.csv(bonf_sig_df, file.path(output_dir, "Bonfsig_allcelltype_result.csv"), row.names=FALSE)
write.csv(fdr_sig_df, file.path(output_dir, "FDRsig_allcelltype_result.csv"), row.names=FALSE)
write.csv(all_sig_df, file.path(output_dir, "BonFDRsig_allcelltype_result.csv"), row.names=FALSE)
cat("Total:", nrow(final_df), "FDR:", nrow(fdr_sig_df), "Bonferroni:", nrow(bonf_sig_df), "union:", nrow(all_sig_df), "\n")
