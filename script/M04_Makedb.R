# M04db: build downstream databases from the fixed 24-column M03 output.
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(dplyr))
required <- function(key) {
  x <- Sys.getenv(key)
  if (!nzchar(x)) stop("Required environment: ", key)
  x
}
summary_dir <- required("MC_M03_SUMMARY_DIR")
weights_dir <- required("MC_M03_WEIGHTS_DIR")
db_dir <- required("MC_DB_DIR")
celltypes <- strsplit(required("MC_CELLTYPES"), "[[:space:]]+")[[1]]
chr_text <- strsplit(required("MC_CHROMS"), ",", fixed=TRUE)[[1]]
if (any(!grepl("^[0-9]+$", chr_text))) stop("MC_CHROMS must be CSV integers")
chr_list <- sort(as.integer(chr_text))
if (anyDuplicated(chr_list) || any(!chr_list %in% 1:22) || anyDuplicated(celltypes) ||
    any(!grepl("^[A-Za-z0-9_.-]+$", celltypes))) stop("Invalid chromosome/celltype selection")
model_set <- required("MC_TWAS_MODEL_SET")
if (!model_set %in% c("unfiltered", "filtered-signif")) stop("Invalid MC_TWAS_MODEL_SET")
population <- Sys.getenv("MC_TRAIN_POPULATION", "unspecified")
if (!nzchar(population)) stop("Declare MC_TRAIN_POPULATION or unspecified")
build <- required("MC_GENOME_BUILD")
expected <- strsplit("gene_name gene_id gene_type alpha n_snps_in_window n_snps_in_model lambda_min_mse test_R2_avg test_R2_sd cv_R2_avg cv_R2_sd in_sample_R2 nested_cv_fisher_pval rho_avg rho_se rho_zscore rho_avg_squared zscore_pval cv_rho_avg cv_rho_se cv_rho_avg_squared cv_zscore_est cv_zscore_pval cv_pval_est", " ")[[1]]
read_summary <- function(path) {
  widths <- count.fields(path, sep="\t", quote="", comment.char="")
  if (length(widths) < 2 || any(widths != 24L))
    stop("Expected fixed 24-column M03 summary; regenerate old training output: ", path)
  dat <- read.delim(path, check.names=FALSE, quote="", comment.char="", fill=FALSE,
                    colClasses=c(rep("character", 3), rep("numeric", 21)))
  if (!identical(names(dat), expected)) stop("Unknown summary header: ", path)
  counts <- as.matrix(dat[c("n_snps_in_window", "n_snps_in_model")])
  if (anyNA(dat[1:3]) || any(!nzchar(as.matrix(dat[1:3]))) ||
      any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", dat$gene_type)) ||
      anyNA(dat$alpha) || any(dat$alpha != 0.5) ||
      any(!is.finite(counts)) || any(counts < 0 | counts %% 1 != 0))
    stop("Invalid model identifiers, alpha or SNP counts: ", path)
  dat
}
indexes <- function(conn, table) {
  dbExecute(conn, "CREATE INDEX weights_rsid ON weights (rsid)")
  dbExecute(conn, "CREATE INDEX weights_gene ON weights (gene)")
  dbExecute(conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
  dbExecute(conn, paste0("CREATE INDEX gene_model_summary ON ", table, " (gene)"))
}
dir.create(db_dir, recursive=TRUE, showWarnings=FALSE)
for (ct in celltypes) {
  filename <- function(dir, chr, suffix) file.path(dir, paste0("Model_training_", ct, "_chr", chr, "_", suffix, ".txt"))
  parts <- lapply(chr_list, function(chr) read_summary(filename(summary_dir, chr, "model_summaries")))
  model_summaries <- bind_rows(parts) %>% rename(gene=gene_name)
  if (anyDuplicated(model_summaries$gene) || anyDuplicated(model_summaries$gene_id)) stop("Duplicate model gene mapping")
  expression_keys <- model_summaries$gene
  stable_keys <- model_summaries$gene_id
  model_summaries$gene <- stable_keys
  tiss <- lapply(chr_list, function(chr) {
    x <- read.table(filename(summary_dir, chr, "summary"), header=TRUE)
    if (!all(c("n_samples", "chrom", "cv_seed") %in% names(x)) || nrow(x)!=1 ||
        is.na(x$chrom) || x$chrom != chr) stop("Invalid tissue summary")
    x
  })
  tiss_summary <- bind_rows(tiss)
  n_samples <- unique(tiss_summary$n_samples)
  if (length(n_samples)!=1 || !is.finite(n_samples) || n_samples < 2 || n_samples %% 1 != 0) stop("Inconsistent/invalid n_samples")
  weights <- bind_rows(lapply(chr_list, function(chr) {
    w <- read.table(filename(weights_dir, chr, "weights"), header=TRUE, stringsAsFactors=FALSE,
                    colClasses=c(rep("character",5),"numeric"))
    if (!identical(names(w), c("gene_name", "rsid", "varID", "ref", "alt", "beta"))) stop("Unknown weights schema")
    pieces <- strsplit(w$varID, "_", fixed=TRUE)
    if (any(lengths(pieces)!=4)) stop("Invalid canonical VarID")
    for (i in seq_len(nrow(w))) {
      a <- pieces[[i]]
      if (a[1] != as.character(chr) || !grepl("^[1-9][0-9]*$", a[2]) ||
          !a[3] %in% c("A","C","G","T") || !a[4] %in% c("A","C","G","T") ||
          a[3]==a[4] || a[3]!=w$ref[i] || a[4]!=w$alt[i]) stop("VarID/allele/chromosome mismatch")
    }
    if (anyNA(w) || any(!nzchar(w$rsid)) || !is.numeric(w$beta) || any(!is.finite(w$beta))) stop("Invalid weights")
    w
  })) %>% rename(gene=gene_name)
  if (!nrow(weights)) stop("No nonzero weights across selected chromosomes: ", ct)
  if (any(!weights$gene %in% expression_keys)) stop("Unmatched weight genes")
  weights$gene <- stable_keys[match(weights$gene, expression_keys)]
  mapping <- unique(weights[c("rsid", "varID", "ref", "alt")])
  if (anyDuplicated(mapping$rsid) || anyDuplicated(mapping$varID) || anyDuplicated(weights[c("gene", "varID")])) stop("Variant mapping collision/duplicate weights")
  weights$rsid <- weights$varID
  significant <- model_summaries %>% filter(zscore_pval < 0.05, rho_avg_squared > 0.01)
  model_sets <- list(unfiltered=model_summaries)
  if (model_set == "filtered-signif") model_sets$filtered_signif <- significant
  targets <- file.path(db_dir, paste0(ct, "_models_", names(model_sets), ".db"))
  if (any(file.exists(targets))) stop("Refusing existing database output: ", ct)
  construction <- tiss_summary %>% select(chrom, cv_seed) %>% rename(chromosome=chrom)
  sample_info <- data.frame(n_samples=n_samples, population=population)
  for (i in seq_along(model_sets)) {
    models <- model_sets[[i]]
    extra <- models %>% rename(pred.perf.R2=rho_avg_squared, genename=gene_id, pred.perf.pval=zscore_pval, n.snps.in.model=n_snps_in_model)
    extra$genename <- expression_keys[match(models$gene, stable_keys)]
    extra$pred.perf.qval <- rep(NA, nrow(extra))
    w <- weights %>% filter(gene %in% models$gene) %>% rename(eff_allele=alt, ref_allele=ref, weight=beta)
    conn <- dbConnect(SQLite(), targets[i])
    dbWriteTable(conn, "extra", extra)
    dbWriteTable(conn, "construction", construction)
    dbWriteTable(conn, "sample_info", sample_info)
    dbWriteTable(conn, "weights", w)
    indexes(conn, "extra")
    dbDisconnect(conn)
  }
  write.table(data.frame(gene_name=expression_keys, gene_id=stable_keys),
              file.path(db_dir,paste0(ct,"_gene_sources.tsv")), sep="	", quote=FALSE, row.names=FALSE)
  write.table(mapping, file.path(db_dir, paste0(ct, "_variant_sources.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(data.frame(key=c("genome_build", "train_population", "summary_format", "chromosomes"),
                        value=c(build, population, "standard", paste(chr_list, collapse=","))),
              file.path(db_dir, paste0(ct, "_metadata.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  cat(ct, ":", nrow(model_summaries), "models;", nrow(significant), "significant;", nrow(weights), "weights\n")
}
