# Helper library for LINCS fgsea enrichment and perturbation annotation.

parse_lincs_gmt <- function(gmt_file) {
  if (!file.exists(gmt_file)) {
    stop("GMT file not found: ", gmt_file)
  }
  lines <- readLines(gmt_file, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) {
    stop("GMT file is empty: ", gmt_file)
  }
  parts <- strsplit(lines, "\t", fixed = TRUE)
  valid <- vapply(parts, function(x) length(x) >= 3L && nzchar(x[1]), logical(1))
  parts <- parts[valid]
  if (length(parts) == 0) {
    stop("No valid pathways parsed from GMT file: ", gmt_file)
  }
  pathways <- lapply(parts, function(x) {
    genes <- x[-(1:2)]
    genes <- genes[nzchar(genes)]
    unique(genes)
  })
  ids <- vapply(parts, function(x) x[1], character(1))
  descs <- vapply(parts, function(x) x[2], character(1))
  names(pathways) <- ids
  attr(pathways, "gmt_desc") <- setNames(descs, ids)
  pathways
}

infer_direction_by_zscore <- function(df, z_col = "zscore") {
  if (!z_col %in% colnames(df)) {
    stop("Missing z-score column: ", z_col)
  }
  z <- suppressWarnings(as.numeric(df[[z_col]]))
  z <- z[is.finite(z) & z != 0]
  if (length(z) == 0) {
    stop("No valid non-zero z-score values found.")
  }
  if (all(z > 0)) return("up")
  if (all(z < 0)) return("down")
  stop("Input contains both positive and negative z-scores; split by direction first.")
}

run_lincs_fgsea <- function(deg_table,
                            direction = c("up", "down"),
                            pathways,
                            gene_symbol_col = "gene_name",
                            zscore_col = "zscore",
                            padj_col = NULL,
                            zscore_cutoff = 0.25,
                            padj_cutoff = 0.05,
                            min_size = 10,
                            max_size = 500,
                            cores = 4,
                            remove_ribo = FALSE,
                            sig_col = c("pval", "padj"),
                            sig_cutoff = 0.05) {
  direction <- match.arg(direction)
  sig_col <- match.arg(sig_col)

  for (col in c(gene_symbol_col, zscore_col)) {
    if (!col %in% colnames(deg_table)) {
      stop("Cannot find required column in TWAS table: ", col)
    }
  }

  deg_df <- deg_table
  deg_df[[zscore_col]] <- suppressWarnings(as.numeric(deg_df[[zscore_col]]))
  deg_df <- deg_df[!is.na(deg_df[[gene_symbol_col]]) & nzchar(trimws(deg_df[[gene_symbol_col]])) & is.finite(deg_df[[zscore_col]]), , drop = FALSE]

  if (!is.null(padj_col) && nzchar(padj_col) && padj_col %in% colnames(deg_df)) {
    deg_df[[padj_col]] <- suppressWarnings(as.numeric(deg_df[[padj_col]]))
    deg_df <- deg_df[is.finite(deg_df[[padj_col]]) & deg_df[[padj_col]] < padj_cutoff, , drop = FALSE]
  }

  if (direction == "up") {
    deg_df <- deg_df[deg_df[[zscore_col]] > zscore_cutoff, , drop = FALSE]
    deg_df <- deg_df[order(deg_df[[zscore_col]], decreasing = TRUE), , drop = FALSE]
  } else {
    deg_df <- deg_df[deg_df[[zscore_col]] < -zscore_cutoff, , drop = FALSE]
    deg_df <- deg_df[order(deg_df[[zscore_col]], decreasing = FALSE), , drop = FALSE]
  }

  if (remove_ribo && nrow(deg_df) > 0) {
    keep <- !grepl("^(RPL|RPS|MRPL|MRPS|MT-)", deg_df[[gene_symbol_col]])
    deg_df <- deg_df[keep, , drop = FALSE]
  }

  genelist <- deg_df[[zscore_col]]
  names(genelist) <- trimws(deg_df[[gene_symbol_col]])
  genelist <- genelist[!duplicated(names(genelist))]

  if (length(genelist) < min_size) {
    empty_res <- data.frame(
      pathway = character(0), pval = numeric(0), padj = numeric(0),
      log2err = numeric(0), ES = numeric(0), NES = numeric(0),
      size = integer(0), leadingEdge = character(0),
      stringsAsFactors = FALSE
    )
    return(list(deg_used = deg_df, genelist = genelist, fgsea_all = empty_res, fgsea_sig = empty_res))
  }

  fgsea_res <- fgsea::fgsea(
    pathways = pathways,
    stats = genelist,
    minSize = min_size,
    maxSize = max_size,
    nproc = cores
  )
  fgsea_res <- as.data.frame(fgsea_res, stringsAsFactors = FALSE)
  if (nrow(fgsea_res) > 0) {
    fgsea_res$leadingEdge <- vapply(
      fgsea_res$leadingEdge,
      function(x) paste(x, collapse = ", "),
      character(1)
    )
  } else {
    fgsea_res$leadingEdge <- character(0)
  }

  sig_mask <- !is.na(fgsea_res[[sig_col]]) & fgsea_res[[sig_col]] < sig_cutoff
  fgsea_sig <- fgsea_res[sig_mask, , drop = FALSE]
  if (nrow(fgsea_sig) > 0) {
    fgsea_sig <- fgsea_sig[order(fgsea_sig$NES, decreasing = TRUE), , drop = FALSE]
  }

  list(
    deg_used = deg_df,
    genelist = genelist,
    fgsea_all = fgsea_res,
    fgsea_sig = fgsea_sig
  )
}

load_pert_annotation <- function(anno_files) {
  anno_files <- anno_files[nzchar(anno_files)]
  if (length(anno_files) == 0) return(NULL)
  dfs <- lapply(seq_along(anno_files), function(i) {
    path <- anno_files[i]
    if (!file.exists(path)) {
      stop("Perturbation annotation file not found: ", path)
    }
    df <- data.table::fread(path, sep = "\t", data.table = FALSE, showProgress = FALSE)
    keep_cols <- intersect(c("pert_id", "pert_iname", "pert_type", "pubchem_cid"), colnames(df))
    if (!"pert_id" %in% keep_cols) {
      stop("Missing pert_id column in perturbation file: ", path)
    }
    df <- df[, keep_cols, drop = FALSE]
    if (!"pert_iname" %in% colnames(df)) df$pert_iname <- NA_character_
    if (!"pert_type" %in% colnames(df)) df$pert_type <- NA_character_
    if (!"pubchem_cid" %in% colnames(df)) df$pubchem_cid <- NA_character_
    df$pert_id <- as.character(df$pert_id)
    df$pert_iname <- as.character(df$pert_iname)
    df$pert_type <- as.character(df$pert_type)
    df$pubchem_cid <- as.character(df$pubchem_cid)
    df$anno_source <- paste0("source_", i)
    df[, c("pert_id", "pert_iname", "pert_type", "pubchem_cid", "anno_source"), drop = FALSE]
  })
  combined <- do.call(rbind, dfs)
  combined[!duplicated(combined$pert_id), , drop = FALSE]
}

annotate_lincs_fgsea <- function(fgsea_res, pert_anno = NULL, gmt_desc = NULL, pathway_col = "pathway") {
  if (nrow(fgsea_res) == 0) {
    fgsea_res$pert_iname <- character(0)
    fgsea_res$pert_type <- character(0)
    fgsea_res$pubchem_cid <- character(0)
    fgsea_res$anno_source <- character(0)
    return(fgsea_res)
  }
  fgsea_res[[pathway_col]] <- as.character(fgsea_res[[pathway_col]])
  if (!is.null(pert_anno)) {
    idx <- match(fgsea_res[[pathway_col]], pert_anno$pert_id)
    fgsea_res$pert_iname <- pert_anno$pert_iname[idx]
    fgsea_res$pert_type <- pert_anno$pert_type[idx]
    fgsea_res$pubchem_cid <- pert_anno$pubchem_cid[idx]
    fgsea_res$anno_source <- pert_anno$anno_source[idx]
  } else {
    fgsea_res$pert_iname <- NA_character_
    fgsea_res$pert_type <- NA_character_
    fgsea_res$pubchem_cid <- NA_character_
    fgsea_res$anno_source <- "gmt"
  }
  missing_name <- is.na(fgsea_res$pert_iname) | !nzchar(trimws(fgsea_res$pert_iname))
  if (any(missing_name)) {
    if (!is.null(gmt_desc)) {
      desc_vals <- unname(gmt_desc[fgsea_res[[pathway_col]][missing_name]])
      use_desc <- !is.na(desc_vals) & nzchar(trimws(desc_vals))
      fgsea_res$pert_iname[missing_name][use_desc] <- desc_vals[use_desc]
    }
    still_missing <- is.na(fgsea_res$pert_iname) | !nzchar(trimws(fgsea_res$pert_iname))
    fgsea_res$pert_iname[still_missing] <- fgsea_res[[pathway_col]][still_missing]
  }
  fgsea_res
}
