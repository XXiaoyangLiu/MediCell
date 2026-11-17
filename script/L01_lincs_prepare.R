# L01: Prepare directional TWAS gene lists, normalized TWAS z-score table,
# and gene x cell-type expression matrix for LINCS enrichment and scoring.
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

pick_col <- function(cols, explicit, candidates, label) {
  if (nzchar(explicit)) {
    if (!explicit %in% cols) {
      stop("Specified ", label, " column not found: ", explicit,
           " (available: ", paste(cols, collapse = ", "), ")")
    }
    return(explicit)
  }
  hit <- candidates[candidates %in% cols]
  if (length(hit) == 0) {
    stop("Could not auto-detect ", label, " column (checked: ",
         paste(candidates, collapse = ", "), "; available: ",
         paste(cols, collapse = ", "), ")")
  }
  hit[1]
}

read_twas_table <- function(path, sheet = "") {
  if (!file.exists(path)) stop("TWAS input file not found: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package 'readxl' is required to read Excel TWAS files: ", path)
    }
    sheets <- readxl::excel_sheets(path)
    target_sheet <- if (nzchar(sheet)) sheet else sheets[1]
    if (!target_sheet %in% sheets) {
      stop("Excel sheet not found: ", target_sheet,
           " (available: ", paste(sheets, collapse = ", "), ")")
    }
    df <- as.data.frame(readxl::read_excel(path, sheet = target_sheet), stringsAsFactors = FALSE)
    known <- c("gene_name", "Symbol", "Gene", "symbol", "zscore", "Zscore", "Z")
    if (!any(known %in% colnames(df))) {
      for (sk in 1:4) {
        df_try <- as.data.frame(readxl::read_excel(path, sheet = target_sheet, skip = sk),
                                stringsAsFactors = FALSE)
        if (any(known %in% colnames(df_try))) {
          df <- df_try
          break
        }
      }
    }
    return(df)
  }
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(data.table::fread(path, data.table = FALSE, showProgress = FALSE, check.names = FALSE))
  }
  if (ext == "csv") {
    return(read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE))
  }
  read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
             check.names = FALSE, quote = "\"", comment.char = "")
}

l01_dir <- req_env("MC_LINCS_L01_DIR")
trait <- req_env("MC_LINCS_TRAIT")
label <- req_env("MC_LINCS_LABEL")
direction <- req_env("MC_LINCS_DIRECTION")
twas_file <- req_env("MC_LINCS_TWAS_FILE")
twas_sheet <- opt_env("MC_LINCS_TWAS_SHEET")
cell_types_raw <- opt_env("MC_CELLTYPES")
gene_col_in <- opt_env("MC_LINCS_GENE_COL")
z_col_in <- opt_env("MC_LINCS_Z_COL")
padj_col_in <- opt_env("MC_LINCS_PADJ_COL")
ct_col_in <- opt_env("MC_LINCS_CT_COL")
exp_matrix_in <- opt_env("MC_LINCS_EXP_MATRIX")
enigma_rdata <- opt_env("MC_ENIGMA_RDATA")
gene_pos_file <- opt_env("MC_GENE_POS")

dir.create(l01_dir, recursive = TRUE, showWarnings = FALSE)

twas_df <- read_twas_table(twas_file, twas_sheet)
if (nrow(twas_df) == 0) stop("TWAS table is empty: ", twas_file)

cols <- colnames(twas_df)
gene_col <- pick_col(cols, gene_col_in, c("gene_name", "Symbol", "Gene", "symbol"), "gene symbol")
z_col <- pick_col(cols, z_col_in, c("zscore", "Zscore", "Z"), "z-score")
ct_col <- pick_col(cols, ct_col_in, c("cell_type", "Celltype", "CellType", "celltype"), "cell type")

padj_col <- ""
if (nzchar(padj_col_in)) {
  padj_col <- pick_col(cols, padj_col_in, c(padj_col_in), "adjusted p-value")
} else {
  cand <- c("fdr_p", "FDR P.value", "padj", "p_val_adj", "FDR", "bonferroni_p")
  hit <- cand[cand %in% cols]
  if (length(hit) > 0) padj_col <- hit[1]
}

twas_df[[gene_col]] <- trimws(as.character(twas_df[[gene_col]]))
twas_df[[ct_col]] <- trimws(as.character(twas_df[[ct_col]]))
twas_df[[z_col]] <- suppressWarnings(as.numeric(twas_df[[z_col]]))
if (nzchar(padj_col)) {
  twas_df[[padj_col]] <- suppressWarnings(as.numeric(twas_df[[padj_col]]))
}

twas_df <- twas_df[
  !is.na(twas_df[[gene_col]]) & nzchar(twas_df[[gene_col]]) &
    !is.na(twas_df[[ct_col]]) & nzchar(twas_df[[ct_col]]) &
    is.finite(twas_df[[z_col]]),
  ,
  drop = FALSE
]

min_cts <- suppressWarnings(as.integer(opt_env("MC_LINCS_MIN_CTS", "2")))
if (is.na(min_cts) || min_cts < 1L) min_cts <- 2L
cell_types <- if (nzchar(cell_types_raw)) strsplit(trimws(cell_types_raw), "[[:space:]]+")[[1]] else character(0)
if (length(cell_types) > 0) {
  twas_df <- twas_df[twas_df[[ct_col]] %in% cell_types, , drop = FALSE]
  if (nrow(twas_df) == 0) {
    stop("No TWAS rows match requested --cell-types: ", paste(cell_types, collapse = ", "))
  }
}
n_twas_cts <- length(unique(twas_df[[ct_col]]))
if (min_cts >= 2L && n_twas_cts < min_cts) {
  stop("LINCS scoring requires at least ", min_cts,
       " cell types in TWAS input for cross-cell-type normalization; found ",
       n_twas_cts, " (", paste(unique(twas_df[[ct_col]]), collapse = ", "), ").")
}

norm_twas <- data.frame(
  Gene = twas_df[[gene_col]],
  CellType = twas_df[[ct_col]],
  Zscore = twas_df[[z_col]],
  stringsAsFactors = FALSE
)
norm_twas_path <- file.path(l01_dir, paste0(trait, "_", label, "_twas_all.tsv"))
write.table(norm_twas, file = norm_twas_path, sep = "\t", row.names = FALSE, quote = FALSE)

if (direction == "pos") {
  dir_df <- twas_df[twas_df[[z_col]] > 0, , drop = FALSE]
} else if (direction == "neg") {
  dir_df <- twas_df[twas_df[[z_col]] < 0, , drop = FALSE]
} else {
  stop("Invalid MC_LINCS_DIRECTION: ", direction)
}

std_dir <- data.frame(
  Symbol = dir_df[[gene_col]],
  Zscore = dir_df[[z_col]],
  Celltype = dir_df[[ct_col]],
  stringsAsFactors = FALSE
)
if (nzchar(padj_col)) {
  std_dir$FDR_P <- dir_df[[padj_col]]
}
dir_csv_path <- file.path(l01_dir, paste0(trait, "_", label, "_", direction, ".csv"))
write.csv(std_dir, file = dir_csv_path, row.names = FALSE)
message("Wrote directional TWAS table (", nrow(std_dir), " rows) -> ", dir_csv_path)

exp_out_path <- file.path(l01_dir, paste0(label, "_exp_matrix.csv"))
if (nzchar(exp_matrix_in)) {
  if (!file.exists(exp_matrix_in)) stop("Expression matrix file not found: ", exp_matrix_in)
  exp_df <- read.csv(exp_matrix_in, header = TRUE, na.strings = c("", "NA"),
                     row.names = 1, check.names = FALSE)
  if (length(cell_types) > 0) {
    missing_ct <- setdiff(cell_types, colnames(exp_df))
    if (length(missing_ct) > 0) {
      stop("Cell types missing from --exp-matrix: ", paste(missing_ct, collapse = ", "))
    }
    exp_df <- exp_df[, cell_types, drop = FALSE]
  }
  if (min_cts >= 2L && ncol(exp_df) < min_cts) {
    stop("LINCS scoring requires at least ", min_cts,
         " cell types in expression matrix; found ", ncol(exp_df), ".")
  }
  write.csv(exp_df, file = exp_out_path, row.names = TRUE)
  message("Prepared expression matrix (", nrow(exp_df), " genes x ", ncol(exp_df), " cell types) -> ", exp_out_path)
} else {
  if (!nzchar(enigma_rdata) || !file.exists(enigma_rdata)) {
    stop("Provide --exp-matrix or a valid M01 ENIGMA RData file (--enigma-rdata).")
  }
  env_e <- new.env(parent = emptyenv())
  load(enigma_rdata, envir = env_e)
  if (!exists("enigma.cse.test", envir = env_e)) {
    stop("Object 'enigma.cse.test' not found in ", enigma_rdata)
  }
  cse <- get("enigma.cse.test", envir = env_e)
  if (length(dim(cse)) != 3L) {
    stop("Expected 3D array enigma.cse.test in ", enigma_rdata)
  }
  if (length(cell_types) > 0) {
    missing_ct <- setdiff(cell_types, dimnames(cse)[[3]])
    if (length(missing_ct) > 0) {
      stop("Cell types missing from enigma.cse.test: ", paste(missing_ct, collapse = ", "))
    }
    cse <- cse[, , cell_types, drop = FALSE]
  }
  exp_mat <- apply(cse, c(1, 3), mean, na.rm = TRUE)
  gene_ids <- rownames(exp_mat)

  id_to_sym <- NULL
  if (nzchar(gene_pos_file)) {
    if (!file.exists(gene_pos_file)) stop("gene_pos file not found: ", gene_pos_file)
    gp <- read.csv(gene_pos_file, stringsAsFactors = FALSE, check.names = FALSE)
    if (all(c("gene_id", "gene_name") %in% colnames(gp))) {
      id_to_sym <- setNames(gp$gene_name, gp$gene_id)
    }
  } else if (all(c("gene_id", "gene_name") %in% colnames(twas_df))) {
    id_to_sym <- setNames(twas_df$gene_name, twas_df$gene_id)
  }

  if (!is.null(id_to_sym)) {
    direct_hits <- sum(gene_ids %in% names(id_to_sym))
    if (direct_hits > 0) {
      mapped <- unname(id_to_sym[gene_ids])
    } else {
      stripped_ids <- sub("\\.[0-9]+$", "", gene_ids)
      names(id_to_sym) <- sub("\\.[0-9]+$", "", names(id_to_sym))
      id_to_sym <- id_to_sym[!duplicated(names(id_to_sym))]
      mapped <- if (any(stripped_ids %in% names(id_to_sym))) unname(id_to_sym[stripped_ids]) else NULL
    }
    if (!is.null(mapped)) {
      keep <- !is.na(mapped) & nzchar(trimws(mapped))
      exp_mat <- exp_mat[keep, , drop = FALSE]
      mapped <- trimws(mapped[keep])
      if (anyDuplicated(mapped)) {
        summed <- rowsum(exp_mat, group = mapped, reorder = FALSE)
        counts <- as.numeric(table(mapped)[rownames(summed)])
        exp_mat <- summed / counts
      } else {
        rownames(exp_mat) <- mapped
      }
    }
  }

  if (min_cts >= 2L && ncol(exp_mat) < min_cts) {
    stop("LINCS scoring requires at least ", min_cts,
         " cell types in aggregated expression matrix; found ", ncol(exp_mat), ".")
  }
  write.csv(exp_mat, file = exp_out_path, row.names = TRUE)
  message("Aggregated ENIGMA expression matrix (", nrow(exp_mat), " genes x ", ncol(exp_mat), " cell types) -> ", exp_out_path)
}
