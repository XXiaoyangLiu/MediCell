# Shared upstream input boundaries; statistical calculations stay in their callers.
mc_choice <- function(name, choices, default = choices[1]) {
  value <- Sys.getenv(name, default)
  if (!value %in% choices) stop(name, " must be one of: ", paste(choices, collapse = " | "))
  value
}
mc_number <- function(name, default, lower = -Inf, upper = Inf, integer = FALSE) {
  raw <- Sys.getenv(name, default)
  value <- suppressWarnings(as.numeric(raw))
  if (length(value) != 1L || !is.finite(value) || value < lower || value > upper ||
      (integer && (value != floor(value) || value > .Machine$integer.max)))
    stop(name, " must be a finite ", if (integer) "integer" else "number",
         " in [", lower, ", ", upper, "]")
  if (integer) as.integer(value) else value
}
mc_ids <- function(ids, label) {
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    stop(label, " must be nonempty, complete and unique (no duplicate IDs)")
  invisible(ids)
}
mc_sample_key <- function(ids, label) {
  mc_ids(ids, label)
  policy <- mc_choice("MC_SAMPLE_ID_POLICY", c("exact", "legacy-rnames"))
  key <- if (policy == "legacy-rnames") make.names(ids, unique = FALSE) else ids
  if (anyDuplicated(key)) stop(label, ": legacy-rnames make.names collision; IDs must map uniquely")
  key
}
mc_sample_match <- function(expected, available, label) {
  target <- mc_sample_key(expected, "Expression/bulk sample IDs")
  source <- mc_sample_key(available, label)
  idx <- match(target, source)
  if (anyNA(idx)) stop(label, ": missing sample match for ", paste(expected[is.na(idx)], collapse = ", "),
                       "; MC_SAMPLE_ID_POLICY expects exact | legacy-rnames (legacy is explicit)")
  idx
}
mc_numeric <- function(x, label, nonnegative = FALSE) {
  if (!is.numeric(x) || any(!is.finite(x)) || (nonnegative && any(x < 0)))
    stop(label, " must contain only finite", if (nonnegative) " nonnegative" else " numeric", " values")
  invisible(x)
}
mc_matrix <- function(file, label, sep = "") {
  x <- read.table(file, header = TRUE, row.names = 1, sep = sep,
                           check.names = FALSE, stringsAsFactors = FALSE,
                           colClasses = "character", comment.char = "")
  x[] <- lapply(x, type.convert, as.is = TRUE)
  x <- as.matrix(x)
  if (!all(dim(x) > 0)) stop(label, " must have at least one row and column")
  mc_ids(rownames(x), paste(label, "row IDs")); mc_ids(colnames(x), paste(label, "column IDs"))
  mc_numeric(x, label)
  x
}
mc_celltypes <- function(available = NULL) {
  types <- strsplit(trimws(Sys.getenv("MC_CELLTYPES")), "[[:space:]]+")[[1]]
  mc_ids(types, "MC_CELLTYPES (space-separated names)")
  if (!is.null(available) && !all(types %in% available))
    stop("MC_CELLTYPES must name available cell types: ", paste(available, collapse = ", "))
  types
}
mc_load_enigma_inputs <- function(bulk_file, ref_file, fractions_file) {
  # Installed ENIGMA 0.1.5 create_ENIGMA: genes x samples / genes x cell types.
  bulk <- mc_matrix(bulk_file, "Bulk", sep = "\t")
  ref <- mc_matrix(ref_file, "Reference")
  mc_numeric(bulk, "Bulk", TRUE); mc_numeric(ref, "Reference", TRUE)
  if (mc_choice("MC_SAMPLE_ID_POLICY", c("exact", "legacy-rnames")) == "legacy-rnames") {
    mc_sample_key(colnames(bulk), "Bulk sample IDs")
    colnames(bulk) <- make.names(colnames(bulk), unique = FALSE)
  }
  fractions <- read.table(fractions_file, header = TRUE, check.names = FALSE,
                          stringsAsFactors = FALSE, colClasses = "character", comment.char = "")
  mc_ids(colnames(fractions), "Fractions column names")
  if (!"Mixture" %in% colnames(fractions)) stop("Fractions require a Mixture sample ID column")
  mc_ids(fractions$Mixture, "Fractions Mixture IDs")
  if (!all(colnames(ref) %in% colnames(fractions)))
    stop("Fractions must contain every reference cell type: ", paste(colnames(ref), collapse = ", "))
  idx <- mc_sample_match(colnames(bulk), fractions$Mixture, "Fractions Mixture IDs")
  fractions <- as.matrix(fractions[idx, colnames(ref), drop = FALSE])
  suppressWarnings(storage.mode(fractions) <- "double")
  rownames(fractions) <- colnames(bulk)
  mc_numeric(fractions, "Fractions", TRUE)
  if (!length(intersect(rownames(bulk), rownames(ref)))) stop("Bulk and reference require shared gene IDs")
  list(bulk = bulk, ref = ref, fractions = fractions)
}
