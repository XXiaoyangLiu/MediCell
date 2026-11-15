# Explicit standard/legacy GWAS boundary. Entrypoint stays in script/ root.
self <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))
helper <- file.path(dirname(self), "helpers")
profile <- Sys.getenv("MC_GWAS_PROFILE")
manifest <- Sys.getenv("MC_GWAS_MANIFEST")
if ((nzchar(profile) + nzchar(manifest)) != 1L) stop("Choose exactly one MC_GWAS_PROFILE or MC_GWAS_MANIFEST")
if (nzchar(profile) && !profile %in% c("sle-legacy-pbmc", "sle-legacy-kidney")) stop("Unsupported GWAS profile")
out <- Sys.getenv("MC_GWAS_OUT")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
if (nzchar(profile)) {
  warning("LEGACY REPRODUCTION ONLY: 183 BETA/SE misprojection, awk whitespace column loss, reversed ALLELE labels, invalid SAMPLESIZE command, historical OR/P formula and all-autosome selection are intentionally retained")
  source(file.path(helper, "legacy_sle_align.R"))
} else {
  source(file.path(helper, "gwas_normalize.R"))
}
