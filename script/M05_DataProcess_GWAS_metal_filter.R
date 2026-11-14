library(data.table)
profile <- Sys.getenv("MC_GWAS_PROFILE")
manifest <- Sys.getenv("MC_GWAS_MANIFEST")
if ((nzchar(profile) + nzchar(manifest)) != 1L) stop("Choose exactly one GWAS mode")
if (nzchar(profile) && !profile %in% c("sle-legacy-pbmc", "sle-legacy-kidney")) stop("Unsupported profile")
policy <- Sys.getenv("MC_GWAS_FILTER_POLICY", if (nzchar(profile)) "legacy" else "none")
if (!policy %in% c("legacy", "none") || (nzchar(profile) && policy != "legacy")) stop("Invalid filter policy")
a <- fread(Sys.getenv("MC_METAL_TBL"))
required <- c("MarkerName", "Allele1", "Allele2", "Direction", "HetISq", "HetPVal")
if (!all(required %in% names(a)) || nrow(a) == 0) stop("Invalid/empty METAL result")
if (profile == "sle-legacy-pbmc" || !nzchar(profile)) {
  a[, `:=`(Allele1=toupper(Allele1), Allele2=toupper(Allele2))]
}
if (policy == "legacy") {
  a <- a[!(grepl("\\+", Direction) & grepl("\\-", Direction))]
  a <- a[!(HetISq > 50 & HetPVal < 0.05)]
}
if (nrow(a) == 0) stop("No variants remain after GWAS filtering")
fwrite(a, Sys.getenv("MC_GWAS_FILT"), sep="\t")
