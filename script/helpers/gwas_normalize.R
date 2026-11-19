# Standard GWAS contract: no strand guessing, liftover, or formula evaluation.
# Unsupported rows are dropped with explicit counts; declared SE<=0 stays a hard error.
library(data.table)
fail <- function(...) stop(..., call.=FALSE)
read_tsv <- function(p, na.strings=c("NA", "#NA"), ...) fread(p, sep="\t", header=TRUE, quote="", na.strings=na.strings, ...)
manifest_path <- Sys.getenv("MC_GWAS_MANIFEST")
if (!startsWith(manifest_path, "/")) fail("manifest path must be absolute")
m <- read_tsv(manifest_path, colClasses="character", na.strings=NULL)
fields <- strsplit("study_id path build effect_scale chr_col pos_col ea_col oa_col effect_mode effect_col se_col p_col eaf_col n_mode n_value n_col", " ")[[1]]
if (!identical(names(m), fields)) fail("manifest columns must be: ", paste(fields, collapse="\t"))
if (!nrow(m) || anyDuplicated(m$study_id) || any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", m$study_id))) fail("study_id must be unique safe labels")
build <- Sys.getenv("MC_GENOME_BUILD")
if (!build %in% c("GRCh37", "GRCh38") || any(m$build != build)) fail("build mismatch (no liftover)")
if (length(unique(m$effect_scale)) != 1L || any(!m$effect_scale %in% c("linear", "log_odds"))) fail("effect_scale must agree (linear or log_odds)")
if (any(!m$effect_mode %in% c("beta_se", "or_se", "or_p")) || any(m$effect_mode != "beta_se" & m$effect_scale != "log_odds")) fail("effect_mode incompatible with effect_scale")
chrs <- strsplit(gsub(",", " ", Sys.getenv("MC_CHROMS", "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22")), " +")[[1]]
if (!length(chrs) || any(!chrs %in% as.character(1:22))) fail("MC_CHROMS must select autosomes")
panel <- read_tsv(Sys.getenv("MC_SNP_ANNOT"))
if (!all(c("chromosome", "Pos", "VariantID", "ref_vcf", "alt_vcf", "RSID") %in% names(panel))) fail("panel columns missing")
panel <- panel[as.character(chromosome) %in% chrs, .(chr=as.character(chromosome), pos=as.numeric(Pos), SNP=as.character(VariantID), REF=toupper(ref_vcf), ALT=toupper(alt_vcf))]
if (!nrow(panel)) fail("empty panel for selected chromosomes")
if (anyNA(panel) || any(!is.finite(panel$pos) | panel$pos <= 0 | panel$pos != floor(panel$pos)) ||
    any(!nzchar(panel$SNP))) fail("missing/invalid panel coordinates or keys")
dup <- duplicated(panel[, .(chr,pos)]) | duplicated(panel$SNP)
if (any(dup)) fail("duplicate panel locus/key (ambiguous reference)")
bad <- !panel$REF %in% c("A","C","G","T") | !panel$ALT %in% c("A","C","G","T") | panel$REF == panel$ALT
panel_dropped <- sum(bad); panel <- panel[!bad]
canonical <- paste(panel$chr, sprintf("%.0f", panel$pos), panel$REF, panel$ALT, sep="_")
if (any(panel$SNP != canonical)) fail("panel VariantID must be canonical chr_pos_ref_alt")
if (!nrow(panel)) fail("no usable panel variants")
pal <- function(a,b) paste0(a,b) %in% c("AT","TA","CG","GC")
out <- Sys.getenv("MC_GWAS_OUT")
studies <- data.table(study_id=character(), path=character())
counts <- list()
for (i in seq_len(nrow(m))) {
  s <- as.list(m[i])
  if (!startsWith(s$path, "/")) fail("study path must be absolute: ", s$study_id)
  if (any(!nzchar(unlist(s[c("chr_col","pos_col","ea_col","oa_col","effect_col")])))) fail("missing required column mapping")
  if (s$effect_mode %in% c("beta_se","or_se") && !nzchar(s$se_col)) fail("SE mapping required")
  if (s$effect_mode == "or_p" && !nzchar(s$p_col)) fail("P mapping required for explicit or_p mode")
  if (!s$n_mode %in% c("none","constant","column") || (s$n_mode=="constant" && (!nzchar(s$n_value) || nzchar(s$n_col))) || (s$n_mode=="column" && (!nzchar(s$n_col) || nzchar(s$n_value))) || (s$n_mode=="none" && (nzchar(s$n_value)||nzchar(s$n_col)))) fail("invalid n_mode/n_value/n_col")
  mapped <- unlist(s[c("chr_col","pos_col","ea_col","oa_col","effect_col","se_col","p_col","eaf_col","n_col")], use.names=FALSE)
  cols <- unique(mapped[nzchar(mapped)])
  header <- names(read_tsv(s$path, na.strings=NULL, nrows=0))
  if (any(!cols %in% header)) fail("missing mapped study columns: ", s$study_id)
  a <- read_tsv(s$path, na.strings=NULL, select=cols)
  input_n <- nrow(a)
  a <- a[as.character(a[[s$chr_col]]) %in% chrs]
  number <- function(col) suppressWarnings(as.numeric(a[[col]]))
  optional_number <- function(col, label, lower, upper, strict_lower=FALSE) {
    if (!nzchar(col)) return(rep(NA_real_, nrow(a)))
    raw <- trimws(as.character(a[[col]]))
    missing <- is.na(raw) | raw %in% c("", "NA", "#NA", ".")
    value <- suppressWarnings(as.numeric(raw))
    invalid <- !missing & (!is.finite(value) | value < lower | value > upper |
                          (strict_lower & value == lower))
    if (any(invalid)) fail("invalid nonmissing ", label, ": ", s$study_id)
    value[missing] <- NA_real_
    value
  }
  p_value <- optional_number(s$p_col, "P", 0, 1, s$effect_mode == "or_p")
  frequency <- optional_number(s$eaf_col, "EAF", 0, 1)
  sample_n <- if (s$n_mode == "column") optional_number(s$n_col, "N", 0, Inf, TRUE) else NA_real_
  if (s$n_mode == "constant") {
    sample_n <- suppressWarnings(as.numeric(s$n_value))
    if (!is.finite(sample_n) || sample_n <= 0) fail("constant N must be finite and positive")
  }
  # Attach every numeric column before any row drops, so counts and data stay aligned.
  d <- data.table(chr=as.character(a[[s$chr_col]]), pos=number(s$pos_col),
                  EA=toupper(as.character(a[[s$ea_col]])), OA=toupper(as.character(a[[s$oa_col]])),
                  BETA=number(s$effect_col),
                  SE=if (nzchar(s$se_col)) number(s$se_col) else NA_real_,
                  P=p_value, FREQ=frequency, N=rep(sample_n, length.out=nrow(a)))
  if (s$effect_mode != "or_p" && any(!is.finite(d$SE) | d$SE <= 0))
    fail("SE must be finite and positive: ", s$study_id)
  keep <- is.finite(d$pos) & d$pos > 0 & d$pos == floor(d$pos)
  dropped_pos <- sum(!keep); d <- d[keep]
  keep <- !is.na(d$EA) & !is.na(d$OA) & d$EA %in% c("A","C","G","T") & d$OA %in% c("A","C","G","T") & d$EA != d$OA
  dropped_alleles <- sum(!keep); d <- d[keep]
  keep <- is.finite(d$BETA)
  dropped_beta <- sum(!keep); d <- d[keep]
  if (s$effect_mode != "beta_se") {
    keep <- d$BETA > 0
    dropped_or_p <- sum(!keep); d <- d[keep]
    d[, BETA:=log(BETA)]
  } else dropped_or_p <- 0L
  if (anyDuplicated(d[, .(chr,pos)])) fail("duplicate study locus (including conflicting alleles/effects)")
  n_pal <- sum(pal(d$EA, d$OA)); d <- d[!pal(EA,OA)]
  if (s$effect_mode == "or_p") {
    keep <- is.finite(d$P) & d$P > 0 & d$P <= 1
    dropped_or_p <- dropped_or_p + sum(!keep)
    d <- d[keep]
    # lower tail + log P avoids 1-P/2 cancellation, including 1e-300.
    z <- abs(qnorm(log(d$P) - log(2), lower.tail=TRUE, log.p=TRUE))
    d[, SE:=abs(BETA)/z]
    keep <- is.finite(d$SE) & d$SE > 0
    dropped_or_p <- dropped_or_p + nrow(d) - sum(keep); d <- d[keep]
  } else if (any(!is.finite(d$SE) | d$SE <= 0)) {
    fail("SE must be finite and positive: ", s$study_id)
  }
  invalid_p <- 0L; invalid_freq <- 0L; invalid_n <- 0L
  before_join <- nrow(d)
  d <- merge(d, panel, by=c("chr","pos"), sort=TRUE)
  keep <- (d$EA==d$ALT & d$OA==d$REF) | (d$EA==d$REF & d$OA==d$ALT)
  unmatched <- before_join - nrow(d) + sum(!keep)
  d <- d[keep]
  flip <- d$EA==d$REF
  d[flip, `:=`(BETA=-BETA, FREQ=1-FREQ)]
  normalized <- d[, .(SNP,A1=ALT,A2=REF,FREQ,BETA,SE,P,N)]
  if (!nrow(normalized)) fail("no matched variants: ", s$study_id)
  dest <- file.path(out, paste0(s$study_id, "_normalized.tsv"))
  fwrite(normalized, dest, sep="\t", na="NA", quote=FALSE)
  studies <- rbind(studies, data.table(study_id=s$study_id, path=normalizePath(dest)))
  counts[[i]] <- data.table(study_id=s$study_id, effect_mode=s$effect_mode, input=input_n, selected=nrow(a), dropped_pos, dropped_alleles, dropped_beta, dropped_or_p, palindromic_dropped=n_pal, unmatched, matched=nrow(normalized), flipped=sum(flip), invalid_p, invalid_freq, invalid_n, panel_dropped=0L)
  rm(a,d,normalized); gc(FALSE)
}
fwrite(studies, file.path(out,"gwas_studies.tsv"), sep="\t", quote=FALSE)
cc <- rbindlist(counts)
if (panel_dropped) cc <- rbind(cc, data.table(study_id="_panel", effect_mode="", input=0, selected=0, dropped_pos=0, dropped_alleles=0, dropped_beta=0, dropped_or_p=0, palindromic_dropped=0, unmatched=0, matched=0, flipped=0, invalid_p=0, invalid_freq=0, invalid_n=0, panel_dropped=panel_dropped))
fwrite(cc, file.path(out,"gwas_normalization_counts.tsv"), sep="\t")
single_out <- Sys.getenv("MC_GWAS_SINGLE_OUT")
if (nzchar(single_out)) {
  if (nrow(m) != 1L) fail("single-study direct GWAS output requires exactly 1 study in manifest")
  norm1 <- fread(studies$path[1], sep="\t", header=TRUE)
  pval <- norm1$P
  miss_p <- is.na(pval)
  if (any(miss_p)) {
    pval[miss_p] <- 2 * pnorm(-abs(norm1$BETA[miss_p] / norm1$SE[miss_p]))
  }
  direct <- data.table(
    MarkerName=norm1$SNP,
    Allele1=norm1$A1,
    Allele2=norm1$A2,
    Effect=norm1$BETA,
    StdErr=norm1$SE,
    `P-value`=pval
  )
  fwrite(direct, single_out, sep="\t", na="NA", quote=FALSE)
}
