#!/usr/bin/env bash
# ==============================================================================
# MediCell Manuscript Workflow Runner
# ==============================================================================
#
# WHAT THIS SCRIPT DOES:
#   Runs the complete MediCell TWAS (`medicell.py`) and LINCS therapeutic
#   prioritization (`medicell.py lincs`) analyses reported in the manuscript.
#   Only MediCell pipeline workflows are included (no plotting or unrelated
#   downstream scripts).
#
#   Workflow scope (GRCh37 autosomes 1-22):
#     1. Stage 1 (`prepare`):
#        - Exports kidney covariates (Glom, Tub) from QTLData.RData.
#        - Removes multiallelic / duplicate (chr, pos) loci from SNP panels,
#          per-chromosome SNP annotations, and GWAS summary statistics so all
#          inputs satisfy MediCell strict biallelic contracts.
#        - Writes resolved absolute-path GWAS manifest TSVs.
#     2. Stage 2 (`twas`, 3 tissues x 9 traits = 27 TWAS runs):
#        - Tissues : PBMC (6 cell types), Glom (9 cell types), Tub (12 cell types)
#        - Traits  : SLE (3-study EUR+EAS meta-analysis, `EUR_EAS` LD),
#                    eGFR (GCST90691925, `EUR` LD),
#                    7 UKB traits: Glu, UA, uAlb, uCr, uK, uNa, SBP (`EUR` LD)
#        - Compute reuse: each tissue trains expression deconvolution, PEER
#          covariates, and PredictDB models once on `SLE`, then symlinks shared
#          M01/M02/M03/M06 directories for `eGFR` (which builds `EUR` LD once)
#          and reuses both models and `M05_ld/cov_EUR` across the 7 UKB traits.
#     3. Stage 3 (`lincs`, 3 tissues x 2 traits = 6 LINCS runs):
#        - Prioritizes therapeutic compounds (`fgsea` against Broad L1000) for
#          `SLE` (`--direction pos`) and `eGFR` (`--direction neg`) per tissue.
#
# PREREQUISITES:
#   - Conda environments required by medicell.py: `R4.4.3`, `r3.6.3`, `imlatools`
#   - Unpacked LINCS reference data under `data/ref/lincs/`
#   - Access to cohort & GWAS source files under LIUXY_ROOT (/node1/liuxy)
#     and ZHAOQ_ROOT (/node1/zhaoq)
#
# HOW TO RUN (QUICK START):
#   # 1) Inspect planned commands without modifying files or launching jobs:
#   bash Results/run.sh --dry-run
#
#   # 2) Run the full manuscript pipeline in the background (recommended):
#   nohup bash Results/run.sh --cores 16 > Results/run.log 2>&1 &
#
#   # 3) Run a subset of tissues or a single stage:
#   bash Results/run.sh --tissue PBMC --stage twas --cores 16
#   bash Results/run.sh --tissue Glom,Tub --stage lincs --cores 8
#
# OUTPUT LAYOUT & RESUME BEHAVIOR:
#   <out-dir>/                  (default: Results/out, git-ignored)
#     prepared/                 staged covariates, biallelic SNP/GWAS files, manifests
#     ld_dosage/                extracted 1000G dosages per tissue and ancestry panel
#     runs/<Tissue>/<Trait>/    M01-M08 TWAS outputs and (for SLE/eGFR) L01-L03 LINCS outputs
#
#   Safe to re-run: completed prepared files, TWAS summaries
#   (`M08_summary/allcelltype_result.csv`), and LINCS score matrices
#   (`L03_drug_scores/<Tissue>_<Trait>_<Dir>.drug_celltype_score.csv`) are
#   detected and skipped automatically.
# ==============================================================================
set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

# ------------------------------------------------------------------------------
# 1. CLI parsing & validation
# ------------------------------------------------------------------------------
CORES="${MC_CORES:-16}"
OUT_ROOT="$ROOT/Results/out"
TISSUE_ARG="PBMC,Glom,Tub"
STAGE="all"
DRY_RUN=0

usage() {
    cat <<USAGE
Usage: bash Results/run.sh [OPTIONS]

Run the manuscript MediCell TWAS and LINCS pipelines (no plotting code).

Options:
  --cores N       CPU cores per step (default: \$MC_CORES or 16)
  --out-dir DIR   Output root directory (default: Results/out)
  --tissue LIST   Comma-separated tissues from PBMC,Glom,Tub (default: PBMC,Glom,Tub)
  --stage STAGE   Stage to run: all | prepare | twas | lincs (default: all)
  --dry-run       Print planned commands without creating files or executing jobs
  -h, --help      Show this help message

Examples:
  bash Results/run.sh --dry-run
  nohup bash Results/run.sh --cores 16 > Results/run.log 2>&1 &
  bash Results/run.sh --tissue PBMC --stage twas --cores 16
  bash Results/run.sh --stage lincs --cores 8
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cores)   CORES="${2:?--cores requires a value}"; shift 2 ;;
        --out-dir) OUT_ROOT="${2:?--out-dir requires a value}"; shift 2 ;;
        --tissue)  TISSUE_ARG="${2:?--tissue requires a value}"; shift 2 ;;
        --stage)   STAGE="${2:?--stage requires a value}"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

[[ "$CORES" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid --cores (must be positive integer): $CORES" >&2; exit 1; }
case "$STAGE" in
    all|prepare|twas|lincs) ;;
    *) echo "Invalid --stage (expected all|prepare|twas|lincs): $STAGE" >&2; exit 1 ;;
esac

SELECTED_TISSUES=()
IFS=, read -r -a _raw_tissues <<< "$TISSUE_ARG"
for _t in "${_raw_tissues[@]}"; do
    case "$_t" in
        PBMC|Glom|Tub)
            _seen=0
            for _prev in ${SELECTED_TISSUES[@]+"${SELECTED_TISSUES[@]}"}; do
                [[ "$_prev" == "$_t" ]] && _seen=1
            done
            [[ "$_seen" -eq 0 ]] && SELECTED_TISSUES+=("$_t")
            ;;
        *)
            echo "Invalid tissue in --tissue (expected PBMC, Glom, or Tub): $_t" >&2
            exit 1
            ;;
    esac
done
[[ "${#SELECTED_TISSUES[@]}" -gt 0 ]] || { echo "--tissue must select at least one tissue" >&2; exit 1; }

case "$OUT_ROOT" in
    /*) ;;
    *) OUT_ROOT="$ROOT/$OUT_ROOT" ;;
esac

PREP_DIR="$OUT_ROOT/prepared"
RUNS_DIR="$OUT_ROOT/runs"
LD_DIR="$OUT_ROOT/ld_dosage"

# ------------------------------------------------------------------------------
# 2. Manuscript data roots & shared model/LINCS parameters
# ------------------------------------------------------------------------------
LIUXY_ROOT="${LIUXY_ROOT:-/node1/liuxy}"
ZHAOQ_ROOT="${ZHAOQ_ROOT:-/node1/zhaoq}"

G1000_VCF_PATTERN="$LIUXY_ROOT/index/1000G_2/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
SAMPLES_EUR_EAS="$LIUXY_ROOT/index/1000G/EUR_EAS_samples.txt"
SAMPLES_EUR="$LIUXY_ROOT/index/1000G/EUR_samples.txt"
GENE_POS="$LIUXY_ROOT/index/gene_pos_grch37.82.csv"

UKB_TRAITS=(Glu UA uAlb uCr uK uNa SBP)

# Shared ENIGMA / normalization / elastic-net / PEER thresholds across all 3 tissues
SHARED_MODEL_ARGS=(
    --alpha 0.1               # ENIGMA trace-norm regularization weight
    --tpm-min 0.1             # Minimum deconvolved TPM threshold per sample
    --min-pct 0.01            # Retain genes expressed in >= 1% of samples
    --min-samples 5           # Floor on expressing sample count
    --min-level-count 5       # Minimum level size for categorical covariate dummies
    --maf 0.01                # Cis-genotype MAF filter for PredictDB training
    --maf-policy filter
)

# Shared fgsea & cell-type drug prioritization settings across SLE and eGFR
SHARED_LINCS_ARGS=(
    --twas-padj-cutoff 0.05   # Retain FDR-significant TWAS genes
    --twas-zscore-cutoff 0.25 # Minimum |TWAS z-score| filter
    --lincs-min-size 10       # fgsea gene-set minimum overlap
    --lincs-max-size 500      # fgsea gene-set maximum size
    --lincs-sig-col pval      # Filter enrichments by nominal p-value
    --lincs-sig-cutoff 0.05
    --lincs-nes-filter neg    # Keep negatively enriched signatures (disease reversal)
    --lincs-rank-by pval
    --top-n 50                # Top perturbations scored per run
)

# ------------------------------------------------------------------------------
# 3. Runtime helpers
# ------------------------------------------------------------------------------
log() {
    printf "[%s] [Results] %s\n" "$(date "+%F %T")" "$*"
}

run_cmd() {
    if [[ "$DRY_RUN" -eq 1 ]]; then
        printf "[DRY-RUN]"
        printf " %q" "$@"
        printf "\n"
        return 0
    fi
    "$@"
}

# Portable R 4.4.3 runner: prefers $HOME/miniforge3/envs/R4.4.3/bin/Rscript,
# falls back to `conda run -n R4.4.3`, then PATH `Rscript`.
run_r4() {
    local direct="$HOME/miniforge3/envs/R4.4.3/bin/Rscript"
    if [[ -x "$direct" ]]; then
        "$direct" "$@"
        return
    fi
    local conda_bin=""
    if [[ -n "${CONDA_EXE:-}" && -x "${CONDA_EXE}" ]]; then
        conda_bin="$CONDA_EXE"
    elif [[ -x "$HOME/miniforge3/bin/conda" ]]; then
        conda_bin="$HOME/miniforge3/bin/conda"
    elif command -v conda >/dev/null 2>&1; then
        conda_bin="$(command -v conda)"
    fi
    if [[ -n "$conda_bin" ]]; then
        "$conda_bin" run -n R4.4.3 --no-capture-output Rscript "$@"
    elif command -v Rscript >/dev/null 2>&1; then
        Rscript "$@"
    else
        echo "Cannot locate R 4.4.3 Rscript binary" >&2
        exit 1
    fi
}

# Symlink immutable upstream tissue outputs into downstream trait run directories
# so ENIGMA, PEER, PredictDB training, SQLite DBs, and EUR LD run only once.
link_shared_stages() {
    local src_dir="$1" dst_dir="$2" include_eur_ld="${3:-0}"
    if [[ "$DRY_RUN" -eq 1 ]]; then
        if [[ "$include_eur_ld" -eq 1 ]]; then
            echo "[DRY-RUN] link_shared_stages $src_dir -> $dst_dir (+cov_EUR)"
        else
            echo "[DRY-RUN] link_shared_stages $src_dir -> $dst_dir"
        fi
        return 0
    fi
    mkdir -p "$dst_dir"
    local sub
    for sub in M01_deconvolution M02_expression M02_covariates M03_models M06_databases; do
        ln -sfn "$src_dir/$sub" "$dst_dir/$sub"
    done
    if [[ "$include_eur_ld" -eq 1 ]]; then
        mkdir -p "$dst_dir/M05_ld"
        ln -sfn "$src_dir/M05_ld/cov_EUR" "$dst_dir/M05_ld/cov_EUR"
    fi
}

# ------------------------------------------------------------------------------
# 4. Stage 1: Prepare covariates, biallelic SNP annotations, and GWAS manifests
# ------------------------------------------------------------------------------
prepare_inputs() {
    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "[DRY-RUN] prepare_inputs -> $PREP_DIR"
        return 0
    fi

    log "Staging covariates, biallelic SNP annotations, and GWAS manifests under $PREP_DIR"
    mkdir -p "$PREP_DIR/covariates" "$PREP_DIR/snp_annot/pbmc" "$PREP_DIR/snp_annot/glom" \
             "$PREP_DIR/snp_annot/tub" "$PREP_DIR/gwas" "$PREP_DIR/manifests"

    MC_PREP_CORES="$CORES" MC_PREP_DIR="$PREP_DIR" LIUXY_ROOT="$LIUXY_ROOT" ZHAOQ_ROOT="$ZHAOQ_ROOT" \
    run_r4 - <<'RS'
options(scipen = 999)
suppressPackageStartupMessages(library(data.table))

prep_dir <- Sys.getenv("MC_PREP_DIR")
liuxy    <- Sys.getenv("LIUXY_ROOT")
zhaoq    <- Sys.getenv("ZHAOQ_ROOT")
cores    <- suppressWarnings(as.integer(Sys.getenv("MC_PREP_CORES", "8")))
if (is.na(cores) || cores < 1L) cores <- 8L
setDTthreads(min(cores, 16L))

cov_dir <- file.path(prep_dir, "covariates")
snp_dir <- file.path(prep_dir, "snp_annot")
gwas_dir <- file.path(prep_dir, "gwas")
man_dir <- file.path(prep_dir, "manifests")

# A. Export kidney known covariates from QTLData.RData.
# Note: rownames(df) match bulk/CIBERSORTx/genotype sample IDs; df$Sample differs and is ignored.
export_qtl_cov <- function(rdata_path, obj_name, cols, out_path) {
  if (file.exists(out_path) && file.info(out_path)$size > 0) return(invisible(NULL))
  e <- new.env(parent = emptyenv())
  load(rdata_path, envir = e)
  df <- get(obj_name, envir = e)
  out <- data.frame(Sample = rownames(df), df[, cols, drop = FALSE],
                    stringsAsFactors = FALSE, check.names = FALSE)
  tmp <- paste0(out_path, ".tmp")
  fwrite(out, tmp, sep = "\t", quote = FALSE, na = "NA", compress = "none")
  file.rename(tmp, out_path)
}

export_qtl_cov(
  file.path(zhaoq, "scTWAS_Kidney/predictDB/kidney_KS/glom/QTLData.RData"),
  "AAglomPhe2",
  c("Age", "RIN", "PC1", "PC2", "PC3", "PC4", "PC5", "Sclerosis", "Pairend", "Gender", "Batch"),
  file.path(cov_dir, "glom_covariates.tsv")
)
export_qtl_cov(
  file.path(zhaoq, "scTWAS_Kidney/predictDB/kidney_KS/tub/QTLData.RData"),
  "AAtubPhe2",
  c("Age", "RIN", "PC1", "PC2", "PC3", "PC4", "PC5", "Fibrosis", "Pairend", "Gender", "Batch"),
  file.path(cov_dir, "tub_covariates.tsv")
)

# B. Deduplicate SNP panels and chr1..22 annotations by (chr, pos).
# MediCell ld_contracts.py and gwas_normalize.R fail closed on multiallelic/duplicate loci;
# dropping positions with >1 record leaves unambiguous biallelic SNPs only.
dedup_tsv <- function(src, dst, chr_col, pos_col) {
  if (file.exists(dst) && file.info(dst)$size > 0) return(invisible(NULL))
  dt <- fread(src, sep = "\t", header = TRUE, quote = "", na.strings = NULL, colClasses = "character")
  key <- paste0(dt[[chr_col]], ":", dt[[pos_col]])
  dt <- dt[!(duplicated(key) | duplicated(key, fromLast = TRUE))]
  tmp <- paste0(dst, ".tmp")
  fwrite(dt, tmp, sep = "\t", quote = FALSE, na = "", compress = "none")
  file.rename(tmp, dst)
  rm(dt, key); gc(FALSE)
}

dedup_tsv(
  file.path(liuxy, "chplot/combind_plasma_B/genotype/usefile/snp_annotation.txt"),
  file.path(snp_dir, "pbmc_panel.txt"), "chromosome", "Pos"
)
dedup_tsv(
  file.path(liuxy, "test_peer_glom/exp0.2/predictDB/usefile/new_gt/snp_annotation.txt"),
  file.path(snp_dir, "kidney_panel.txt"), "chromosome", "Pos"
)

chr_patterns <- list(
  pbmc = file.path(liuxy, "MediCell/pbmc/predictdb/usefile/split_chr/snp_annot/snp_annot.chr%d.txt"),
  glom = file.path(liuxy, "test_peer_glom/exp0.2/predictDB/usefile/new_gt/snp_annot/snp_annot.chr%d.txt"),
  tub  = file.path(liuxy, "MediCell/kidney/tub/predictdb_2/usefile/snp_annot/split_chr/snp_annot.chr%d.txt")
)
for (tis in names(chr_patterns)) {
  for (chr in 1:22) {
    dedup_tsv(sprintf(chr_patterns[[tis]], chr),
              file.path(snp_dir, tis, sprintf("snp_annot.chr%d.txt", chr)),
              "chromosome", "pos")
  }
}

# C. Clean multiallelic/duplicate autosome loci from GWAS tables & emit manifests.
autosomes <- as.character(1:22)
bases <- c("A", "C", "G", "T")

make_spec <- function(study_id, path, effect_scale, chr_col = "chromosome", pos_col = "base_pair_location",
                      ea_col = "effect_allele", oa_col = "other_allele", effect_mode = "beta_se",
                      effect_col = "beta", se_col = "standard_error", p_col = "p_value",
                      eaf_col = "effect_allele_frequency", n_mode = "none", n_value = "", n_col = "",
                      dedup = TRUE) {
  list(study_id = study_id, path = path, build = "GRCh37", effect_scale = effect_scale,
       chr_col = chr_col, pos_col = pos_col, ea_col = ea_col, oa_col = oa_col,
       effect_mode = effect_mode, effect_col = effect_col, se_col = se_col, p_col = p_col,
       eaf_col = eaf_col, n_mode = n_mode, n_value = n_value, n_col = n_col, dedup = dedup)
}

prep_gwas <- function(sp) {
  if (!isTRUE(sp$dedup)) return(sp$path)
  out_path <- file.path(gwas_dir, paste0(sp$study_id, "_uniq.tsv.gz"))
  if (file.exists(out_path) && file.info(out_path)$size > 0) return(normalizePath(out_path))
  cols <- unique(unlist(sp[c("chr_col", "pos_col", "ea_col", "oa_col", "effect_col", "se_col", "p_col", "eaf_col", "n_col")],
                        use.names = FALSE))
  cols <- cols[nzchar(cols)]
  dt <- fread(sp$path, sep = "\t", header = TRUE, quote = "", na.strings = NULL, select = cols)
  dt <- dt[as.character(dt[[sp$chr_col]]) %in% autosomes]
  pos_num <- suppressWarnings(as.numeric(dt[[sp$pos_col]]))
  eff_num <- suppressWarnings(as.numeric(dt[[sp$effect_col]]))
  ea <- toupper(as.character(dt[[sp$ea_col]]))
  oa <- toupper(as.character(dt[[sp$oa_col]]))
  positive_eff <- sp$effect_mode != "beta_se"
  valid <- is.finite(pos_num) & pos_num > 0 & pos_num == floor(pos_num) &
           ea %in% bases & oa %in% bases & ea != oa &
           is.finite(eff_num) & (!positive_eff | eff_num > 0)
  dt <- dt[valid]
  key <- paste0(dt[[sp$chr_col]], ":", format(pos_num[valid], scientific = FALSE, trim = TRUE))
  dt <- dt[!(duplicated(key) | duplicated(key, fromLast = TRUE))]
  tmp <- paste0(out_path, ".tmp.gz")
  fwrite(dt, tmp, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
  file.rename(tmp, out_path)
  rm(dt, key, pos_num, eff_num, ea, oa, valid); gc(FALSE)
  normalizePath(out_path)
}

stage_manifest <- function(name, specs) {
  fields <- c("study_id", "path", "build", "effect_scale", "chr_col", "pos_col",
              "ea_col", "oa_col", "effect_mode", "effect_col", "se_col", "p_col",
              "eaf_col", "n_mode", "n_value", "n_col")
  rows <- lapply(specs, function(sp) {
    sp$path <- prep_gwas(sp)
    sp[fields]
  })
  df <- rbindlist(rows)[, ..fields]
  man_path <- file.path(man_dir, paste0(name, ".tsv"))
  tmp <- paste0(man_path, ".tmp")
  fwrite(df, tmp, sep = "\t", quote = FALSE, na = "", compress = "none")
  file.rename(tmp, man_path)
}

# 1) SLE: 3-study GRCh37 meta-analysis (GCST011096 already has zero duplicate positions)
stage_manifest("SLE", list(
  make_spec("GCST011096", file.path(liuxy, "check/other_SLE_GWAS/GCST011096/GCST011096_buildGRCh37.tsv"),
            "log_odds", n_mode = "constant", n_value = "25268", dedup = FALSE),
  make_spec("GCST90018917", file.path(liuxy, "check/other_SLE_GWAS/GCST90018917/GCST90018917_buildGRCh37.tsv.gz"),
            "log_odds", n_mode = "constant", n_value = "659165"),
  make_spec("GCST90476183", file.path(liuxy, "check/other_SLE_GWAS/GCST90476183/GCST90476183_hg19.tsv.gz"),
            "log_odds", chr_col = "chr_hg19", pos_col = "pos_hg19",
            effect_mode = "or_p", effect_col = "odds_ratio", se_col = "",
            n_mode = "column", n_col = "n")
))

# 2) eGFR: Karczewski GCST90691925 (raw file uses neg_log_10_p_value; empty p_col causes
#    gwas_normalize.R to compute two-sided Wald p-values directly from beta/standard_error)
stage_manifest("eGFR", list(
  make_spec("GCST90691925", file.path(liuxy, "add_plot/allimmu_TWAS/new/usefile/eGFR_Karczewski/GCST90691925.tsv.gz"),
            "linear", p_col = "")
))

# 3) 7 UKB metabolic/renal traits (GRCh37 hg19 lifted coordinates)
for (tr in c("Glu", "UA", "uAlb", "uCr", "uK", "uNa", "SBP")) {
  stage_manifest(tr, list(
    make_spec(tr, file.path(liuxy, sprintf("MediCell/trait_compare/metabolic/usefile/%s/%s_hg19.tsv.gz", tr, tr)),
              "linear", chr_col = "chr_hg19", pos_col = "pos_hg19",
              n_mode = "column", n_col = "n")
  ))
}
RS
    log "Input staging complete."
}

# ------------------------------------------------------------------------------
# 5. Tissue-specific cohorts and model configurations
# ------------------------------------------------------------------------------
set_tissue_config() {
    local tissue="$1"
    TISSUE_CAT_COVS=""
    TISSUE_SAMPLE_POLICY="exact"
    TISSUE_GENE_ANNOT_FMT="standard"
    case "$tissue" in
        PBMC)
            # OneK1K PBMC (971 samples, 6 immune cell types)
            # legacy-rnames aligns X-prefixed CIBERSORTx IDs (X1_1) with bulk/genotype IDs (1_1).
            # legacy-name-is-id supplies gene_id from gene_name (absent in PBMC gene_annot.parsed.txt).
            # Manuscript PBMC SLE used unfiltered models; PBMC EUR traits used filtered-signif.
            TISSUE_CELL_TYPES="B,CD4T,CD8T,NK,Mono,DC"
            TISSUE_PEER_NK=2
            TISSUE_PEER_KEEP=2
            TISSUE_SLE_MODEL_SET="unfiltered"
            TISSUE_SAMPLE_POLICY="legacy-rnames"
            TISSUE_GENE_ANNOT_FMT="legacy-name-is-id"
            TISSUE_BULK_TPM="$LIUXY_ROOT/chplot/combind_plasma_B/enigma/usefile/pbmc_6_CellType_971samples_pseudobulk_TPM.txt"
            TISSUE_ENIGMA_REF="$LIUXY_ROOT/chplot/combind_plasma_B/enigma/usefile/TPM_ENIGMA_ref_6CT.matrx"
            TISSUE_CBS_FRAX="$ZHAOQ_ROOT/scTWAS_PBMC/ENIGMA/new_result_54genes/CIBERSORTx_Job129_Adjusted_pbmc_54genes.txt"
            TISSUE_COVARIATES="$LIUXY_ROOT/sc_eQTL_science2022/OneK1K_sex_age_10PCs.txt"
            TISSUE_GENE_ANNOT="$LIUXY_ROOT/MediCell/pbmc/predictdb/usefile/gene/gene_annot.parsed.txt"
            TISSUE_GT_PATTERN="$LIUXY_ROOT/MediCell/pbmc/predictdb/usefile/split_chr/genotype/genotype.chr{chr}.txt"
            TISSUE_SNP_ANNOT_PATTERN="$PREP_DIR/snp_annot/pbmc/snp_annot.chr{chr}.txt"
            TISSUE_SNP_ANNOT_FILE="$PREP_DIR/snp_annot/pbmc_panel.txt"
            ;;
        Glom)
            # Human kidney glomerulus (303 samples, 9 cell types, 60 PEER factors / 15 kept)
            TISSUE_CELL_TYPES="Bcell,Endo,MC,MNP,NK,NKT,Podo,Th_cell,CD8T"
            TISSUE_PEER_NK=60
            TISSUE_PEER_KEEP=15
            TISSUE_SLE_MODEL_SET="filtered-signif"
            TISSUE_CAT_COVS="Pairend,Gender,Batch"
            TISSUE_BULK_TPM="$LIUXY_ROOT/chplot/kidney_new/usefile/glom.bulk.EXP.tpm.Name.Unique.ProtCoding.txt"
            TISSUE_ENIGMA_REF="$LIUXY_ROOT/chplot/kidney_new/usefile/TPM_ENIGMA_ref_withoutNeu_Mast_pDC_ProtCodinng.matrix"
            TISSUE_CBS_FRAX="$ZHAOQ_ROOT/scTWAS_Kidney/ENIGMA/max_L2_TPM/glom/input/CIBERSORTx_Job65_Adjusted_quan.txt"
            TISSUE_COVARIATES="$PREP_DIR/covariates/glom_covariates.tsv"
            TISSUE_GENE_ANNOT="$LIUXY_ROOT/MediCell/kidney/glom/predictdb/usefile/gene_annot.parsed.txt"
            TISSUE_GT_PATTERN="$LIUXY_ROOT/test_peer_glom/exp0.2/predictDB/usefile/new_gt/geno/genotype.chr{chr}.txt"
            TISSUE_SNP_ANNOT_PATTERN="$PREP_DIR/snp_annot/glom/snp_annot.chr{chr}.txt"
            TISSUE_SNP_ANNOT_FILE="$PREP_DIR/snp_annot/kidney_panel.txt"
            ;;
        Tub)
            # Human kidney tubule (356 samples, 12 cell types, 60 PEER factors / 11 kept)
            TISSUE_CELL_TYPES="Bcell,CD8T,CDICA,CDICB,DCT,LOH,MNP,NK,NKT,PT,Th_cell,Plasmacytoid_DC"
            TISSUE_PEER_NK=60
            TISSUE_PEER_KEEP=11
            TISSUE_SLE_MODEL_SET="filtered-signif"
            TISSUE_CAT_COVS="Pairend,Gender,Batch"
            TISSUE_BULK_TPM="$LIUXY_ROOT/chplot/kidney_new/usefile/tub.bulk.EXP.tpm.Name.Unique.ProtCoding.txt"
            TISSUE_ENIGMA_REF="$LIUXY_ROOT/chplot/kidney_new/usefile/TPM_ENIGMA_ref_withoutNeu_Mast_ProtCoding.matrx"
            TISSUE_CBS_FRAX="$ZHAOQ_ROOT/scTWAS_Kidney/ENIGMA/max_L2_TPM/tub/input/pre/CIBERSORTx_Job67_Adjusted_quan.txt"
            TISSUE_COVARIATES="$PREP_DIR/covariates/tub_covariates.tsv"
            TISSUE_GENE_ANNOT="$LIUXY_ROOT/MediCell/kidney/glom/predictdb/usefile/gene_annot.parsed.txt"
            TISSUE_GT_PATTERN="$LIUXY_ROOT/MediCell/kidney/tub/predictdb_2/usefile/genotype/split_chr/genotype.chr{chr}.txt"
            TISSUE_SNP_ANNOT_PATTERN="$PREP_DIR/snp_annot/tub/snp_annot.chr{chr}.txt"
            TISSUE_SNP_ANNOT_FILE="$PREP_DIR/snp_annot/kidney_panel.txt"
            ;;
        *)
            echo "Unsupported tissue: $tissue" >&2
            exit 1
            ;;
    esac
}

# ------------------------------------------------------------------------------
# 6. Stage 2: Cell-type TWAS per tissue (SLE -> eGFR -> 7 UKB traits)
# ------------------------------------------------------------------------------
run_tissue_twas() {
    local tissue="$1"
    set_tissue_config "$tissue"

    local sle_dir="$RUNS_DIR/$tissue/SLE"
    local egfr_dir="$RUNS_DIR/$tissue/eGFR"

    local -a common_args=(
        --cores "$CORES"
        --cell-types "$TISSUE_CELL_TYPES"
        --training-population EUR
        --genome-build GRCh37
        --snp-annot-file "$TISSUE_SNP_ANNOT_FILE"
        --gene-pos "$GENE_POS"
    )
    local -a eur_common_args=(
        --twas-model-set filtered-signif
        --skip-gwas-meta
        --gwas-population EUR
        --ld-panel-id EUR
    )
    local -a extra_args=()
    if [[ -n "$TISSUE_CAT_COVS" ]]; then
        extra_args+=(--cat-covs "$TISSUE_CAT_COVS")
    fi

    # 2a. SLE: full 12-step pipeline (builds M01-M06 models once) + EUR_EAS meta-analysis
    #     (`--gwas-filter-policy legacy` filters discordant METAL directions and HetISq>50 & HetPVal<0.05)
    if [[ ! -s "$sle_dir/M08_summary/allcelltype_result.csv" ]]; then
        log "Running TWAS: $tissue | SLE (full pipeline, EUR_EAS meta)"
        run_cmd python3 -u medicell.py --cell-type-twas \
            --output-dir "$sle_dir" \
            "${common_args[@]}" \
            "${SHARED_MODEL_ARGS[@]}" \
            --peer-nk "$TISSUE_PEER_NK" \
            --peer-keep "$TISSUE_PEER_KEEP" \
            --sample-id-policy "$TISSUE_SAMPLE_POLICY" \
            --gene-annot-format "$TISSUE_GENE_ANNOT_FMT" \
            --twas-model-set "$TISSUE_SLE_MODEL_SET" \
            --bulk-tpm "$TISSUE_BULK_TPM" \
            --enigma-ref "$TISSUE_ENIGMA_REF" \
            --cibersortx-frax "$TISSUE_CBS_FRAX" \
            --covariates "$TISSUE_COVARIATES" \
            --snp-annot-pattern "$TISSUE_SNP_ANNOT_PATTERN" \
            --genotype-pattern "$TISSUE_GT_PATTERN" \
            --gene-annot "$TISSUE_GENE_ANNOT" \
            --gwas-manifest "$PREP_DIR/manifests/SLE.tsv" \
            --gwas-id SLE --gwas-filter-policy legacy \
            --gwas-population EUR_EAS --ld-panel-id EUR_EAS \
            --g1000-vcf-pattern "$G1000_VCF_PATTERN" \
            --g1000-samples "$SAMPLES_EUR_EAS" \
            --ld-dosage-dir "$LD_DIR/${tissue}_EUR_EAS" \
            ${extra_args[@]+"${extra_args[@]}"}
    else
        log "Skipping TWAS: $tissue | SLE (already complete)"
    fi

    # 2b. eGFR: reuse tissue expression/models from SLE, extract/build EUR LD once
    #     (`--steps` is mutually exclusive with `--cell-type-twas`)
    if [[ ! -s "$egfr_dir/M08_summary/allcelltype_result.csv" ]]; then
        log "Running TWAS: $tissue | eGFR (reuse models, build EUR LD)"
        link_shared_stages "$sle_dir" "$egfr_dir" 0
        run_cmd python3 -u medicell.py \
            --steps alignGWAS extractLD buildLD runTWAS summarizeTWAS \
            --output-dir "$egfr_dir" \
            "${common_args[@]}" \
            "${eur_common_args[@]}" \
            --snp-annot-pattern "$TISSUE_SNP_ANNOT_PATTERN" \
            --gwas-manifest "$PREP_DIR/manifests/eGFR.tsv" \
            --gwas-id eGFR \
            --g1000-vcf-pattern "$G1000_VCF_PATTERN" \
            --g1000-samples "$SAMPLES_EUR" \
            --ld-dosage-dir "$LD_DIR/${tissue}_EUR"
    else
        log "Skipping TWAS: $tissue | eGFR (already complete)"
    fi

    # 2c. 7 UKB traits: reuse both tissue models and EUR LD covariances built on eGFR
    local trait
    for trait in "${UKB_TRAITS[@]}"; do
        local tr_dir="$RUNS_DIR/$tissue/$trait"
        if [[ -s "$tr_dir/M08_summary/allcelltype_result.csv" ]]; then
            log "Skipping TWAS: $tissue | $trait (already complete)"
            continue
        fi
        log "Running TWAS: $tissue | $trait (reuse models + EUR LD)"
        link_shared_stages "$egfr_dir" "$tr_dir" 1
        run_cmd python3 -u medicell.py \
            --steps alignGWAS runTWAS summarizeTWAS \
            --output-dir "$tr_dir" \
            "${common_args[@]}" \
            "${eur_common_args[@]}" \
            --gwas-manifest "$PREP_DIR/manifests/${trait}.tsv" \
            --gwas-id "$trait"
    done
}

# ------------------------------------------------------------------------------
# 7. Stage 3: LINCS drug prioritization (SLE pos & eGFR neg)
# ------------------------------------------------------------------------------
# Reads M01_deconvolution and M08_summary/allcelltype_result.csv from each TWAS
# run dir; `--lincs-nes-filter neg` prioritizes compounds reversing disease TWAS
# direction (up-regulated risk genes in SLE [`pos`]; reduced eGFR genes [`neg`]).
run_lincs() {
    local tissue="$1" trait="$2" direction="$3"
    local run_dir="$RUNS_DIR/$tissue/$trait"
    local final_csv="$run_dir/L03_drug_scores/${tissue}_${trait}_${direction}.drug_celltype_score.csv"
    if [[ -s "$final_csv" ]]; then
        log "Skipping LINCS: $tissue | $trait ($direction) (already complete)"
        return 0
    fi
    log "Running LINCS: $tissue | $trait (direction=$direction)"
    run_cmd python3 -u medicell.py lincs \
        --output-dir "$run_dir" \
        --cores "$CORES" \
        --trait "$trait" \
        --tissue "$tissue" \
        --direction "$direction" \
        --gene-pos "$GENE_POS" \
        "${SHARED_LINCS_ARGS[@]}"
}

# ------------------------------------------------------------------------------
# 8. Entrypoint
# ------------------------------------------------------------------------------
main() {
    local tissue

    if [[ "$STAGE" == "all" || "$STAGE" == "prepare" ]]; then
        prepare_inputs
    fi

    if [[ "$STAGE" == "all" || "$STAGE" == "twas" ]]; then
        for tissue in "${SELECTED_TISSUES[@]}"; do
            run_tissue_twas "$tissue"
        done
    fi

    if [[ "$STAGE" == "all" || "$STAGE" == "lincs" ]]; then
        for tissue in "${SELECTED_TISSUES[@]}"; do
            run_lincs "$tissue" SLE pos
            run_lincs "$tissue" eGFR neg
        done
    fi

    log "Finished stage=$STAGE (tissues=${SELECTED_TISSUES[*]})."
}

main
