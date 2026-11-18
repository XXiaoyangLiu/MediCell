#!/usr/bin/env bash
# End-to-end MediCell runner for manuscript TWAS and LINCS analyses.
# Scope: ENIGMA/PEER/PredictDB/GWAS/LD/TWAS across 3 tissues
# (PBMC, Glom, Tub) x 9 GRCh37 traits (SLE, eGFR, Glu, UA, uAlb, uCr, uK, uNa, SBP),
# followed by LINCS prioritization (lincs) for SLE (pos) and eGFR (neg).
# No plotting or non-MediCell downstream analyses are included.
set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CORES="${MC_CORES:-16}"
OUT_ROOT="$ROOT/Results/out"
DRY_RUN=0

usage() {
    cat <<USAGE
Usage: bash Results/run.sh [--cores N] [--out-dir DIR] [--dry-run] [-h|--help]

Runs all manuscript MediCell TWAS and LINCS workflows:
  Tissues : PBMC, Glom, Tub
  Traits  : SLE (EUR+EAS meta), eGFR (EUR), Glu, UA, uAlb, uCr, uK, uNa, SBP (EUR)
  LINCS   : SLE (direction=pos) and eGFR (direction=neg) across PBMC, Glom, Tub
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cores) CORES="${2:?--cores requires a value}"; shift 2 ;;
        --out-dir) OUT_ROOT="${2:?--out-dir requires a value}"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

case "$OUT_ROOT" in
    /*) ;;
    *) OUT_ROOT="$ROOT/$OUT_ROOT" ;;
esac

PREP_DIR="$OUT_ROOT/prepared"
RUNS_DIR="$OUT_ROOT/runs"
LD_DIR="$OUT_ROOT/ld_dosage"

# Manuscript raw paths
G1000_VCF_PATTERN="/node1/liuxy/index/1000G_2/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
SAMPLES_EUR_EAS="/node1/liuxy/index/1000G/EUR_EAS_samples.txt"
SAMPLES_EUR="/node1/liuxy/index/1000G/EUR_samples.txt"
GENE_POS="/node1/liuxy/index/gene_pos_grch37.82.csv"

UKB_TRAITS=(Glu UA uAlb uCr uK uNa SBP)
TISSUES=(PBMC Glom Tub)

find_r4_bin() {
    if [[ -x "/public/home/dongb/miniforge3/envs/R4.4.3/bin/Rscript" ]]; then
        echo "/public/home/dongb/miniforge3/envs/R4.4.3/bin/Rscript"
    elif command -v Rscript >/dev/null 2>&1; then
        command -v Rscript
    else
        echo "Cannot locate R4.4.3 Rscript binary" >&2
        exit 1
    fi
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

# Stage clean inputs compatible with MediCell strict contracts:
# 1) Export Glom/Tub covariate tables from QTLData.RData
# 2) Filter SNP panels and per-chr annotations to unambiguous biallelic (chr, pos) loci
# 3) Filter autosome biallelic GWAS inputs to unique (chr, pos) loci and emit manifests
prepare_inputs() {
    mkdir -p "$PREP_DIR/covariates" "$PREP_DIR/snp_annot/pbmc" "$PREP_DIR/snp_annot/glom" "$PREP_DIR/snp_annot/tub" \
             "$PREP_DIR/gwas" "$PREP_DIR/manifests"

    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "[DRY-RUN] prepare_inputs -> $PREP_DIR"
        return 0
    fi

    local r4_bin
    r4_bin="$(find_r4_bin)"

    # 1. Glom & Tub covariates from QTLData.RData (rownames are sample IDs matching bulk/CIBERSORTx/genotypes)
    "$r4_bin" - "$PREP_DIR/covariates" <<'RS'
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
cov_dir <- args[1]

export_qtl_cov <- function(rdata_path, obj_name, cols, out_path) {
  if (file.exists(out_path)) return(invisible(NULL))
  e <- new.env(parent = emptyenv())
  load(rdata_path, envir = e)
  df <- get(obj_name, envir = e)
  out <- data.frame(Sample = rownames(df), df[, cols, drop = FALSE],
                    stringsAsFactors = FALSE, check.names = FALSE)
  write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

export_qtl_cov(
  "/node1/zhaoq/scTWAS_Kidney/predictDB/kidney_KS/glom/QTLData.RData",
  "AAglomPhe2",
  c("Age", "RIN", "PC1", "PC2", "PC3", "PC4", "PC5", "Sclerosis", "Pairend", "Gender", "Batch"),
  file.path(cov_dir, "glom_covariates.tsv")
)
export_qtl_cov(
  "/node1/zhaoq/scTWAS_Kidney/predictDB/kidney_KS/tub/QTLData.RData",
  "AAtubPhe2",
  c("Age", "RIN", "PC1", "PC2", "PC3", "PC4", "PC5", "Fibrosis", "Pairend", "Gender", "Batch"),
  file.path(cov_dir, "tub_covariates.tsv")
)
RS

    # 2. Deduplicate biallelic (chr, pos) loci across SNP panels and chr1..22 annotations
    python3 - "$PREP_DIR/snp_annot" <<'PY'
import csv
from collections import Counter
from pathlib import Path
import sys

out_root = Path(sys.argv[1])

def dedup_tsv(src: Path, dst: Path, chr_col: str, pos_col: str) -> None:
    if dst.exists() and dst.stat().st_size > 0:
        return
    counts = Counter()
    with src.open("r", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            counts[(row[chr_col], row[pos_col])] += 1
    tmp = dst.with_suffix(dst.suffix + ".tmp")
    with src.open("r", newline="") as fin, tmp.open("w", newline="") as fout:
        r = csv.DictReader(fin, delimiter="\t")
        w = csv.DictWriter(fout, fieldnames=r.fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in r:
            if counts[(row[chr_col], row[pos_col])] == 1:
                w.writerow(row)
    tmp.replace(dst)

dedup_tsv(
    Path("/node1/liuxy/chplot/combind_plasma_B/genotype/usefile/snp_annotation.txt"),
    out_root / "pbmc_panel.txt", "chromosome", "Pos",
)
dedup_tsv(
    Path("/node1/liuxy/test_peer_glom/exp0.2/predictDB/usefile/new_gt/snp_annotation.txt"),
    out_root / "kidney_panel.txt", "chromosome", "Pos",
)

chr_patterns = {
    "pbmc": "/node1/liuxy/MediCell/pbmc/predictdb/usefile/split_chr/snp_annot/snp_annot.chr{chr}.txt",
    "glom": "/node1/liuxy/test_peer_glom/exp0.2/predictDB/usefile/new_gt/snp_annot/snp_annot.chr{chr}.txt",
    "tub": "/node1/liuxy/MediCell/kidney/tub/predictdb_2/usefile/snp_annot/split_chr/snp_annot.chr{chr}.txt",
}
for tissue, pat in chr_patterns.items():
    for c in range(1, 23):
        dedup_tsv(Path(pat.format(chr=c)), out_root / tissue / f"snp_annot.chr{c}.txt", "chromosome", "pos")
PY

    # 3. Clean duplicate/multiallelic autosome loci from GWAS inputs & emit manifest TSVs
    MC_PREP_CORES="$CORES" "$r4_bin" - "$PREP_DIR/gwas" "$PREP_DIR/manifests" <<'RS'
options(scipen = 999)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
gwas_dir <- args[1]
man_dir <- args[2]
cores <- suppressWarnings(as.integer(Sys.getenv("MC_PREP_CORES", "8")))
if (is.na(cores) || cores < 1L) cores <- 8L
setDTthreads(min(cores, 16L))

autosomes <- as.character(1:22)
bases <- c("A", "C", "G", "T")
prep_gwas <- function(raw_path, out_path, chr_col, pos_col, ea_col, oa_col, eff_col, keep_cols, positive_eff = FALSE) {
  if (file.exists(out_path) && file.info(out_path)$size > 0) return(normalizePath(out_path))
  cols <- unique(keep_cols[nzchar(keep_cols)])
  dt <- fread(raw_path, sep = "\t", header = TRUE, quote = "", na.strings = NULL, select = cols)
  dt <- dt[as.character(dt[[chr_col]]) %in% autosomes]
  pos_num <- suppressWarnings(as.numeric(dt[[pos_col]]))
  eff_num <- suppressWarnings(as.numeric(dt[[eff_col]]))
  ea <- toupper(as.character(dt[[ea_col]]))
  oa <- toupper(as.character(dt[[oa_col]]))
  valid <- is.finite(pos_num) & pos_num > 0 & pos_num == floor(pos_num) &
           ea %in% bases & oa %in% bases & ea != oa &
           is.finite(eff_num) & (!positive_eff | eff_num > 0)
  dt <- dt[valid]
  key <- paste0(dt[[chr_col]], ":", format(pos_num[valid], scientific = FALSE, trim = TRUE))
  dt <- dt[!(duplicated(key) | duplicated(key, fromLast = TRUE))]
  tmp <- paste0(out_path, ".tmp.gz")
  fwrite(dt, tmp, sep = "\t", quote = FALSE, na = "NA")
  file.rename(tmp, out_path)
  rm(dt, key, pos_num, eff_num, ea, oa, valid); gc(FALSE)
  normalizePath(out_path)
}

write_manifest <- function(name, rows) {
  fields <- c("study_id", "path", "build", "effect_scale", "chr_col", "pos_col",
              "ea_col", "oa_col", "effect_mode", "effect_col", "se_col", "p_col",
              "eaf_col", "n_mode", "n_value", "n_col")
  df <- rbindlist(rows)[, ..fields]
  fwrite(df, file.path(man_dir, paste0(name, ".tsv")), sep = "\t", quote = FALSE)
}

# SLE (3-study GRCh37 meta-analysis)
p096 <- "/node1/liuxy/check/other_SLE_GWAS/GCST011096/GCST011096_buildGRCh37.tsv"
p917 <- prep_gwas(
  "/node1/liuxy/check/other_SLE_GWAS/GCST90018917/GCST90018917_buildGRCh37.tsv.gz",
  file.path(gwas_dir, "GCST90018917_uniq.tsv.gz"),
  "chromosome", "base_pair_location", "effect_allele", "other_allele", "beta",
  c("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "p_value", "effect_allele_frequency")
)
p183 <- prep_gwas(
  "/node1/liuxy/check/other_SLE_GWAS/GCST90476183/GCST90476183_hg19.tsv.gz",
  file.path(gwas_dir, "GCST90476183_uniq.tsv.gz"),
  "chr_hg19", "pos_hg19", "effect_allele", "other_allele", "odds_ratio",
  c("chr_hg19", "pos_hg19", "effect_allele", "other_allele", "odds_ratio", "p_value", "effect_allele_frequency", "n"),
  positive_eff = TRUE
)
write_manifest("SLE", list(
  list(study_id="GCST011096", path=p096, build="GRCh37", effect_scale="log_odds",
       chr_col="chromosome", pos_col="base_pair_location", ea_col="effect_allele", oa_col="other_allele",
       effect_mode="beta_se", effect_col="beta", se_col="standard_error", p_col="p_value",
       eaf_col="effect_allele_frequency", n_mode="constant", n_value="25268", n_col=""),
  list(study_id="GCST90018917", path=p917, build="GRCh37", effect_scale="log_odds",
       chr_col="chromosome", pos_col="base_pair_location", ea_col="effect_allele", oa_col="other_allele",
       effect_mode="beta_se", effect_col="beta", se_col="standard_error", p_col="p_value",
       eaf_col="effect_allele_frequency", n_mode="constant", n_value="659165", n_col=""),
  list(study_id="GCST90476183", path=p183, build="GRCh37", effect_scale="log_odds",
       chr_col="chr_hg19", pos_col="pos_hg19", ea_col="effect_allele", oa_col="other_allele",
       effect_mode="or_p", effect_col="odds_ratio", se_col="", p_col="p_value",
       eaf_col="effect_allele_frequency", n_mode="column", n_value="", n_col="n")
))

# eGFR (Karczewski GCST90691925, GRCh37)
p_egfr <- prep_gwas(
  "/node1/liuxy/add_plot/allimmu_TWAS/new/usefile/eGFR_Karczewski/GCST90691925.tsv.gz",
  file.path(gwas_dir, "GCST90691925_uniq.tsv.gz"),
  "chromosome", "base_pair_location", "effect_allele", "other_allele", "beta",
  c("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency")
)
write_manifest("eGFR", list(
  list(study_id="GCST90691925", path=p_egfr, build="GRCh37", effect_scale="linear",
       chr_col="chromosome", pos_col="base_pair_location", ea_col="effect_allele", oa_col="other_allele",
       effect_mode="beta_se", effect_col="beta", se_col="standard_error", p_col="",
       eaf_col="effect_allele_frequency", n_mode="none", n_value="", n_col="")
))

# 7 UKB metabolic/renal traits (GRCh37 hg19 coordinates)
for (tr in c("Glu", "UA", "uAlb", "uCr", "uK", "uNa", "SBP")) {
  raw_tr <- sprintf("/node1/liuxy/MediCell/trait_compare/metabolic/usefile/%s/%s_hg19.tsv.gz", tr, tr)
  p_tr <- prep_gwas(
    raw_tr,
    file.path(gwas_dir, paste0(tr, "_uniq.tsv.gz")),
    "chr_hg19", "pos_hg19", "effect_allele", "other_allele", "beta",
    c("chr_hg19", "pos_hg19", "effect_allele", "other_allele", "beta", "standard_error", "p_value", "effect_allele_frequency", "n")
  )
  write_manifest(tr, list(
    list(study_id=tr, path=p_tr, build="GRCh37", effect_scale="linear",
         chr_col="chr_hg19", pos_col="pos_hg19", ea_col="effect_allele", oa_col="other_allele",
         effect_mode="beta_se", effect_col="beta", se_col="standard_error", p_col="p_value",
         eaf_col="effect_allele_frequency", n_mode="column", n_value="", n_col="n")
  ))
}
RS
}

set_tissue_config() {
    local tissue="$1"
    TISSUE_CAT_COVS=""
    TISSUE_SAMPLE_POLICY="exact"
    TISSUE_GENE_ANNOT_FMT="standard"
    case "$tissue" in
        PBMC)
            TISSUE_CELL_TYPES="B,CD4T,CD8T,NK,Mono,DC"
            TISSUE_PEER_NK=2
            TISSUE_PEER_KEEP=2
            TISSUE_SLE_MODEL_SET="unfiltered"
            TISSUE_SAMPLE_POLICY="legacy-rnames"
            TISSUE_GENE_ANNOT_FMT="legacy-name-is-id"
            TISSUE_BULK_TPM="/node1/liuxy/chplot/combind_plasma_B/enigma/usefile/pbmc_6_CellType_971samples_pseudobulk_TPM.txt"
            TISSUE_ENIGMA_REF="/node1/liuxy/chplot/combind_plasma_B/enigma/usefile/TPM_ENIGMA_ref_6CT.matrx"
            TISSUE_CBS_FRAX="/node1/zhaoq/scTWAS_PBMC/ENIGMA/new_result_54genes/CIBERSORTx_Job129_Adjusted_pbmc_54genes.txt"
            TISSUE_COVARIATES="/node1/liuxy/sc_eQTL_science2022/OneK1K_sex_age_10PCs.txt"
            TISSUE_GENE_ANNOT="/node1/liuxy/MediCell/pbmc/predictdb/usefile/gene/gene_annot.parsed.txt"
            TISSUE_GT_PATTERN="/node1/liuxy/MediCell/pbmc/predictdb/usefile/split_chr/genotype/genotype.chr{chr}.txt"
            TISSUE_SNP_ANNOT_PATTERN="$PREP_DIR/snp_annot/pbmc/snp_annot.chr{chr}.txt"
            TISSUE_SNP_ANNOT_FILE="$PREP_DIR/snp_annot/pbmc_panel.txt"
            ;;
        Glom)
            TISSUE_CELL_TYPES="Bcell,Endo,MC,MNP,NK,NKT,Podo,Th_cell,CD8T"
            TISSUE_PEER_NK=60
            TISSUE_PEER_KEEP=15
            TISSUE_SLE_MODEL_SET="filtered-signif"
            TISSUE_CAT_COVS="Pairend,Gender,Batch"
            TISSUE_BULK_TPM="/node1/liuxy/chplot/kidney_new/usefile/glom.bulk.EXP.tpm.Name.Unique.ProtCoding.txt"
            TISSUE_ENIGMA_REF="/node1/liuxy/chplot/kidney_new/usefile/TPM_ENIGMA_ref_withoutNeu_Mast_pDC_ProtCodinng.matrix"
            TISSUE_CBS_FRAX="/node1/zhaoq/scTWAS_Kidney/ENIGMA/max_L2_TPM/glom/input/CIBERSORTx_Job65_Adjusted_quan.txt"
            TISSUE_COVARIATES="$PREP_DIR/covariates/glom_covariates.tsv"
            TISSUE_GENE_ANNOT="/node1/liuxy/MediCell/kidney/glom/predictdb/usefile/gene_annot.parsed.txt"
            TISSUE_GT_PATTERN="/node1/liuxy/test_peer_glom/exp0.2/predictDB/usefile/new_gt/geno/genotype.chr{chr}.txt"
            TISSUE_SNP_ANNOT_PATTERN="$PREP_DIR/snp_annot/glom/snp_annot.chr{chr}.txt"
            TISSUE_SNP_ANNOT_FILE="$PREP_DIR/snp_annot/kidney_panel.txt"
            ;;
        Tub)
            TISSUE_CELL_TYPES="Bcell,CD8T,CDICA,CDICB,DCT,LOH,MNP,NK,NKT,PT,Th_cell,Plasmacytoid_DC"
            TISSUE_PEER_NK=60
            TISSUE_PEER_KEEP=11
            TISSUE_SLE_MODEL_SET="filtered-signif"
            TISSUE_CAT_COVS="Pairend,Gender,Batch"
            TISSUE_BULK_TPM="/node1/liuxy/chplot/kidney_new/usefile/tub.bulk.EXP.tpm.Name.Unique.ProtCoding.txt"
            TISSUE_ENIGMA_REF="/node1/liuxy/chplot/kidney_new/usefile/TPM_ENIGMA_ref_withoutNeu_Mast_ProtCoding.matrx"
            TISSUE_CBS_FRAX="/node1/zhaoq/scTWAS_Kidney/ENIGMA/max_L2_TPM/tub/input/pre/CIBERSORTx_Job67_Adjusted_quan.txt"
            TISSUE_COVARIATES="$PREP_DIR/covariates/tub_covariates.tsv"
            TISSUE_GENE_ANNOT="/node1/liuxy/MediCell/kidney/glom/predictdb/usefile/gene_annot.parsed.txt"
            TISSUE_GT_PATTERN="/node1/liuxy/MediCell/kidney/tub/predictdb_2/usefile/genotype/split_chr/genotype.chr{chr}.txt"
            TISSUE_SNP_ANNOT_PATTERN="$PREP_DIR/snp_annot/tub/snp_annot.chr{chr}.txt"
            TISSUE_SNP_ANNOT_FILE="$PREP_DIR/snp_annot/kidney_panel.txt"
            ;;
        *) echo "Unsupported tissue: $tissue" >&2; exit 1 ;;
    esac
}

link_shared_stages() {
    local src_dir="$1" dst_dir="$2" include_eur_ld="${3:-0}"
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

run_tissue_twas() {
    local tissue="$1"
    set_tissue_config "$tissue"

    local sle_dir="$RUNS_DIR/$tissue/SLE"
    local egfr_dir="$RUNS_DIR/$tissue/eGFR"

    # 1. SLE: full tissue pipeline + EUR_EAS meta-analysis TWAS
    if [[ ! -s "$sle_dir/M08_summary/allcelltype_result.csv" ]]; then
        local -a extra_args=()
        if [[ -n "$TISSUE_CAT_COVS" ]]; then
            extra_args+=(--cat-covs "$TISSUE_CAT_COVS")
        fi
        run_cmd python3 -u medicell.py --cell-type-twas \
            --output-dir "$sle_dir" \
            --cores "$CORES" \
            --cell-types "$TISSUE_CELL_TYPES" \
            --alpha 0.1 --tpm-min 0.1 --min-pct 0.01 --min-samples 5 \
            --peer-nk "$TISSUE_PEER_NK" --peer-keep "$TISSUE_PEER_KEEP" --min-level-count 5 \
            --sample-id-policy "$TISSUE_SAMPLE_POLICY" \
            --gene-annot-format "$TISSUE_GENE_ANNOT_FMT" \
            --training-population EUR --maf 0.01 --maf-policy filter \
            --twas-model-set "$TISSUE_SLE_MODEL_SET" \
            --genome-build GRCh37 \
            --bulk-tpm "$TISSUE_BULK_TPM" \
            --enigma-ref "$TISSUE_ENIGMA_REF" \
            --cibersortx-frax "$TISSUE_CBS_FRAX" \
            --covariates "$TISSUE_COVARIATES" \
            --snp-annot-pattern "$TISSUE_SNP_ANNOT_PATTERN" \
            --genotype-pattern "$TISSUE_GT_PATTERN" \
            --gene-annot "$TISSUE_GENE_ANNOT" \
            --snp-annot-file "$TISSUE_SNP_ANNOT_FILE" \
            --gwas-manifest "$PREP_DIR/manifests/SLE.tsv" \
            --gwas-id SLE --gwas-filter-policy legacy \
            --gwas-population EUR_EAS --ld-panel-id EUR_EAS \
            --g1000-vcf-pattern "$G1000_VCF_PATTERN" \
            --g1000-samples "$SAMPLES_EUR_EAS" \
            --ld-dosage-dir "$LD_DIR/${tissue}_EUR_EAS" \
            --gene-pos "$GENE_POS" \
            ${extra_args[@]+"${extra_args[@]}"}
    fi

    # 2. eGFR: reuse tissue expression/models, build EUR LD once, run single-study TWAS
    if [[ ! -s "$egfr_dir/M08_summary/allcelltype_result.csv" ]]; then
        if [[ "$DRY_RUN" -eq 0 ]]; then
            link_shared_stages "$sle_dir" "$egfr_dir" 0
        else
            echo "[DRY-RUN] link_shared_stages $sle_dir -> $egfr_dir"
        fi
        run_cmd python3 -u medicell.py \
            --steps alignGWAS extractLD buildLD runTWAS summarizeTWAS \
            --output-dir "$egfr_dir" \
            --cores "$CORES" \
            --cell-types "$TISSUE_CELL_TYPES" \
            --training-population EUR --twas-model-set filtered-signif \
            --genome-build GRCh37 \
            --snp-annot-pattern "$TISSUE_SNP_ANNOT_PATTERN" \
            --snp-annot-file "$TISSUE_SNP_ANNOT_FILE" \
            --gwas-manifest "$PREP_DIR/manifests/eGFR.tsv" \
            --gwas-id eGFR --skip-gwas-meta \
            --gwas-population EUR --ld-panel-id EUR \
            --g1000-vcf-pattern "$G1000_VCF_PATTERN" \
            --g1000-samples "$SAMPLES_EUR" \
            --ld-dosage-dir "$LD_DIR/${tissue}_EUR" \
            --gene-pos "$GENE_POS"
    fi

    # 3. 7 UKB traits: reuse tissue models + EUR LD covariances built for eGFR
    local trait
    for trait in "${UKB_TRAITS[@]}"; do
        local tr_dir="$RUNS_DIR/$tissue/$trait"
        if [[ -s "$tr_dir/M08_summary/allcelltype_result.csv" ]]; then
            continue
        fi
        if [[ "$DRY_RUN" -eq 0 ]]; then
            link_shared_stages "$egfr_dir" "$tr_dir" 1
        else
            echo "[DRY-RUN] link_shared_stages $egfr_dir -> $tr_dir (+cov_EUR)"
        fi
        run_cmd python3 -u medicell.py \
            --steps alignGWAS runTWAS summarizeTWAS \
            --output-dir "$tr_dir" \
            --cores "$CORES" \
            --cell-types "$TISSUE_CELL_TYPES" \
            --training-population EUR --twas-model-set filtered-signif \
            --genome-build GRCh37 \
            --snp-annot-file "$TISSUE_SNP_ANNOT_FILE" \
            --gwas-manifest "$PREP_DIR/manifests/${trait}.tsv" \
            --gwas-id "$trait" --skip-gwas-meta \
            --gwas-population EUR --ld-panel-id EUR \
            --gene-pos "$GENE_POS"
    done
}

run_lincs() {
    local tissue="$1" trait="$2" direction="$3"
    local run_dir="$RUNS_DIR/$tissue/$trait"
    local final_csv="$run_dir/L03_drug_scores/${tissue}_${trait}_${direction}.drug_celltype_score.csv"
    if [[ -s "$final_csv" ]]; then
        return 0
    fi
    run_cmd python3 -u medicell.py lincs \
        --output-dir "$run_dir" \
        --cores "$CORES" \
        --trait "$trait" \
        --tissue "$tissue" \
        --direction "$direction" \
        --gene-pos "$GENE_POS" \
        --twas-padj-cutoff 0.05 \
        --twas-zscore-cutoff 0.25 \
        --lincs-min-size 10 \
        --lincs-max-size 500 \
        --lincs-sig-col pval \
        --lincs-sig-cutoff 0.05 \
        --lincs-nes-filter neg \
        --lincs-rank-by pval \
        --top-n 50
}

main() {
    prepare_inputs
    local tissue
    for tissue in "${TISSUES[@]}"; do
        run_tissue_twas "$tissue"
    done
    for tissue in "${TISSUES[@]}"; do
        run_lincs "$tissue" SLE pos
        run_lincs "$tissue" eGFR neg
    done
}

main
