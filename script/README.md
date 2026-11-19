# MediCell Scripts (`script/`)

This directory holds the pipeline step entrypoints executed by `medicell.py`
and the helper modules they invoke under `script/helpers/`. Entrypoint scripts
receive configuration via `MC_*` environment variables (and positional arguments
where noted) rather than hardcoded local paths.

## Directory Layout

- `script/`: top-level step entrypoints (`M01_*`–`M08_*` for cell-type TWAS, `L01_*`–`L03_*` for
  LINCS therapeutic prioritization).
- `script/helpers/`: modular helper libraries, validation contracts, parallel
  runners, and sub-step utilities called by the entrypoint scripts.

---

## Cell-Type TWAS Pipeline Entrypoints

Invoked by `python3 medicell.py` (using `--cell-type-twas`, `--steps`, `--<step>`, or
positional `<step>`).

| Step (`--<step>`) | Script | Conda Env | Parallelism (`--cores`) | Purpose |
|---|---|---|---|---|
| `inferExp` | `M01_ENIGMA.R` | `R4.4.3` | BLAS threads (`OMP`/`OPENBLAS`) | Run `ENIGMA_trace_norm` deconvolution on bulk TPM using single-cell reference profiles and CIBERSORTx fractions; save `enigma.cse.test`. |
| `normExp` | `M02_DataProcess_exp_fordb.R` | `r3.6.3` | `mclapply` across cell types | Filter genes (`TPM > MC_TPM_MIN` in `>= max(ceil(MC_MIN_PCT * N), MC_MIN_SAMPLES)` samples) and apply per-gene rank-based inverse normal transformation (INT). |
| `peerCov` | `M02_DataProcess_covariance_for_db.R` | `r3.6.3` | `mclapply` across cell types | Auto-detect covariate layout, dummy-encode `MC_CAT_COVS`, estimate `MC_PEER_NK` PEER hidden factors per cell type, and retain `MC_PEER_KEEP` factors for training. |
| `trainModel` | `M03_ModelTraining.R` | `R4.4.3` | Outer `(ct x chr)` pool + inner `doParallel` | Train nested-CV elastic-net GReX models (`alpha = 0.5`, `+/- 1 Mb` cis window) for one `(cell_type, chr)` pair via `helpers/gtex_v7_nested_cv_elnet.R`. |
| `buildDB` | `M04_Makedb.R` | `R4.4.3` | `mclapply` across cell types | Validate 24-column model summaries and canonical `VarID` weights, then assemble `{ct}_models_unfiltered.db` (retained for `buildLD`) and gene-filtered `{ct}_models_filtered_signif.db` (`rho_avg_squared > 0.01`, `zscore_pval < 0.05`). |
| `alignGWAS` | `M05_DataProcess_GWAS_align4METAL.R` | `R4.4.3` | Single-process | Dispatch GWAS alignment to `helpers/gwas_normalize.R` (standard manifest mode) or `helpers/legacy_sle_align.R` (explicit historical SLE profiles). |
| `metaGWAS` *(opt.)* | `M05_DataProcess_GWAS_metal.sh` | `R4.4.3` | Single-process | Dispatch METAL inverse-variance meta-analysis via `helpers/gwas_control.py meta`. |
| `filterGWAS` *(opt.)* | `M05_DataProcess_GWAS_metal_filter.R` | `R4.4.3` | `data.table` | Normalize METAL alleles and optionally apply legacy heterogeneity/direction filtering (`MC_GWAS_FILTER_POLICY=legacy` vs `none`). |
| `extractLD` | `M06_LD.sh` | `imlatools` | `parallel_run.py` across chromosomes | Extract biallelic SNPs from 1000G VCFs with pre-subset AF (`bcftools +fill-tags`) and align dosages to per-chromosome SNP annotations. |
| `buildLD` | `M06_Makeld.sh` | `imlatools` | `parallel_run.py` across cell types | Run MetaXcan `M01_covariances_correlations.py` per cell type on isolated chromosome dosages and validate outputs. |
| `runTWAS` | `M07_TWAS.sh` | `imlatools` | `parallel_run.py` across cell types | Run `SPrediXcan.py` per cell type against the selected weight DB (`unfiltered` or `filtered-signif`) and validate TWAS output schemas. |
| `summarizeTWAS` | `M08_ResultStatistic.R` | `R4.4.3` | `mclapply` across cell types | Map `gene_id` to `gene_name`, compute per-cell-type BH FDR and Bonferroni adjusted p-values, and write the 4 summary CSVs. |

---

## LINCS Therapeutic Prioritization Entrypoints

Invoked by `python3 medicell.py lincs` (or `python3 medicell.py --lincs`).

| Step (`--<step>`) | Script | Conda Env | Parallelism (`--cores`) | Purpose |
|---|---|---|---|---|
| `prepLincs` | `L01_lincs_prepare.R` | auto (`R4.4.3` / `r_4.6.1`) | Single-process | Read TWAS results (`.csv`/`.tsv`/`.xlsx`), split by z-score sign (`pos`/`neg`), write normalized `Gene/CellType/Zscore` table, and prepare gene-symbol x cell-type expression matrix (from `--exp-matrix` or `inferExp` ENIGMA `.RData`). |
| `enrichLincs` | `L02_lincs_enrichment.R` | auto (`R4.4.3` / `r_4.6.1`) | `fgsea(..., nproc)` (+ parallel `pos`/`neg` with `--direction both`) | Run `fgsea` enrichment on directional ranked TWAS genes against a GMT library (default `data/ref/lincs/LINCS_PRL_top250_human.gmt`), attach optional L1000 Phase 1/2 perturbation metadata, and filter hits. |
| `scoreLincs` | `L03_lincs_drug_score.sh` | auto (`R4.4.3` / `r_4.6.1`) | Parallel `pos`/`neg` with `--direction both` | Orchestrate the 4 LINCS scoring helpers in `script/helpers/` to filter top perturbations, extract leading-edge TWAS z-scores, and write drug x cell-type `RankScore` matrices (`.RData` and `.csv`). |

---

## Helper Modules (`script/helpers/`)

| File | Caller(s) | Purpose |
|---|---|---|
| `input_contracts.R` | `M01_ENIGMA.R`, `M02_DataProcess_exp_fordb.R`, `M02_DataProcess_covariance_for_db.R`, `gtex_v7_nested_cv_elnet.R` | Input matrix/ID/numeric validators (`mc_matrix`, `mc_celltypes`, `mc_sample_match`, `mc_number`, `mc_choice`, `mc_load_enigma_inputs`). |
| `gtex_v7_nested_cv_elnet.R` | `M03_ModelTraining.R` | Nested cross-validation elastic-net (`glmnet`) training core with TSS-window validation/collapse, `[0, 2]` dosage and `MC_MAF`/`MC_MAF_POLICY` enforcement, writing 24-column summaries, weights, and SNP covariances. |
| `gwas_normalize.R` | `M05_DataProcess_GWAS_align4METAL.R` | Manifest-driven standard GWAS normalization (strict build check, `beta_se`/`or_se`/`or_p` handling, no silent strand guessing). |
| `gwas_control.py` | `M05_DataProcess_GWAS_metal.sh` | Generates METAL control scripts, runs `metal`, and validates output variant counts/schema for standard and legacy modes. |
| `legacy_sle_align.R` | `M05_DataProcess_GWAS_align4METAL.R` | Isolated legacy adapter reproducing historical SLE GWAS alignment (retained solely for `--gwas-profile` reproduction). |
| `legacy_sle_metal.sh` | `gwas_control.py` | Isolated legacy adapter reproducing historical SLE METAL formatting (retained solely for `--gwas-profile` reproduction). |
| `ld_contracts.py` | `M06_LD.sh`, `M06_Makeld.sh` | Extracts SNP site lists, validates sample order, converts `bcftools` genotypes to dosages, aligns dosages to SNP annotations, and checks LD outputs. |
| `parallel_run.py` | `M06_LD.sh`, `M06_Makeld.sh`, `M07_TWAS.sh` | Portable Python `ThreadPoolExecutor` runner for executing exported Bash worker functions in parallel with fail-fast error propagation. |
| `lincs_fgsea_lib.R` | `L02_lincs_enrichment.R` | GMT parser, z-score direction inference, `fgsea` wrapper, and L1000 Phase 1/2 metadata annotator. |
| `lincsTissueSpecific.R` | `L03_lincs_drug_score.sh` | Filters LINCS hits (e.g., excluding `^BRD-\|trt_sh.css`), ranks by `pval` or `padj`, and exports top-N perturbations and `leadingEdge` gene lists. |
| `formatLeadingSpecific.pl` | `L03_lincs_drug_score.sh` | Splits `pert_iname` / `leadingEdge` tables into a deduplicated gene pool and `drug -> gene` pairs. |
| `formatZscoreSpecific.pl` | `L03_lincs_drug_score.sh` | Extracts TWAS z-scores for leading-edge genes per cell type (supports named CSV/TSV headers as well as legacy 13-column whitespace tables). |
| `getZscoreSpecificMatrix.R` | `L03_lincs_drug_score.sh` | Computes cross-cell-type expression specificity (`expZscore`, requiring `>= 2` overlapping cell types) and sums `|TWAS_zscore| * expZscore` over leading-edge genes for each drug. |

