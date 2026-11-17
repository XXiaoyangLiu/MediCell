# MediCell

MediCell is an integrative TWAS framework that combines single-cell reference
profiles with bulk transcriptomic and genotype data to identify
cell-type-dependent genetic effects. It deconvolves bulk RNA-seq into
cell-type-resolved expression via CIBERSORTx + ENIGMA, trains
cell-type-specific genetically regulated expression (GReX) models with
elastic-net, and performs S-PrediXcan association against GWAS summary
statistics. Part 2 (`medicell.py lincs`) integrates cell-type TWAS associations
and expression profiles with LINCS L1000 perturbation signatures (`fgsea`) to
prioritize candidate drugs per cell type.

Analysis is configured by explicit cell types, PEER factors and model-selection
parameters. No tissue name or tissue-specific presets are required.

## Installation

Requirements: a Linux server with conda (a recent Miniconda/Miniforge; the
installer auto-detects `~/miniforge3/bin/conda` or `$CONDA_EXE`, since system
condas on some clusters are too old) and network access to conda-forge,
bioconda and GitHub.

```bash
git clone git@github.com:delect-aj/MediCell.git
cd MediCell
bash setup.sh
```

`setup.sh` is idempotent (safe to re-run at any time) and installs:

| Component | Description |
|---|---|
| conda env `R4.4.3` | R 4.4.3 with ENIGMA, glmnet, RSQLite and other model-building packages |
| conda env `r3.6.3` | legacy R 3.6.3 with PEER (PEER requires it) |
| conda env `imlatools` | Python MetaXcan runtime + bcftools |
| PEER | R package, built from source (github.com/PMBio/peer) |
| ENIGMA | R package, installed from github.com/WWXkenmo/ENIGMA |
| MetaXcan | cloned to `tools/MetaXcan` |
| METAL | compiled to `tools/METAL/metal` |

To install a lab-internal ENIGMA version instead of the public one:

```bash
bash setup.sh --enigma-src /path/to/ENIGMA
```

## Usage

```bash
python3 medicell.py --output-dir run_out --all \
    --cores 4 \
    --cell-types B,CD4T,CD8T,NK,Mono,DC \
    --peer-nk 2 --peer-keep 2 \
    --training-population EUR --twas-model-set unfiltered \
    --bulk-tpm           data/bulk_TPM.txt \
    --enigma-ref         data/sc_reference_TPM.matrx \
    --cibersortx-frax    data/CIBERSORTx_fractions.txt \
    --covariates         data/covariates.txt \
    --snp-annot-pattern  data/snp_annot.chr{chr}.txt \
    --genotype-pattern   data/genotype.chr{chr}.txt \
    --gene-annot         data/gene_annot.parsed.txt \
    --snp-annot-file     data/snp_annotation.txt \
    --gwas-manifest      data/gwas_manifest.tsv \
    --g1000-vcf-pattern  data/1000G/ALL.chr{chr}.phase3...genotypes.vcf.gz \
    --g1000-samples      data/EUR_EAS_samples.txt \
    --gene-pos           data/gene_pos_grch37.82.csv
```

Every input path pattern takes a `{chr}` placeholder that the driver expands
per chromosome (1–22 by default). The covariates file is required for `peerCov` (`M02cov`, format below).

### Run modes

Use `--all` for the complete pipeline (or omit step selectors). To run a
single step or a subset, use the descriptive step flags (`--inferExp`,
`--summarizeTWAS`, etc.), positional subcommand (`medicell.py inferExp ...`),
or `--steps`:

```bash
# Run only ENIGMA expression inference (M01)
python3 medicell.py --inferExp --output-dir run_out --cell-types B,CD4T \
    --bulk-tpm data/bulk_TPM.txt --enigma-ref data/sc_reference_TPM.matrx \
    --cibersortx-frax data/CIBERSORTx_fractions.txt

# Run only the TWAS summary step (M08), using existing TWAS results in run_out
python3 medicell.py --summarizeTWAS --output-dir run_out --cell-types B \
    --gene-pos data/gene_pos_grch37.82.csv

# Combine multiple step flags (or use --steps normExp,peerCov)
python3 medicell.py --normExp --peerCov --output-dir run_out --cell-types B \
    --peer-nk 2 --peer-keep 2 --covariates data/covariates.txt
```

`--all`, `--steps`, and `--<step>` flags are mutually exclusive (combining
multiple `--<step>` flags is supported). Selected steps always execute in
pipeline order, regardless of argument order; duplicates run only once.
Only selected steps run: prerequisites are **not** run automatically.

Canonical step names (all legacy `M01`–`M08` / `M06` / `M04*` / `M05*` names
remain valid aliases in `--steps`):

| Step (`--<step>`) | Legacy alias | Script | Required data arguments |
|---|---|---|---|
| `inferExp` | `M01` | `M01_ENIGMA.R` | `--bulk-tpm`, `--enigma-ref`, `--cibersortx-frax` |
| `normExp` | `M02exp` | `M02_DataProcess_exp_fordb.R` | (uses `inferExp` output) |
| `peerCov` | `M02cov` | `M02_DataProcess_covariance_for_db.R` | `--covariates` |
| `trainModel` | `M03` | `M03_ModelTraining.R` | `--snp-annot-pattern`, `--genotype-pattern`, `--gene-annot` |
| `buildDB` | `M04db`, `M06` | `M04_Makedb.R` | (uses `trainModel` output) |
| `alignGWAS` | `M05align`, `M04align` | `M05_DataProcess_GWAS_align4METAL.R` | `--snp-annot-file`, `--gwas-dir` + GWAS mode |
| `metaGWAS` | `M05metal`, `M04metal` | `M05_DataProcess_GWAS_metal.sh` | (uses `alignGWAS` output) |
| `filterGWAS` | `M05filter`, `M04filter` | `M05_DataProcess_GWAS_metal_filter.R` | (uses `metaGWAS` output) |
| `extractLD` | `M06ld`, `M05ld` | `M06_LD.sh` | `--snp-annot-pattern`; also `--g1000-vcf-pattern`, `--g1000-samples` unless `--ld-dosage-dir` |
| `buildLD` | `M06cov`, `M05cov` | `M06_Makeld.sh` | (uses `buildDB` + `extractLD` output) |
| `runTWAS` | `M07` | `M07_TWAS.sh` | (uses `buildDB` + `buildLD` + `filterGWAS` / `--gwas-summary`) |
| `summarizeTWAS` | `M08` | `M08_ResultStatistic.R` | `--gene-pos` |

Other steps use upstream results and need no additional original data
arguments. For multiple steps, provide the union of their required inputs.
Existing tool and environment options still apply.

### GWAS modes

`alignGWAS` / `metaGWAS` / `filterGWAS` (`M05align`/`M05metal`/`M05filter`) require an **explicit** mode; changing
`--gwas-id` alone never silently switches what data is analyzed:

- `--gwas-profile sle-legacy-pbmc | sle-legacy-kidney` — bit-exact
  reproduction of the historical SLE implementation. **This profile
  reproduces known statistical errors of the old code** (the third study's
  BETA/SE column mis-projection, the swapped effect-allele declaration in
  METAL, whitespace-parsing column shifts). It exists so historical results
  remain reproducible; use it only for that. Its output names are fixed
  (`SLE_metal_096_917_183…`) and it forces `--gwas-filter-policy legacy`.
- `--gwas-manifest FILE` — the corrected standard path. `FILE` is a TSV with
  one row per GWAS study and the exact columns:

  ```
  study_id  path  build  effect_scale  chr_col  pos_col  ea_col  oa_col
  effect_mode  effect_col  se_col  p_col  eaf_col  n_mode  n_value  n_col
  ```

  `path` must be absolute; `effect_scale` is `linear` or `log_odds`;
  `effect_mode` is `beta_se`, `or_se` (SE of the log-OR) or `or_p` (SE
  back-computed from OR and P); `n_mode` is `none`, `constant` (`n_value`)
  or `column` (`n_col`). Unused optional columns are left empty. The driver
  prints a warning when the legacy profile is active; standard runs get the
  `<gwas_id>_meta` output prefix and filter policy `none` by default.

The active mode is recorded in `M04_gwas/gwas_run.json`, so re-running
`metaGWAS`/`filterGWAS`/`runTWAS` in the same output directory recovers it. `runTWAS` can also
run on any external summary with `--gwas-summary FILE` (cannot be combined
with the GWAS steps). `--gwas-id` is a label used in standard-mode output
names and reports only.

### Explicit analysis configuration

- `--cell-types` is always required and must match the reference cell-type names.
- `--peer-nk` and `--peer-keep` are required only for `peerCov` (`M02cov`). Use `0`/`0`
  for known covariates only; no PEER factors are inferred in this mode.
- `--twas-model-set unfiltered|filtered-signif` is required for `buildDB` (`M04db`) and `runTWAS` (`M07`).
- `--training-population` is optional metadata for `buildDB` (default `unspecified`).

Sample IDs default to exact matching and gene annotations default to the
standard schema. Old inputs requiring R-style sample names or a missing
`gene_id` column must explicitly select `--sample-id-policy legacy-rnames`
or `--gene-annot-format legacy-name-is-id`, respectively. Use separate output
directories for separate cohorts. `--ld-panel-id` names the LD panel independently
of the training population.

For comparison with the original analysis, these were the former parameter sets;
they are no longer selected automatically:

| Original analysis | Cell types | PEER estimated / kept | TWAS model set |
|---|---|---|---|
| PBMC | `B,CD4T,CD8T,NK,Mono,DC` | 2 / 2 | `unfiltered` |
| Glomerulus | `Bcell,Endo,MC,MNP,NK,NKT,Podo,Th_cell,CD8T` | 60 / 15 | `filtered-signif` |
| Tubule | `Bcell,CD8T,CDICA,CDICB,DCT,LOH,MNP,NK,NKT,PT,Th_cell,Plasmacytoid_DC` | 60 / 11 | `filtered-signif` |

Matching these settings alone does not guarantee identical results: input files,
GWAS mode, filtering, sample-ID policy and training seed must also match. `summarizeTWAS` (`M08`) omits
the `tissue` and `peer` label columns; association statistics and per-cell-type
corrections are unchanged by removing these labels.

### Parameters

`--cell-types` and `--output-dir` are always required. Data inputs are required
according to the selected steps above; tuning options retain their defaults:

| Parameter | Description |
|---|---|
| `--output-dir` | run directory; all results are written here |
| `--bulk-tpm` | bulk TPM expression matrix (genes x samples) |
| `--enigma-ref` | single-cell reference matrix (TPM, aggregated by cell type) |
| `--cibersortx-frax` | CIBERSORTx adjusted cell-fraction output |
| `--covariates` | covariates txt, first column = sample ID, then covariate columns (required for `peerCov`; format below) |
| `--cat-covs` | comma-separated categorical covariate columns to dummy-encode (default none) |
| `--min-level-count` | drop categorical levels with fewer samples (default 5) |
| `--peer-nk` | PEER hidden factors per cell type; explicit value required for `peerCov` |
| `--peer-keep` | PEER factors kept; explicit value required for `peerCov`; must be ≤ `--peer-nk` |
| `--snp-annot-pattern` | per-chromosome SNP annotation, with `{chr}` placeholder |
| `--genotype-pattern` | per-chromosome genotype dosage, with `{chr}` placeholder |
| `--gene-annot` | gene annotation (TSS-collapsed `start == end`, or with `strand`/`tss` column to resolve TSS) |
| `--snp-annot-file` | single-file genotype panel SNP annotation (for GWAS align) |
| `--gwas-dir` | directory with the three SLE GWAS files: `GCST011096_buildGRCh37.tsv`, `GCST90018917_buildGRCh37.tsv.gz`, `GCST90476183_hg19.tsv.gz` (legacy profiles) |
| `--g1000-vcf-pattern` | 1000 Genomes Phase 3 VCF pattern, with `{chr}` placeholder |
| `--g1000-samples` | 1000G sample list for the LD panel (EUR+EAS by default; pass an EUR-only list for a European LD panel) |
| `--gene-pos` | gene position annotation CSV (GRCh37, with gene_id/gene_name columns) |

Input contracts:

| Parameter | Description |
|---|---|
| `--sample-id-policy` | `exact` (default) or explicitly selected `legacy-rnames` |
| `--cov-layout` | `auto`, `samples` or `covariates` (default auto) |
| `--gene-annot-format` | `standard` (default) or explicitly selected `legacy-name-is-id` |
| `--genome-build` | `GRCh37` (default) or `GRCh38`; declared for this cohort's coordinates, no liftover is performed |
| `--training-population` | optional population label stored with the trained models (default `unspecified`) |
| `--ld-panel-id` | LD panel id used in dosage/covariance file names (default `EUR_EAS`) |
| `--twas-model-set` | `unfiltered` or `filtered-signif`; required for `buildDB`/`runTWAS` |
| `--maf` | MAF threshold for `trainModel` genotypes in `[0, 0.5]` (default `0.01`; manuscript used `0.05`) |
| `--maf-policy` | `filter` (drop variants `< --maf`) or `verify` (fail if input genotype violates `--maf`) |

#### Covariates file format

Tab- or space-separated with a header. Two layouts are auto-detected:

- **Recommended**: rows = samples, first column = sample ID, then one covariate
  column each. Numeric columns pass through as-is.
- **Legacy** (OneK1K-style): rows = covariates, columns = samples, first column
  = covariate names.

Columns listed in `--cat-covs` are treated as categorical and dummy-encoded,
and levels with fewer than `--min-level-count` samples are dropped.

```
id	age	sex	pc1	...	RIN	Sclerosis	Batch
sample1	54	0	0.01	...	7.2	0.15	3
sample2	47	1	-0.02	...	6.8	0.04	1
```
To reuse a phenotype table stored in an RData file (e.g. the kidney
`QTLData.RData` with `AAglomPhe2`/`AAtubPhe2`), export it once:

```r
load("QTLData.RData")
# put numeric covariates first, then the categorical ones you will list in --cat-covs
write.table(cbind(id = rownames(AAglomPhe2), AAglomPhe2),
            "glom_covariates.txt", sep = "\t", quote = FALSE)
```

then run with `--covariates glom_covariates.txt --cat-covs Pairend,Gender,Batch`.

Run control:

| Parameter | Description |
|---|---|
| `--alpha` | ENIGMA alpha (default 0.1) |
| `--tpm-min` | `normExp` gene filter: min inferred TPM (default 0.1) |
| `--min-pct` | `normExp` gene filter: min sample fraction, 0–1 (default 0.01) |
| `--min-samples` | `normExp` gene filter: min sample count floor (default 5) |
| `--all` | run all 12 steps; also the default when no mode is specified |
| `--steps STEP [STEP ...]` | run a subset in pipeline order; space/comma-separated; old `M01`–`M08` names are aliases |
| `--<step>` | single-step flags (e.g. `--inferExp`, `--summarizeTWAS`), combinable across steps |
| `--chromosomes` | comma-separated chromosomes (default 1–22) |
| `--cell-types` | explicit space/comma-separated cell-type list; required |
| `--cores` | unified parallel cores/workers across all pipeline steps (default 4) |
| `--seed` | explicit `trainModel` seed; forces serial inner training (`--training-cores 1`) |
| `--training-cores` | inner `trainModel` cores per `(cell_type, chr)` job (defaults from `--cores`) |
| `--ld-dosage-dir` | reuse pre-extracted 1000G dosages, skip `extractLD` (`M06ld`) extraction |
| `--gwas-id` | trait label (required for `--gwas-manifest` and standalone `summarizeTWAS` without `gwas_run.json`; auto-recovered from `gwas_run.json` when present) |
| `--gwas-summary` | explicit filtered GWAS summary for `runTWAS` (`M07`) |

Environment/tools (defaults match `setup.sh`):

| Parameter | Description |
|---|---|
| `--r-env` / `--peer-env` / `--py-env` | conda envs (default `R4.4.3` / `r3.6.3` / `imlatools`) |
| `--metaxcan-dir` | MetaXcan checkout (default `tools/MetaXcan`) |
| `--metal-bin` | METAL binary (default `tools/METAL/metal`) |
| `--bcftools` / `--bcf-plugins` | bcftools executable / plugin dir |

`trainModel` (`M03`) always uses the bundled training core, `script/helpers/gtex_v7_nested_cv_elnet.R`.

### Pipeline steps

Executed in dependency order (canonical names; legacy aliases in parentheses):

| Step | What it does |
|---|---|
| `inferExp` (`M01`) | ENIGMA trace-norm deconvolution with CIBERSORTx fractions |
| `normExp` (`M02exp`) | filter low-expression genes (TPM > `--tpm-min` in `max(--min-pct x N, --min-samples)` samples), per-gene inverse-normal transform per cell type |
| `peerCov` (`M02cov`) | PEER factors + covariates for model training (parallelized across cell types) |
| `trainModel` (`M03`) | nested-CV elastic-net GReX models (parallelized across cell types x chromosomes and genes) |
| `buildDB` (`M04db`, `M06`) | SQLite weight databases for S-PrediXcan (parallelized across cell types) |
| `alignGWAS` (`M05align`, `M04align`) | align GWAS studies to the genotype panel (allele flips, effect direction) |
| `metaGWAS` (`M05metal`, `M04metal`) | METAL inverse-variance meta-analysis |
| `filterGWAS` (`M05filter`, `M04filter`) | post-meta filtering (policy-dependent) |
| `extractLD` (`M06ld`, `M05ld`) | 1000G dosages aligned to the panel annotation (parallelized across chromosomes) |
| `buildLD` (`M06cov`, `M05cov`) | LD covariance matrices per cell type (parallelized across cell types) |
| `runTWAS` (`M07`) | S-PrediXcan association per cell type (parallelized across cell types) |
| `summarizeTWAS` (`M08`) | per-cell-type FDR / Bonferroni result summary |

`trainModel` dominates total runtime (roughly proportional to cell types x
chromosomes). `--cores` parallelizes both `(cell_type, chromosome)` jobs and
inner gene training, and `--steps` / `--chromosomes` can also split disjoint
subsets across separate sessions into the same output directory.

`buildDB` (`M04db`) accepts only the current 24-column `trainModel` summaries with
explicit `gene_name` and `gene_id`. Old 23-column summaries are rejected; rerun
`trainModel` rather than relying on automatic column repair. The significant-model thresholds
remain `zscore_pval < 0.05` and `rho_avg_squared > 0.01`. No intermediate
`{ct}_models.db` is created. Existing cached databases are not deleted.

### Outputs

```
{output-dir}/
  M01_deconvolution/    enigma.cse.trace_alpha_{alpha}.RData    deconvolved expression array
                        myenigma_alpha_{alpha}.RData           fitted ENIGMA object
  M02_expression/       {ct}_transformed_expression.txt         filtered + INT expression, samples x genes, sample IDs retained
  M02_covariates/       {ct}_covariates.txt                     covariates used by M03
                        {ct}_peer_res.RData                     raw PEER factors
  M03_models/summary/   Model_training_{ct}_chr{n}_model_summaries.txt
                        Model_training_{ct}_chr{n}_summary.txt
  M03_models/weights/   Model_training_{ct}_chr{n}_weights.txt
  M06_databases/        {ct}_models_unfiltered.db               all models; required for LD covariance and unfiltered TWAS
                        {ct}_models_filtered_signif.db          only when --twas-model-set filtered-signif is selected
  M04_gwas/             gwas_run.json                           recorded GWAS mode (profile or manifest)
                        SLE_metal_096_917_1831.tbl              legacy METAL meta-analysis result
                        SLE_metal_096_917_183_filt.tsv          legacy filtered meta statistics
                        {gwas_id}_meta1.tbl                     standard-mode METAL result
                        {gwas_id}_meta_filt.tsv                 standard-mode filtered statistics
  M05_ld/dosage/        chr{n}_{panel}.txt.gz                   1000G panel dosages
  M05_ld/cov_{panel}/   chr{n}.txt.gz                           dosages aligned to the panel
                        {ct}.txt.gz                             LD covariance per cell type
  M07_twas/             {ct}.csv                                per-cell-type S-PrediXcan associations
                        {ct}.log
  M08_summary/          allcelltype_result.csv                   all associations with FDR / Bonferroni columns
                        Bonfsig_allcelltype_result.csv          Bonferroni-significant genes
                        FDRsig_allcelltype_result.csv           FDR-significant genes
                        BonFDRsig_allcelltype_result.csv        union of both
```

`{ct}` = cell type, `{n}` = chromosome, `{alpha}` = the `--alpha` value,
`{panel}` = `--ld-panel-id`. Historical output directory names are unchanged
so existing caches remain reusable. Called helper scripts live in
`script/helpers/`; the M0-step entry scripts stay in `script/`.

## Reproducing historical results

The historical SLE numbers (TWAS, summaries) come from the legacy GWAS
profile, which preserves the old implementation's known errors. Reproduce
them with `--gwas-profile sle-legacy-pbmc` (or `sle-legacy-kidney`) and the
original inputs; do not mix legacy intermediates with standard-mode ones.
The corrected standard path (`--gwas-manifest`) can legitimately produce
different statistics for the same studies — that is expected, and the
difference is the fix.

## Part 2: LINCS drug prioritization (`medicell.py lincs`)

Part 2 is invoked separately via `medicell.py lincs` (or `medicell.py --lincs`).
It connects cell-type-resolved TWAS results and cell-type expression profiles
with LINCS L1000 perturbation signatures (`fgsea`) to score candidate compounds
for each cell type:

1. **`prepLincs` (`L01`, `script/L01_lincs_prepare.R`)**: Splits TWAS summary
   statistics by z-score sign (`pos`: `zscore > 0`, `neg`: `zscore < 0`) and
   builds a gene-symbol x cell-type expression matrix (either from `--exp-matrix`
   or by averaging `inferExp` ENIGMA `enigma.cse.test` across samples per cell type).
2. **`enrichLincs` (`L02`, `script/L02_lincs_enrichment.R`)**: Runs `fgsea` on
   the ranked directional TWAS genes against `data/lincs/LINCS_PRL_top250_human.gmt`
   (or a custom `--lincs-gmt`), annotates hits with L1000 Phase 1/2 perturbation
   metadata (`GSE92742` and `GSE70138` in `data/lincs/`), and filters hits by
   significance and NES sign (default `NES < 0`, `pval < 0.05`).
3. **`scoreLincs` (`L03`, `script/L03_lincs_drug_score.sh`)**: Selects the top
   perturbations per run, extracts their leading-edge genes and per-cell-type
   TWAS z-scores, computes cross-cell-type expression specificity (`expZscore`),
   and sums `RankScore = |TWAS_zscore| * expZscore` over each drug's
   leading-edge genes to produce a drug x cell-type prioritization matrix.

Called helper scripts reside in `script/helpers/` (`lincs_fgsea_lib.R`,
`lincsTissueSpecific.R`, `formatLeadingSpecific.pl`, `formatZscoreSpecific.pl`,
`getZscoreSpecificMatrix.R`). Because `setup.sh` provisions only Part 1 envs
(`R4.4.3`, `r3.6.3`, `imlatools`), Part 2 additionally requires an R environment
with `fgsea`, `dplyr`, `tidyr`, `data.table`, and `readxl` (for `.xlsx` inputs),
selected automatically (`R4.4.3` / `r_4.6.1`) or via `--lincs-env <env>`. Large
GMT files (`data/lincs/*.gmt`, >100 MB) are excluded from Git tracking; stage
them locally under `data/lincs/` or pass `--lincs-gmt`.

### Quick start

```bash
# 1) Continue from a Part 1 output directory (reuses M01 ENIGMA + M08 summary)
python3 medicell.py lincs \
    --output-dir run_out \
    --gene-pos data/gene_pos_grch37.82.csv \
    --label cohort --direction both

# 2) Run standalone from an external TWAS table (CSV/TSV/XLSX) + expression matrix
python3 medicell.py lincs \
    --output-dir lincs_out \
    --twas-file path/to/twas.xlsx --twas-sheet eGFR_pbmc_neg \
    --exp-matrix path/to/exp_matrix.csv \
    --trait eGFR --label pbmc --direction neg

# 3) Run with custom TWAS column names and a custom GMT (skipping Broad pert_info)
python3 medicell.py lincs \
    --output-dir lincs_out \
    --twas-file custom_twas.csv \
    --twas-gene-col Symbol --twas-celltype-col CellType \
    --twas-zscore-col Zscore --twas-padj-col FDR \
    --exp-matrix custom_exp.csv \
    --lincs-gmt custom.gmt \
    --lincs-pert-phase1 "" --lincs-pert-phase2 "" \
    --trait MyTrait --label MyCohort --direction pos

# 4) Score directly from a pre-computed LINCS enrichment CSV (skip L02)
python3 medicell.py lincs \
    --output-dir lincs_out \
    --twas-file path/to/twas.csv --exp-matrix path/to/exp_matrix.csv \
    --lincs-sig-csv path/to/LINCS_sig.csv \
    --trait MyTrait --label MyCohort --direction pos
```

### Key options (`python3 medicell.py lincs --help`)

| Option | Default | Description |
|---|---|---|
| `--output-dir` | required | Output root (may reuse Part 1 `--output-dir` containing `M01`/`M08`). |
| `--cores` | `4` | Unified parallel worker/core budget across all steps (`--training-cores` / `--lincs-cores` still work as overrides). |
| `--steps` / `--<step>` | `prepLincs,enrichLincs,scoreLincs` | Steps to run (`L01,L02,L03` aliases work; defaults to `prepLincs,scoreLincs` with `--lincs-sig-csv`). |
| `--trait` | auto / `trait` | Output trait tag; auto-reads `gwas_id` from `M04_gwas/gwas_run.json` if present. |
| `--label` / `--tissue` | `cohort` | Cohort/compartment label used in output file names (no tissue presets). |
| `--direction` | `pos` | TWAS z-score direction: `pos` (`>0`), `neg` (`<0`), or `both`. |
| `--cell-types` | all | Optional comma/space-separated subset of cell types (`scoreLincs` requires `>= 2`). |
| `--twas-file` | `M08_summary/...` | TWAS summary table (`.csv`, `.tsv`, `.txt`, or `.xlsx` with `--twas-sheet`). |
| `--twas-gene-col`, `--twas-zscore-col`, `--twas-padj-col`, `--twas-celltype-col` | auto-detect | Explicit column names if your TWAS table uses non-standard headers. |
| `--exp-matrix` | from `M01` | Pre-built gene-symbol x cell-type CSV (rows = genes, columns = cell types). |
| `--lincs-gmt` | `data/lincs/LINCS_PRL_top250_human.gmt` | Perturbation/gene-set GMT file (`--lincs-pert-phase1/2 ""` disables extra annotation). |
| `--lincs-sig-col`, `--lincs-sig-cutoff`, `--lincs-nes-filter` | `pval`, `0.05`, `neg` | Significance filter (`pval` or `padj`) and NES direction (`neg`, `pos`, `all`). |
| `--top-n`, `--lincs-rank-by`, `--exclude-pert-pattern` | `50`, `pval`, `^BRD-\|trt_sh\.css` | Top perturbations kept in `L03`; pass `--exclude-pert-pattern ""` to keep all. |

### LINCS outputs

```text
<output-dir>/
  L01_inputs/       {trait}_{label}_{dir}.csv          directional TWAS subset for fgsea
                    {trait}_{label}_twas_all.tsv       normalized Gene/CellType/Zscore table
                    {label}_exp_matrix.csv             gene x cell-type expression matrix
  L02_enrichment/   LINCS_all_{trait}_{label}_{dir}.csv  all fgsea-tested perturbations
                    LINCS_sig_{trait}_{label}_{dir}.csv  filtered LINCS hits
  L03_drug_scores/  {label}_{trait}_{dir}.SpeciGenes.txt
                    {label}_{trait}_{dir}.LeadingGenesPool.txt
                    {label}_{trait}_{dir}.LeadingGenesDrugPairs.txt
                    {label}_{trait}_{dir}.twas_zscore.txt
                    {label}_{trait}_{dir}.drug_celltype_score.RData
                    {label}_{trait}_{dir}.drug_celltype_score.csv   final drug x cell-type score matrix
```

## References

- ENIGMA: deconvolution of bulk RNA-seq with single-cell reference
- MetaXcan / S-PrediXcan: https://github.com/hakyimlab/MetaXcan
- METAL: https://github.com/statgen/Metal
- PEER: https://github.com/PMBio/peer
- 1000 Genomes Project Phase 3
- LINCS L1000 (GSE92742, GSE70138): https://www.lincsproject.org/
- fgsea: https://bioconductor.org/packages/fgsea/
