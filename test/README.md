# MediCell End-to-End Test

This directory contains a lightweight, de-identified chromosome 10 test dataset
and automated verification scripts for the full MediCell workflow (`M01`–`M08`
TWAS pipeline plus `L01`–`L03` LINCS therapeutic prioritization).

## Command

```bash
bash test/run_full.sh
```

`test/run_full.sh` runs in single-core mode (`--cores 1 --training-cores 1`) on
chromosome 10 across two target cell types (`CD4T,B`), verifies all TWAS
outputs via `test/check_results.py`, and then invokes `test/run_lincs.sh` to
verify the LINCS workflow (including the `M01` + `M08` -> `L01`–`L03` handoff
when `data/ref/lincs/` is unpacked).

## Staged chromosome 10 dataset (`test/data/`)

| Item | Summary |
|---|---|
| Discovery cohort | 50 randomly sampled samples (`seed=2026`), anonymized as `S001`–`S050` |
| Bulk and ENIGMA reference | 150 genes; 6 reference cell types (`Mono,NK,CD4T,B,CD8T,DC`) |
| Known covariates | Sex, age, and 10 genotype principal components |
| PEER | Two factors estimated and retained (`--peer-nk 2 --peer-keep 2`) |
| Target cell types & chromosome | `CD4T,B`, chromosome 10 (`--chromosomes 10`) |
| Target genes | 14 protein-coding `chr10` genes (`C10orf128`, `PPA1`, `MGMT`, `AKR1C2`, `C10orf32`, `ADARB2`, `GDI2`, `STOX1`, `ACTA2`, `VCL`, `CCSER2`, `GSTO1`, `NEBL`, `DPYSL4`) |
| Genotypes | 197 non-palindromic `chr10` SNPs |
| GWAS | Two-study GRCh37 `beta_se` manifest (`bentham`, `jiang`, 197 SNPs each) |
| LD reference | 30 randomly sampled 1000 Genomes samples (15 EUR + 15 EAS, `seed=2026`) |
| Training seed / cores | `42` / `1` |
| Input directory size | ~250 KiB (`manifest.json` stores file SHA-256 checksums) |

## Verification checks

- `test/check_results.py` validates `test/data/manifest.json` SHA-256 checksums,
  execution of all 12 pipeline steps (`inferExp` through `summarizeTWAS`),
  non-empty `CD4T` and `B` SQLite databases (`unfiltered`), finite LD
  covariances, non-empty `chr10` TWAS results, and `M08_summary/` tables.
- `test/check_lincs_results.py` checks standalone directional `prepLincs ->
  enrichLincs -> scoreLincs` score symmetry, Broad `LINCS_PRL_top250_human.gmt`
  annotation/scoring (when `data/ref/lincs/` is present), single-cell-type
  rejection, and the full 3-step `M01` + `M08` -> `L01`–`L03` handoff yielding
  a populated `CD4T`/`B` drug-scoring matrix.

These small fixtures are intended purely as an executable software smoke test,
not for biological association inference.
