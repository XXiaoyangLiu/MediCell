# Compact real-data smoke test

This private fixture exercises all 12 MediCell steps, from ENIGMA to the final
TWAS summaries. It does not reuse a fitted ENIGMA object, PEER factors, trained
weights, METAL results or precomputed LD covariance.

## Run

From the repository root on the analysis server:

```bash
bash test/run_full.sh
```

Each invocation creates a new `test/results/full.XXXXXX/` directory, saves the
input manifest, runs `medicell.py --all`, and checks the resulting files. The
launcher stops on failure. Existing runs are never overwritten.

The verified run `test/results/full.tXLrIR/` completed all 12 steps in 110 seconds
on the current server. Runtime on other machines may differ. Its four M08 CSVs
were checked against the preceding successful run, `full.PoIYJp`.

## Fixture size

| Component | Size |
|---|---|
| Discovery cohort | 40 matched real samples |
| Bulk and aggregated reference | 128 genes; six reference cell types |
| Known covariates | Sex, age and six genotype principal components |
| PEER | Two factors estimated and retained |
| Training | CD4T, chromosome 22, one gene: FAM118A / ENSG00000100376 |
| Genotypes | Eight real nonpalindromic SNPs |
| GWAS | Bentham and Jiang studies; eight corresponding records each |
| LD reference | 20 real 1000 Genomes samples: 10 EUR and 10 EAS |
| Coordinates | GRCh37 |
| Training seed / cores | 42 / 1 |
| Input directory | Approximately 121 KiB, including VCF index and provenance |

Forty is the current default training code's minimum accepted sample count
(10-fold cross-validation with at least four samples per fold). It is not a
recommended scientific sample size.

Panels containing 32 and 64 expression genes were also tested. Both produced
zero training weights and stopped at database construction. The 128-gene panel
passed twice with the same 40 samples. This is the smallest successful panel
among the sizes tested, not a claim that 128 genes is an absolute minimum.
Failed runs and their input snapshots are retained under `test/results/` and
`test/attempts/`, respectively.

## Inputs and provenance

`test/data/manifest.json` records source paths, source sizes and modification
times, staged-file SHA-256 checksums and all preparation decisions. Every sample
is selected by a matching ID, never by assuming that unrelated matrix row
orders correspond. The discovery samples are the first 40 shared IDs in the
original bulk header order.

Preparation preserves the source expression and dosage values. It explicitly:

- maps fraction IDs such as `X1_1` to their uniquely matching bulk IDs;
- adds the original Plasma fraction to B to match the six-type reference;
- retains all six reference cell types for M01, while downstream analysis is CD4T;
- preserves the source training annotation's TSS-point coordinates, adding an
  explicit `gene_id` equal to its existing Ensembl identifier;
- selects up to eight largest-absolute-weight SNPs from a pre-existing CD4T
  model, requiring complete real discovery, LD and two-study GWAS coverage;
- uses the pre-existing model only to select the fixture, not as a pipeline input.

The GWAS manifest uses the standard `beta_se` path with a common log-odds scale.
The standard default post-meta policy is `none`; the filtering step still runs.
The tiny staged VCF defines the test LD panel, so its allele frequencies need
not equal those from the complete 1000 Genomes panel.

`prepare_data.py` regenerates the fixture from the authorized server sources.
It refuses to overwrite a nonempty `test/data/` directory. Regeneration requires
access to those sources and is not required for repeated test runs.
`resize_expression.py` was used only to diagnose the smaller, unsuccessful panels.

## Checks and limitations

The launcher runs `check_results.py` after M08. The check verifies all 12 child
steps executed, a newly trained nonempty database, two normalized GWAS studies,
finite covariance and TWAS outputs, complete eight-SNP coverage, the gene-name
mapping, and all four final summaries. To check an existing successful run:

```bash
python3 test/check_results.py test/results/full.tXLrIR
```

The fixture is deliberately selected for an executable smoke test. It is not an
unbiased evaluation set, historical reproduction, association discovery or
proof of predictive accuracy. Multiple-testing corrections here cover only one
tested gene. The small LD panel is not suitable for scientific conclusions.

Warnings retained in the successful run:

- ENIGMA reports that `scater` and `preprocessCore` are unavailable. The selected
  `Normalize=FALSE` execution path nevertheless completes. These warnings were
  not hidden, and this test does not validate a fresh dependency installation.
- MetaXcan reports that GWAS heritability and sample-size calibration inputs
  were not supplied, so its inflation calibration was not performed.

These inputs contain authorized cohort-derived data. Keep `test/data/`,
`test/results/`, `test/attempts/` and logs private; they are ignored by Git.
Do not redistribute them or include them in a portable release without separate
permission. Only test scripts and documentation are synchronized locally.

## LINCS therapeutic prioritization test (`test/run_lincs.sh`)

`test/run_lincs.sh` validates the LINCS therapeutic prioritization pipeline (`prepLincs -> enrichLincs -> scoreLincs`) and is invoked automatically at the end of `test/run_full.sh`:

```bash
# Run standalone LINCS smoke tests
bash test/run_lincs.sh

# Or run the full TWAS + LINCS suite (including inferExp/summarizeTWAS -> LINCS handoff)
bash test/run_full.sh
```

`test/check_lincs_results.py` checks:
- self-contained multi-gene `prepLincs -> enrichLincs -> scoreLincs` (`--direction both --cores 4`) and score symmetry;
- Broad `LINCS_PRL_top250_human.gmt` + Phase 1/2 perturbation annotation join when `data/lincs/` GMT files are present;
- Upstream `inferExp` (`M01`) ENIGMA `.RData` + `summarizeTWAS` (`M08`) handoff when a TWAS output directory is passed.
