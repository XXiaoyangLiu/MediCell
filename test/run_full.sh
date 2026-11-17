#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
mkdir -p test/results
OUT="$(mktemp -d "$ROOT/test/results/full.XXXXXX")"
printf 'Output directory: %s\n' "$OUT"
cp test/data/manifest.json "$OUT/input_manifest.json"
sed "s|test/data/|$ROOT/test/data/|g" test/data/gwas_manifest.tsv > "$OUT/gwas_manifest.tsv"
SECONDS=0
python3 -u medicell.py --all --output-dir "$OUT" \
    --cores 1 --training-cores 1 \
    --cell-types CD4T,B --chromosomes 10 --alpha 0.1 \
    --peer-nk 2 --peer-keep 2 --cov-layout samples \
    --sample-id-policy exact --gene-annot-format standard \
    --training-population EUR --twas-model-set unfiltered \
    --seed 42 --genome-build GRCh37 \
    --bulk-tpm test/data/bulk_tpm.tsv \
    --enigma-ref test/data/reference_tpm.tsv \
    --cibersortx-frax test/data/fractions.tsv \
    --covariates test/data/covariates.tsv \
    --snp-annot-pattern 'test/data/snp_annot.chr{chr}.txt' \
    --genotype-pattern 'test/data/genotype.chr{chr}.txt' \
    --gene-annot test/data/gene_annot.tsv \
    --snp-annot-file test/data/gwas_panel.tsv \
    --gwas-manifest "$OUT/gwas_manifest.tsv" --gwas-id SLE \
    --g1000-vcf-pattern 'test/data/reference.chr{chr}.vcf.gz' \
    --g1000-samples test/data/ld_samples.txt --ld-panel-id MINI_EUR_EAS \
    --gene-pos test/data/gene_names.csv 2>&1 | tee "$OUT/run.log"
python3 test/check_results.py "$OUT"
bash test/run_lincs.sh "$OUT" 2>&1 | tee "$OUT/lincs.log"
printf 'Elapsed seconds: %s\n' "$SECONDS" | tee "$OUT/duration.txt"
