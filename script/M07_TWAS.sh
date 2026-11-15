#!/bin/bash
# M07: explicit model strategy; preserve vendor math and GWAS column contract.
set -euo pipefail
METAXCAN_DIR="${MC_METAXCAN_DIR:?}"
DB_DIR="${MC_DB_DIR:?}"
LD_COV_DIR="${MC_LD_COV_DIR:?}"
GWAS_FILE="${MC_GWAS_FILT:?}"
TWAS_OUT="${MC_TWAS_OUT:?}"
case "${MC_TWAS_MODEL_SET:?Explicit unfiltered or filtered-signif model set required}" in
    unfiltered) DB_SUFFIX="_models_unfiltered.db" ;;
    filtered-signif) DB_SUFFIX="_models_filtered_signif.db" ;;
    *) echo 'MC_TWAS_MODEL_SET must be unfiltered or filtered-signif' >&2; exit 1 ;;
esac
read -r -a CELLTYPES <<< "${MC_CELLTYPES:?}"
[[ -s "$GWAS_FILE" ]] || { echo "Missing GWAS: $GWAS_FILE" >&2; exit 1; }
mkdir -p "$TWAS_OUT"
# Check every selected CT before touching any output/log.
for ct in "${CELLTYPES[@]}"; do
    [[ "$ct" =~ ^[A-Za-z0-9_.-]+$ ]] || exit 1
    for file in "$DB_DIR/$ct$DB_SUFFIX" "$LD_COV_DIR/$ct.txt.gz"; do
        [[ -s "$file" ]] || { echo "Missing TWAS input: $file" >&2; exit 1; }
    done
    [[ ! -e "$TWAS_OUT/$ct.csv" && ! -e "$TWAS_OUT/$ct.log" ]] || { echo "Refusing existing TWAS output/log: $ct" >&2; exit 1; }
done
for ct in "${CELLTYPES[@]}"; do
    echo "$ct | ${MC_TWAS_MODEL_SET}"
    python3 "$METAXCAN_DIR/software/SPrediXcan.py" \
        --model_db_path "$DB_DIR/$ct$DB_SUFFIX" \
        --covariance "$LD_COV_DIR/$ct.txt.gz" \
        --keep_non_rsid --model_db_snp_key VarID \
        --gwas_file "$GWAS_FILE" \
        --snp_column MarkerName --effect_allele_column Allele1 \
        --non_effect_allele_column Allele2 --beta_column Effect \
        --pvalue_column P-value --se_column StdErr --input_pvalue_fix 0 \
        --output_file "$TWAS_OUT/$ct.csv" > "$TWAS_OUT/$ct.log" 2>&1
    python3 - "$TWAS_OUT/$ct.csv" <<'CHECK'
import csv, math, sys
with open(sys.argv[1]) as f:
    reader = csv.DictReader(f)
    if not {'gene','gene_name','zscore','pvalue'} <= set(reader.fieldnames or []):
        sys.exit('Invalid TWAS output schema')
    count = 0
    for row in reader:
        p, z = float(row['pvalue']), float(row['zscore'])
        if not row['gene'] or not math.isfinite(z) or not math.isfinite(p) or not 0 <= p <= 1:
            sys.exit('Invalid TWAS output row')
        count += 1
    if not count:
        sys.exit('Empty TWAS output')
CHECK
done
