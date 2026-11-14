#!/bin/bash
# M06ld: GT non-reference counts, AF before sample subset, explicit panel label.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
HELPER="${HERE}/helpers/ld_contracts.py"
DOSAGE_DIR="${MC_LD_DOSAGE_DIR:?}"
COV_DIR="${MC_LD_COV_DIR:?}"
ANNOT_PATTERN="${MC_SNP_ANNOT_PATTERN:?}"
PANEL="${MC_LD_PANEL_ID:-EUR_EAS}"
[[ "$PANEL" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo 'Invalid LD panel label' >&2; exit 1; }
read -r -a CHROMS <<< "${MC_CHROMS:?}"
for chr in "${CHROMS[@]}"; do
    [[ "$chr" =~ ^([1-9]|1[0-9]|2[0-2])$ ]] || { echo 'MC_CHROMS must contain space-separated autosomes' >&2; exit 1; }
done
mkdir -p "$DOSAGE_DIR" "$COV_DIR"
WORK="$(mktemp -d "${COV_DIR}/.ld.XXXXXX")"
trap 'rm -rf -- "$WORK"' EXIT
if [[ "${MC_LD_SKIP_EXTRACT:-0}" != 1 ]]; then
    BCF="${MC_BCFTOOLS:-bcftools}"
    if [[ -n "${MC_BCF_PLUGINS:-}" ]]; then export BCFTOOLS_PLUGINS="$MC_BCF_PLUGINS"; fi
    for chr in "${CHROMS[@]}"; do
        INPUT="${MC_G1000_VCF_PATTERN:?}"
        INPUT="${INPUT/\{chr\}/$chr}"
        ANNOT="${ANNOT_PATTERN/\{chr\}/$chr}"
        python3 "$HELPER" sites "$ANNOT" "$chr" "$WORK/sites.tsv"
        "$BCF" query -l "$INPUT" > "$WORK/available"
        python3 "$HELPER" samples "${MC_G1000_SAMPLES:?}" "$WORK/available" "$WORK/requested"
        # Do not move fill-tags after subset: that would change historical AF semantics.
        "$BCF" +fill-tags "$INPUT" -Ou -- -t AF | \
            "$BCF" view -m2 -M2 -v snps -T "$WORK/sites.tsv" -S "$WORK/requested" -Ob -o "$WORK/subset.bcf"
        "$BCF" query -l "$WORK/subset.bcf" > "$WORK/actual"
        cmp "$WORK/requested" "$WORK/actual" || { echo 'VCF sample order mismatch' >&2; exit 1; }
        n="$(wc -l < "$WORK/requested")"
        "$BCF" query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF[\t%GT]\n' "$WORK/subset.bcf" | \
            python3 "$HELPER" gt "$WORK/chr${chr}_dosage.gz" "$chr" "$n"
        mv "$WORK/chr${chr}_dosage.gz" "$DOSAGE_DIR/chr${chr}_${PANEL}.txt.gz"
        cp "$WORK/requested" "$DOSAGE_DIR/chr${chr}_${PANEL}.samples.txt"
    done
fi
for chr in "${CHROMS[@]}"; do
    ANNOT="${ANNOT_PATTERN/\{chr\}/$chr}"
    python3 "$HELPER" align "$ANNOT" "$DOSAGE_DIR/chr${chr}_${PANEL}.txt.gz" "$WORK/chr${chr}.txt.gz" "$chr"
    mv "$WORK/chr${chr}.txt.gz" "$COV_DIR/chr${chr}.txt.gz"
    printf 'panel_id\t%s\ngenome_build\t%s\naf_semantics\tbefore_sample_subset\nsample_order\t%s\n' \
        "$PANEL" "${MC_GENOME_BUILD:?}" "${MC_LD_SKIP_EXTRACT:-0}: 1=preextracted order not independently verifiable; 0=VCF/requested verified" > "$COV_DIR/chr${chr}.metadata.tsv"
done
python3 "$HELPER" check-many "$COV_DIR" "${CHROMS[@]}"
