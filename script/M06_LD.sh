#!/usr/bin/env bash
# M06ld: GT non-reference counts, AF before sample subset, explicit panel label.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
HELPER="${HERE}/helpers/ld_contracts.py"
PARALLEL_RUN="${HERE}/helpers/parallel_run.py"
DOSAGE_DIR="${MC_LD_DOSAGE_DIR:?}"
COV_DIR="${MC_LD_COV_DIR:?}"
ANNOT_PATTERN="${MC_SNP_ANNOT_PATTERN:?}"
PANEL="${MC_LD_PANEL_ID:-EUR_EAS}"
CORES="${MC_CORES:-4}"
[[ "$PANEL" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo 'Invalid LD panel label' >&2; exit 1; }
read -r -a CHROMS <<< "${MC_CHROMS:?}"
for chr in "${CHROMS[@]}"; do
    [[ "$chr" =~ ^([1-9]|1[0-9]|2[0-2])$ ]] || { echo 'MC_CHROMS must contain space-separated autosomes' >&2; exit 1; }
done
mkdir -p "$DOSAGE_DIR" "$COV_DIR"
WORK="$(mktemp -d "${COV_DIR}/.ld.XXXXXX")"
trap 'rm -rf -- "$WORK"' EXIT
BCF="${MC_BCFTOOLS:-bcftools}"
if [[ -n "${MC_BCF_PLUGINS:-}" ]]; then export BCFTOOLS_PLUGINS="$MC_BCF_PLUGINS"; fi

export HELPER DOSAGE_DIR COV_DIR ANNOT_PATTERN PANEL WORK BCF
export MC_G1000_VCF_PATTERN MC_G1000_SAMPLES MC_GENOME_BUILD MC_LD_SKIP_EXTRACT

extract_chr() {
    local chr="$1"
    local cwork="$WORK/extract_chr${chr}"
    mkdir -p "$cwork"
    local INPUT="${MC_G1000_VCF_PATTERN:?}"
    INPUT="${INPUT/\{chr\}/$chr}"
    local ANNOT="${ANNOT_PATTERN/\{chr\}/$chr}"
    python3 "$HELPER" sites "$ANNOT" "$chr" "$cwork/sites.tsv"
    "$BCF" query -l "$INPUT" > "$cwork/available"
    python3 "$HELPER" samples "${MC_G1000_SAMPLES:?}" "$cwork/available" "$cwork/requested"
    # Do not move fill-tags after subset: that would change historical AF semantics.
    "$BCF" +fill-tags "$INPUT" -Ou -- -t AF | \
        "$BCF" view -m2 -M2 -v snps -T "$cwork/sites.tsv" -S "$cwork/requested" -Ob -o "$cwork/subset.bcf"
    "$BCF" query -l "$cwork/subset.bcf" > "$cwork/actual"
    cmp "$cwork/requested" "$cwork/actual" || { echo "VCF sample order mismatch (chr${chr})" >&2; exit 1; }
    local n
    n="$(wc -l < "$cwork/requested")"
    "$BCF" query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF[\t%GT]\n' "$cwork/subset.bcf" | \
        python3 "$HELPER" gt "$cwork/chr${chr}_dosage.gz" "$chr" "$n"
    mv "$cwork/chr${chr}_dosage.gz" "$DOSAGE_DIR/chr${chr}_${PANEL}.txt.gz"
    cp "$cwork/requested" "$DOSAGE_DIR/chr${chr}_${PANEL}.samples.txt"
}
export -f extract_chr

align_chr() {
    local chr="$1"
    local cwork="$WORK/align_chr${chr}"
    mkdir -p "$cwork"
    local ANNOT="${ANNOT_PATTERN/\{chr\}/$chr}"
    python3 "$HELPER" align "$ANNOT" "$DOSAGE_DIR/chr${chr}_${PANEL}.txt.gz" "$cwork/chr${chr}.txt.gz" "$chr"
    mv "$cwork/chr${chr}.txt.gz" "$COV_DIR/chr${chr}.txt.gz"
    printf 'panel_id\t%s\ngenome_build\t%s\naf_semantics\tbefore_sample_subset\nsample_order\t%s\n' \
        "$PANEL" "${MC_GENOME_BUILD:?}" "${MC_LD_SKIP_EXTRACT:-0}: 1=preextracted order not independently verifiable; 0=VCF/requested verified" > "$COV_DIR/chr${chr}.metadata.tsv"
}
export -f align_chr

if [[ "${MC_LD_SKIP_EXTRACT:-0}" != 1 ]]; then
    python3 "$PARALLEL_RUN" --cores "$CORES" --fn extract_chr -- "${CHROMS[@]}"
fi
python3 "$PARALLEL_RUN" --cores "$CORES" --fn align_chr -- "${CHROMS[@]}"
python3 "$HELPER" check-many "$COV_DIR" "${CHROMS[@]}"
