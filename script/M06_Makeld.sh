#!/usr/bin/env bash
# M06cov: isolate only this run's selected dosage chromosomes for vendor scanning.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
HELPER="${HERE}/helpers/ld_contracts.py"
PARALLEL_RUN="${HERE}/helpers/parallel_run.py"
METAXCAN_DIR="${MC_METAXCAN_DIR:?}"
DB_DIR="${MC_DB_DIR:?}"
LD_COV_DIR="${MC_LD_COV_DIR:?}"
CORES="${MC_CORES:-4}"
read -r -a CHROMS <<< "${MC_CHROMS:?}"
read -r -a CELLTYPES <<< "${MC_CELLTYPES:?}"
for chr in "${CHROMS[@]}"; do
    [[ "$chr" =~ ^([1-9]|1[0-9]|2[0-2])$ ]] || { echo 'Invalid chromosome selection' >&2; exit 1; }
done
for ct in "${CELLTYPES[@]}"; do
    [[ "$ct" =~ ^[A-Za-z0-9_.-]+$ ]] || exit 1
    [[ -s "$DB_DIR/${ct}_models_unfiltered.db" ]] || { echo "Missing DB: $ct" >&2; exit 1; }
    [[ ! -e "$LD_COV_DIR/$ct.txt.gz" && ! -e "$LD_COV_DIR/$ct.cov.log" ]] || { echo "Refusing existing covariance/log: $ct" >&2; exit 1; }
done
python3 "$HELPER" check-many "$LD_COV_DIR" "${CHROMS[@]}"
WORK="$(mktemp -d "${LD_COV_DIR}/.cov.XXXXXX")"
trap 'rm -rf -- "$WORK"' EXIT
mkdir "$WORK/input" "$WORK/output"
for chr in "${CHROMS[@]}"; do
    ln -s "$(realpath "$LD_COV_DIR/chr${chr}.txt.gz")" "$WORK/input/chr${chr}.txt.gz"
done

export METAXCAN_DIR DB_DIR LD_COV_DIR WORK HELPER
cov_ct() {
    local ct="$1"
    python3 "$METAXCAN_DIR/software/M01_covariances_correlations.py" \
        --weight_db "$DB_DIR/${ct}_models_unfiltered.db" \
        --input_folder "$WORK/input/" --delimiter $'\t' \
        --covariance_output "$WORK/output/$ct.txt.gz" > "$LD_COV_DIR/$ct.cov.log" 2>&1
    python3 "$HELPER" covariance "$WORK/output/$ct.txt.gz"
    mv "$WORK/output/$ct.txt.gz" "$LD_COV_DIR/$ct.txt.gz"
}
export -f cov_ct

python3 "$PARALLEL_RUN" --cores "$CORES" --fn cov_ct -- "${CELLTYPES[@]}"
