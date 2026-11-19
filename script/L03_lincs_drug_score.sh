#!/usr/bin/env bash
# L03: Filter top LINCS perturbations, extract leading-edge gene TWAS z-scores,
# and compute drug x cell-type prioritization matrices.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
HELPERS="${SCRIPT_DIR}/helpers"

: "${MC_LINCS_L03_DIR:?Missing MC_LINCS_L03_DIR}"
: "${MC_LINCS_SIG_CSV:?Missing MC_LINCS_SIG_CSV}"
: "${MC_LINCS_TWAS_TSV:?Missing MC_LINCS_TWAS_TSV}"
: "${MC_LINCS_EXP_CSV:?Missing MC_LINCS_EXP_CSV}"
: "${MC_LINCS_TRAIT:?Missing MC_LINCS_TRAIT}"
: "${MC_LINCS_LABEL:?Missing MC_LINCS_LABEL}"
: "${MC_LINCS_DIRECTION:?Missing MC_LINCS_DIRECTION}"

TOP_N="${MC_LINCS_TOP_N:-50}"
EXCLUDE_PATTERN="${MC_LINCS_EXCLUDE_PATTERN-^BRD-|trt_sh\\.css}"
if [ -z "${EXCLUDE_PATTERN}" ]; then
  EXCLUDE_PATTERN="NONE"
fi
RANK_BY="${MC_LINCS_RANK_BY:-pval}"
TWAS_TAG="${MC_LINCS_TWAS_LABEL:-${MC_LINCS_TRAIT}_${MC_LINCS_LABEL}_${MC_LINCS_DIRECTION}}"
PREFIX="${MC_LINCS_LABEL}_${MC_LINCS_TRAIT}_${MC_LINCS_DIRECTION}"

mkdir -p "${MC_LINCS_L03_DIR}"

SPECI_TXT="${MC_LINCS_L03_DIR}/${PREFIX}.SpeciGenes.txt"
TOP_RDATA="${MC_LINCS_L03_DIR}/${PREFIX}.top_specific.RData"
POOL_TXT="${MC_LINCS_L03_DIR}/${PREFIX}.LeadingGenesPool.txt"
PAIRS_TXT="${MC_LINCS_L03_DIR}/${PREFIX}.LeadingGenesDrugPairs.txt"
ZSCORE_TXT="${MC_LINCS_L03_DIR}/${PREFIX}.twas_zscore.txt"
RES_RDATA="${MC_LINCS_L03_DIR}/${PREFIX}.drug_celltype_score.RData"
RES_CSV="${MC_LINCS_L03_DIR}/${PREFIX}.drug_celltype_score.csv"

Rscript "${HELPERS}/lincsTissueSpecific.R" \
  "${MC_LINCS_SIG_CSV}" \
  "${TWAS_TAG}" \
  "${SPECI_TXT}" \
  "${TOP_RDATA}" \
  "${TOP_N}" \
  "${EXCLUDE_PATTERN}" \
  "${RANK_BY}"

perl "${HELPERS}/formatLeadingSpecific.pl" \
  --input "${SPECI_TXT}" \
  --pool-out "${POOL_TXT}" \
  --pairs-out "${PAIRS_TXT}"

Z_ARGS=(--pool "${POOL_TXT}" --twas "${MC_LINCS_TWAS_TSV}" --out "${ZSCORE_TXT}")
if [ -n "${MC_LINCS_L03_GENE_COL:-}" ]; then
  Z_ARGS+=(--gene-col "${MC_LINCS_L03_GENE_COL}")
fi
if [ -n "${MC_LINCS_L03_CT_COL:-}" ]; then
  Z_ARGS+=(--ct-col "${MC_LINCS_L03_CT_COL}")
fi
if [ -n "${MC_LINCS_L03_Z_COL:-}" ]; then
  Z_ARGS+=(--z-col "${MC_LINCS_L03_Z_COL}")
fi
perl "${HELPERS}/formatZscoreSpecific.pl" "${Z_ARGS[@]}"

Rscript "${HELPERS}/getZscoreSpecificMatrix.R" \
  "${ZSCORE_TXT}" \
  "${MC_LINCS_EXP_CSV}" \
  "${PAIRS_TXT}" \
  "${RES_RDATA}" \
  "${RES_CSV}"
