#!/bin/bash
# Build lean MediCell release archives:
#   release/medicell-software.tar.gz  (extracts to MediCell/)
#   release/medicell-ref-data.tar.gz  (extracts to ref/, placed under MediCell/data/)
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT_DIR="${1:-${REPO_DIR}/release}"
mkdir -p "${OUT_DIR}"

SOFT_TGZ="${OUT_DIR}/medicell-software.tar.gz"
DATA_TGZ="${OUT_DIR}/medicell-ref-data.tar.gz"

# Sanity checks before packing
for req in \
    "${REPO_DIR}/medicell.py" \
    "${REPO_DIR}/README.md" \
    "${REPO_DIR}/envs/R4.4.3.yml" \
    "${REPO_DIR}/envs/r3.6.3.yml" \
    "${REPO_DIR}/envs/imlatools.yml" \
    "${REPO_DIR}/tools/METAL/metal" \
    "${REPO_DIR}/tools/MetaXcan/software/SPrediXcan.py" \
    "${REPO_DIR}/tools/MetaXcan/software/M01_covariances_correlations.py" \
    "${REPO_DIR}/tools/R_libs/R4.4.3/ENIGMA/DESCRIPTION" \
    "${REPO_DIR}/tools/R_libs/r3.6.3/peer/libs/peer.so" \
    "${REPO_DIR}/data/ref/README.md" \
    "${REPO_DIR}/data/ref/lincs/LINCS_PRL_top250_human.gmt" \
    "${REPO_DIR}/data/ref/lincs/LINCS_PRL_top250.gmt" \
    "${REPO_DIR}/data/ref/lincs/GSE92742_Broad_LINCS_pert_info.txt.gz" \
    "${REPO_DIR}/data/ref/lincs/GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz"
do
    [ -e "${req}" ] || { echo "Missing required release file: ${req}" >&2; exit 1; }
done
[ -x "${REPO_DIR}/tools/METAL/metal" ] || { echo "METAL binary not executable" >&2; exit 1; }

STAGE="$(mktemp -d "${TMPDIR:-/tmp}/medicell_release.XXXXXX")"
trap 'rm -rf -- "${STAGE}"' EXIT
PKG="${STAGE}/MediCell"

mkdir -p "${PKG}/envs" "${PKG}/script" "${PKG}/tools/METAL" "${PKG}/tools/MetaXcan" "${PKG}/tools/R_libs" "${PKG}/data/ref"

cp -a "${REPO_DIR}/medicell.py" "${PKG}/"
chmod +x "${PKG}/medicell.py"
cp -a "${REPO_DIR}/README.md" "${PKG}/"
cp -a "${REPO_DIR}/envs/"*.yml "${PKG}/envs/"

rsync -a \
    --exclude '__pycache__' \
    --exclude '*.pyc' \
    --exclude '*.bak' \
    --exclude 'M00_envs' \
    --exclude 'HARDCODING_NOTES.md' \
    "${REPO_DIR}/script/" "${PKG}/script/"

cp -a "${REPO_DIR}/tools/METAL/metal" "${PKG}/tools/METAL/metal"

cp -a "${REPO_DIR}/tools/MetaXcan/LICENSE" "${REPO_DIR}/tools/MetaXcan/README.md" "${PKG}/tools/MetaXcan/"
rsync -a \
    --exclude 'tests' \
    --exclude 'doc' \
    --exclude '__pycache__' \
    --exclude '*.pyc' \
    --exclude '.gitignore' \
    "${REPO_DIR}/tools/MetaXcan/software/" "${PKG}/tools/MetaXcan/software/"

cp -a "${REPO_DIR}/tools/R_libs/R4.4.3" "${PKG}/tools/R_libs/"
cp -a "${REPO_DIR}/tools/R_libs/r3.6.3" "${PKG}/tools/R_libs/"

cp -a "${REPO_DIR}/data/ref/README.md" "${PKG}/data/ref/"

rm -f "${SOFT_TGZ}" "${DATA_TGZ}"
echo "Creating ${SOFT_TGZ} ..."
tar -czf "${SOFT_TGZ}" -C "${STAGE}" MediCell

echo "Creating ${DATA_TGZ} ..."
tar -czf "${DATA_TGZ}" -C "${REPO_DIR}/data" ref

echo "Done:"
ls -lh "${SOFT_TGZ}" "${DATA_TGZ}"
