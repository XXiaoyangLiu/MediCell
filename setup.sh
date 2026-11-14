#!/bin/bash
# MediCell one-click environment setup.
# Creates the three conda environments, installs PEER and ENIGMA (R packages),
# clones MetaXcan and builds METAL.
#
# Usage:
#   bash setup.sh [--prefix DIR] [--enigma-src PATH]
#
# Options:
#   --prefix DIR       installation prefix for tools (default: <repo>/tools)
#   --enigma-src PATH  ENIGMA R package source (directory or tarball).
#                      ENIGMA is a lab-internal package and not publicly
#                      downloadable; provide your lab copy here.
set -e

REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
PREFIX="${REPO_DIR}/tools"
ENIGMA_SRC=""

while [ $# -gt 0 ]; do
    case "$1" in
        --prefix) PREFIX="$2"; shift 2 ;;
        --enigma-src) ENIGMA_SRC="$2"; shift 2 ;;
        *) echo "unknown option: $1"; exit 1 ;;
    esac
done

SRC_DIR="${PREFIX}/src"
mkdir -p "${PREFIX}" "${SRC_DIR}" "${PREFIX}/METAL"

# pick a conda new enough to read modern conda-forge repodata
# (the system conda on some clusters is too old, e.g. 4.6 from EPEL)
if [ -n "${CONDA_EXE}" ] && [ -x "${CONDA_EXE}" ]; then
    CONDA="${CONDA_EXE}"
elif [ -x "${HOME}/miniforge3/bin/conda" ]; then
    CONDA="${HOME}/miniforge3/bin/conda"
else
    CONDA="conda"
fi
echo "using conda: ${CONDA}"

echo "== 1/5 conda environments =="
for env_file in envs/R4.4.3.yml envs/r3.6.3.yml envs/imlatools.yml; do
    env_name="$(grep '^name:' "${REPO_DIR}/${env_file}" | cut -d' ' -f2)"
    if "${CONDA}" env list | awk '{print $1}' | grep -qx "${env_name}"; then
        echo "conda env '${env_name}' already exists, skipping"
    else
        "${CONDA}" env create -f "${REPO_DIR}/${env_file}"
    fi
done
# tip: if conda-forge is slow from your location, configure a mirror first,
# e.g. Tsinghua: https://mirrors.tuna.tsinghua.edu.cn/help/anaconda/

# build tools for PEER/METAL compilation (cmake 4 is incompatible with both
# 2016/2020-era CMakeLists). Installed unconditionally so the env's tools win
# over the old system ones (/bin/cmake 2.8, /usr/local/bin/swig).
"${CONDA}" install -y -n r3.6.3 -c https://conda.anaconda.org/conda-forge swig "cmake<4" zlib

echo "== 2/5 PEER (R package, into r3.6.3) =="
if "${CONDA}" run -n r3.6.3 Rscript -e 'library(peer)' >/dev/null 2>&1; then
    echo "PEER already installed in r3.6.3, skipping"
else
    [ -d "${SRC_DIR}/peer" ] || git clone --depth 1 https://github.com/PMBio/peer.git "${SRC_DIR}/peer"
    rm -rf "${SRC_DIR}/peer/build"
    # cmake builds the C++ library and SWIG R bindings; the R package is
    # assembled in build/R/peer and installed automatically (POST_BUILD step)
    "${CONDA}" run -n r3.6.3 bash -c \
        "cd '${SRC_DIR}/peer' && mkdir -p build && cd build && \
         cmake .. -DBUILD_R_PACKAGE=ON -DBUILD_PYTHON_PACKAGE=OFF -DCMAKE_INSTALL_PREFIX='${SRC_DIR}/peer-install' && \
         make"
fi

echo "== 3/5 ENIGMA (R package, into R4.4.3) =="
if "${CONDA}" run -n R4.4.3 Rscript -e 'library(ENIGMA)' >/dev/null 2>&1; then
    echo "ENIGMA already installed in R4.4.3, skipping"
else
    if [ -z "${ENIGMA_SRC}" ]; then
        # public source: https://github.com/WWXkenmo/ENIGMA
        if [ ! -f "${SRC_DIR}/ENIGMA/DESCRIPTION" ]; then
            rm -rf "${SRC_DIR}/ENIGMA"
            # GitHub HTTPS can be flaky; retry once on failure
            git clone --depth 1 https://github.com/WWXkenmo/ENIGMA.git "${SRC_DIR}/ENIGMA" \
                || git clone --depth 1 https://github.com/WWXkenmo/ENIGMA.git "${SRC_DIR}/ENIGMA"
        fi
        ENIGMA_SRC="${SRC_DIR}/ENIGMA"
    fi
    # dependencies via a scratch environment: installing bioconductor packages
    # directly into R4.4.3 is pathologically slow under strict channel
    # priority (r-base resolution conflicts). An empty env with flexible
    # priority solves in minutes; R packages and shared libs are then copied
    # in (same r-base build, binary compatible).
    export CONDA_CHANNEL_PRIORITY=flexible
    "${CONDA}" env remove -n biocdeps -y >/dev/null 2>&1
    "${CONDA}" create -y -n biocdeps \
        -c https://conda.anaconda.org/bioconda \
        -c https://conda.anaconda.org/conda-forge \
        r-base=4.4.3 r-mass r-ggplot2 r-magrittr r-nnls r-plyr r-purrr r-tidyr \
        r-corpcor r-matrixstats r-e1071 r-cowplot r-tibble r-vctrs \
        bioconductor-s4vectors bioconductor-singlecellexperiment \
        bioconductor-summarizedexperiment bioconductor-sva
    ENVS_DIR="$("${CONDA}" info --base)/envs"
    for p in "${ENVS_DIR}/biocdeps/lib/R/library"/*; do
        name="$(basename "${p}")"
        [ -d "${ENVS_DIR}/R4.4.3/lib/R/library/${name}" ] || cp -r "${p}" "${ENVS_DIR}/R4.4.3/lib/R/library/"
    done
    for f in "${ENVS_DIR}/biocdeps/lib"/lib*.so*; do
        base="$(basename "${f}")"
        [ -e "${ENVS_DIR}/R4.4.3/lib/${base}" ] || cp "${f}" "${ENVS_DIR}/R4.4.3/lib/"
    done
    "${CONDA}" env remove -n biocdeps -y >/dev/null 2>&1
    "${CONDA}" run -n R4.4.3 R CMD INSTALL "${ENIGMA_SRC}"
fi

echo "== 4/5 MetaXcan =="
if [ -d "${PREFIX}/MetaXcan" ]; then
    echo "MetaXcan already present, skipping"
else
    git clone --depth 1 https://github.com/hakyimlab/MetaXcan.git "${PREFIX}/MetaXcan"
fi

echo "== 5/5 METAL =="
if [ -x "${PREFIX}/METAL/metal" ]; then
    echo "METAL already built, skipping"
else
    # the original download site (csg.sph.umich.edu) is offline;
    # build from the official GitHub source mirror instead
    if [ ! -f "${SRC_DIR}/Metal/CMakeLists.txt" ]; then
        rm -rf "${SRC_DIR}/Metal"
        # GitHub HTTPS can be flaky; retry once on failure
        git clone --depth 1 https://github.com/statgen/Metal.git "${SRC_DIR}/Metal" \
            || git clone --depth 1 https://github.com/statgen/Metal.git "${SRC_DIR}/Metal"
    fi
    rm -rf "${SRC_DIR}/Metal/build"
    # cmake does not inherit the env's include/lib paths; CPATH/LIBRARY_PATH
    # are honored by gcc directly, so Main.cpp finds zlib.h and libz
    R363_ENV="$("${CONDA}" info --base)/envs/r3.6.3"
    "${CONDA}" run -n r3.6.3 env \
        CPATH="${R363_ENV}/include" \
        LIBRARY_PATH="${R363_ENV}/lib" \
        bash -c "cd '${SRC_DIR}/Metal' && mkdir -p build && cd build && \
                 cmake -DCMAKE_BUILD_TYPE=Release .. && make"
    # binary lands in build/metal/ (build/bin only appears after make install)
    find "${SRC_DIR}/Metal/build" -type f -name metal -perm -u+x \
        -exec cp {} "${PREFIX}/METAL/metal" \;
fi

echo ""
echo "== setup complete =="
echo "medicell.py defaults match this setup:"
echo "  --r-env R4.4.3 --peer-env r3.6.3 --py-env imlatools"
echo "  --metaxcan-dir ${PREFIX}/MetaXcan"
echo "  --metal-bin ${PREFIX}/METAL/metal"
"${CONDA}" run -n R4.4.3 Rscript -e 'library(ENIGMA)' >/dev/null 2>&1 \
    && echo "  ENIGMA: OK" || echo "  ENIGMA: MISSING (run with --enigma-src)"
"${CONDA}" run -n r3.6.3 Rscript -e 'library(peer)' >/dev/null 2>&1 \
    && echo "  PEER: OK" || echo "  PEER: MISSING"
"${CONDA}" run -n imlatools python -c 'import pandas, numpy, scipy' >/dev/null 2>&1 \
    && echo "  MetaXcan runtime: OK" || echo "  MetaXcan runtime: MISSING"
[ -x "${PREFIX}/METAL/metal" ] && echo "  METAL: OK" || echo "  METAL: MISSING"
