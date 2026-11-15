#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec "${MC_PYTHON:-python3}" "$HERE/helpers/gwas_control.py" meta
