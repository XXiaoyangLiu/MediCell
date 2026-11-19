#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

PART1_OUT="${1:-}"
mkdir -p test/results
OUT="$(mktemp -d "$ROOT/test/results/lincs.XXXXXX")"
printf "LINCS test output directory: %s\n" "$OUT"
SECONDS=0

mkdir -p "$OUT/fixture"

# 1. Self-contained multi-gene fixture (runs with or without large Broad GMTs).
python3 - "$OUT/fixture" <<'PY'
from pathlib import Path
import sys

fx = Path(sys.argv[1])
with (fx / "twas.csv").open("w") as f:
    f.write("gene_name,cell_type,zscore,fdr_p\n")
    for i in range(1, 9):
        f.write(f"G{i},CellA,{3.5 - i * 0.2},0.001\n")
        f.write(f"G{i},CellB,{1.5 - i * 0.1},0.01\n")
        f.write(f"N{i},CellA,{-3.5 + i * 0.2},0.001\n")
        f.write(f"N{i},CellB,{-1.5 + i * 0.1},0.01\n")

with (fx / "exp.csv").open("w") as f:
    f.write(",CellA,CellB\n")
    for i in range(1, 9):
        f.write(f"G{i},{10 + i},{2 + i}\n")
        f.write(f"N{i},{2 + i},{10 + i}\n")

with (fx / "custom.gmt").open("w") as f:
    f.write("PERT1\tdrug_one\tG1\tG2\tG3\tG4\tN1\tN2\tN3\tN4\n")
    f.write("BRD-123\tdrug_brd\tG5\tG6\tG7\tG8\tN5\tN6\tN7\tN8\n")
PY

python3 -u medicell.py lincs \
    --output-dir "$OUT/standalone" \
    --cores 1 \
    --twas-file "$OUT/fixture/twas.csv" \
    --exp-matrix "$OUT/fixture/exp.csv" \
    --lincs-gmt "$OUT/fixture/custom.gmt" \
    --lincs-pert-phase1 "" \
    --lincs-pert-phase2 "" \
    --trait TESTTRAIT --label cohortA --direction both \
    --lincs-min-size 3 --lincs-nes-filter all --lincs-sig-cutoff 1.0 \
    --exclude-pert-pattern "" 2>&1 | tee "$OUT/standalone.log"

# 2. If local Broad LINCS GMT + Phase 1/2 metadata exist, verify real GMT/annotation join.
BROAD_GMT="data/ref/lincs/LINCS_PRL_top250_human.gmt"
PERT1="data/ref/lincs/GSE92742_Broad_LINCS_pert_info.txt.gz"
PERT2="data/ref/lincs/GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz"
if [[ -s "$BROAD_GMT" && -s "$PERT1" && -s "$PERT2" ]]; then
    python3 - "$BROAD_GMT" "$OUT/fixture" <<'PY'
from pathlib import Path
import sys

gmt_path = Path(sys.argv[1])
fx = Path(sys.argv[2])
chosen_genes = []
all_genes = []
with gmt_path.open() as f:
    for idx, line in enumerate(f):
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 15:
            continue
        pid, pname, genes = parts[0], parts[1], parts[2:]
        if not pid.startswith("BRD-") and pname != pid and not chosen_genes:
            chosen_genes = genes[:15]
        if idx < 20:
            for g in genes[:15]:
                if g not in all_genes:
                    all_genes.append(g)
        if chosen_genes and len(all_genes) >= 80:
            break

other_genes = [g for g in all_genes if g not in chosen_genes][:60]
ranked = other_genes + chosen_genes
with (fx / "broad_twas.csv").open("w") as f:
    f.write("gene_name,cell_type,zscore,fdr_p\n")
    n = len(ranked)
    for idx, g in enumerate(ranked):
        z = round(5.0 - idx * (4.5 / max(1, n - 1)), 4)
        f.write(f"{g},CellA,{z},0.001\n")
        f.write(f"{g},CellB,{round(z * 0.5, 4)},0.01\n")

with (fx / "broad_exp.csv").open("w") as f:
    f.write(",CellA,CellB\n")
    for idx, g in enumerate(ranked):
        f.write(f"{g},{2.0 + (idx % 5)},{8.0 - (idx % 5)}\n")
PY

    python3 -u medicell.py lincs \
        --output-dir "$OUT/broad" \
        --cores 1 \
        --twas-file "$OUT/fixture/broad_twas.csv" \
        --exp-matrix "$OUT/fixture/broad_exp.csv" \
        --lincs-gmt "$BROAD_GMT" \
        --lincs-pert-phase1 "$PERT1" \
        --lincs-pert-phase2 "$PERT2" \
        --trait BROADTEST --label cohortB --direction pos \
        --lincs-min-size 10 --lincs-nes-filter neg --lincs-sig-cutoff 0.05 \
        --top-n 20 2>&1 | tee "$OUT/broad.log"
fi

# 3. Verify single-cell-type LINCS scoring is explicitly rejected (no silent all-zero scores).
if python3 -u medicell.py lincs \
    --output-dir "$OUT/single_ct_reject" \
    --twas-file "$OUT/fixture/twas.csv" \
    --exp-matrix "$OUT/fixture/exp.csv" \
    --lincs-gmt "$OUT/fixture/custom.gmt" \
    --lincs-pert-phase1 "" --lincs-pert-phase2 "" \
    --cell-types CellA --trait TESTTRAIT --label cohortA --direction pos \
    >"$OUT/single_ct_reject.log" 2>&1; then
    echo "ERROR: single-cell-type LINCS scoring should have failed" >&2
    exit 1
fi

# 4. If TWAS output and Broad GMTs exist, verify full M01 + M08 -> prepLincs -> enrichLincs -> scoreLincs handoff.
if [[ -n "$PART1_OUT" && -s "$BROAD_GMT" && -s "$PERT1" && -s "$PERT2" ]]; then
    python3 -u medicell.py lincs \
        --output-dir "$PART1_OUT" \
        --cores 1 \
        --gene-pos test/data/gene_names.csv \
        --label pbmc --direction neg \
        --twas-padj-cutoff 1.0 --lincs-min-size 1 \
        --lincs-nes-filter all --lincs-sig-cutoff 1.0 2>&1 | tee "$OUT/part1_handoff.log"
fi

python3 test/check_lincs_results.py "$OUT" "${PART1_OUT:-}"
printf "LINCS elapsed seconds: %s\n" "$SECONDS" | tee "$OUT/duration.txt"
