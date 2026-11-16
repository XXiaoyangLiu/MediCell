#!/usr/bin/env python3
"""Verify outputs of test/run_lincs.sh."""
import csv
import math
from pathlib import Path
import sys


def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


out = Path(sys.argv[1]).resolve()
part1_out = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 and sys.argv[2] else None

# 1. Verify standalone multi-gene LINCS test
st = out / "standalone"
st_log = (out / "standalone.log").read_text()
assert "MediCell LINCS run complete." in st_log and "Traceback" not in st_log

for rel in [
    "L01_inputs/TESTTRAIT_cohortA_pos.csv",
    "L01_inputs/TESTTRAIT_cohortA_neg.csv",
    "L01_inputs/TESTTRAIT_cohortA_twas_all.tsv",
    "L01_inputs/cohortA_exp_matrix.csv",
    "L02_enrichment/LINCS_all_TESTTRAIT_cohortA_pos.csv",
    "L02_enrichment/LINCS_sig_TESTTRAIT_cohortA_pos.csv",
    "L02_enrichment/LINCS_all_TESTTRAIT_cohortA_neg.csv",
    "L02_enrichment/LINCS_sig_TESTTRAIT_cohortA_neg.csv",
    "L03_drug_scores/cohortA_TESTTRAIT_pos.drug_celltype_score.RData",
    "L03_drug_scores/cohortA_TESTTRAIT_neg.drug_celltype_score.RData",
]:
    p = st / rel
    assert p.is_file() and p.stat().st_size > 0, f"Missing or empty: {p}"

pos_scores = read_csv(st / "L03_drug_scores/cohortA_TESTTRAIT_pos.drug_celltype_score.csv")
neg_scores = read_csv(st / "L03_drug_scores/cohortA_TESTTRAIT_neg.drug_celltype_score.csv")
assert len(pos_scores) == 2 and len(neg_scores) == 2
assert {r["pert_iname"] for r in pos_scores} == {"drug_one", "drug_brd"}

pos_map = {r["pert_iname"]: r for r in pos_scores}
neg_map = {r["pert_iname"]: r for r in neg_scores}
for drug in ("drug_one", "drug_brd"):
    for ct in ("CellA", "CellB"):
        pv = float(pos_map[drug][ct])
        nv = float(neg_map[drug][ct])
        assert math.isfinite(pv) and math.isfinite(nv)
        assert math.isclose(pv, -nv, rel_tol=1e-9), f"Asymmetric score for {drug}/{ct}: {pv} vs {nv}"

# 2. Verify Broad GMT + Phase 1/2 annotation run if executed
broad = out / "broad"
if broad.is_dir():
    broad_log = (out / "broad.log").read_text()
    assert "MediCell LINCS run complete." in broad_log and "Traceback" not in broad_log
    sig_rows = read_csv(broad / "L02_enrichment/LINCS_sig_BROADTEST_cohortB_pos.csv")
    assert len(sig_rows) >= 1, "Expected >=1 significant Broad LINCS hit"
    assert all(r.get("pert_iname") and r.get("pert_type") for r in sig_rows), "Missing pert annotation"
    broad_scores = read_csv(broad / "L03_drug_scores/cohortB_BROADTEST_pos.drug_celltype_score.csv")
    assert len(broad_scores) >= 1, "Expected >=1 scored drug in Broad LINCS test"

# 3. Verify Part 1 M01 + M08 -> Part 2 handoff if executed
if part1_out is not None:
    h_log = (out / "part1_handoff.log").read_text()
    assert "MediCell LINCS run complete." in h_log and "Traceback" not in h_log
    for rel in [
        "L01_inputs/SLE_pbmc_pos.csv",
        "L01_inputs/SLE_pbmc_twas_all.tsv",
        "L01_inputs/pbmc_exp_matrix.csv",
        "L02_enrichment/LINCS_all_SLE_pbmc_pos.csv",
        "L02_enrichment/LINCS_sig_SLE_pbmc_pos.csv",
        "L03_drug_scores/pbmc_SLE_pos.drug_celltype_score.RData",
        "L03_drug_scores/pbmc_SLE_pos.drug_celltype_score.csv",
    ]:
        assert (part1_out / rel).is_file(), f"Missing Part 1->2 handoff output: {rel}"

print("PASS: LINCS standalone, annotation/scoring, and handoff checks succeeded")
print("Output:", out)
