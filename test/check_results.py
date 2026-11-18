#!/usr/bin/env python3
"""Check outputs of the chr10 de-identified end-to-end test."""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import sys

root = Path(__file__).resolve().parents[1]
out = Path(sys.argv[1]).resolve()
data = root / "test/data"
manifest = json.loads((data / "manifest.json").read_text())
for name, info in manifest["files"].items():
    assert hashlib.sha256((data / name).read_bytes()).hexdigest() == info["sha256"], name


def table(path, delimiter="\t"):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter=delimiter))


log = (out / "run.log").read_text()
assert "MediCell run complete." in log and "Traceback" not in log

n_samples = int(manifest["discovery_sample_count"])
n_snps = int(manifest["training_variant_count"])
cell_types = manifest["analysis_cell_types"]
assert log.count("[run]") == 11 + len(cell_types), "Not all twelve child steps ran"

counts = table(out / "M04_gwas/gwas_normalization_counts.tsv")
assert {r["study_id"] for r in counts} == {"bentham", "jiang"}
assert all(int(r["matched"]) == n_snps for r in counts)

twas_rows_total = 0
for ct in cell_types:
    summary = table(out / f"M03_models/summary/Model_training_{ct}_chr10_summary.txt")
    assert len(summary) == 1 and int(summary[0]["n_samples"]) == n_samples
    models = table(out / f"M03_models/summary/Model_training_{ct}_chr10_model_summaries.txt")
    assert len(models) >= 1 and all(len(m) == 24 for m in models)
    weights = table(out / f"M03_models/weights/Model_training_{ct}_chr10_weights.txt")
    assert len(weights) >= 1 and all(math.isfinite(float(r["beta"])) for r in weights)
    with sqlite3.connect(str(out / f"M06_databases/{ct}_models_unfiltered.db")) as c:
        assert c.execute("SELECT n_samples,population FROM sample_info").fetchall() == [(n_samples, "EUR")]
        assert c.execute("SELECT count(*) FROM extra").fetchone()[0] >= 1
        assert c.execute("SELECT count(*) FROM weights").fetchone()[0] >= 1
    with sqlite3.connect(str(out / f"M06_databases/{ct}_models_filtered_signif.db")) as cf:
        assert cf.execute("SELECT n_samples,population FROM sample_info").fetchall() == [(n_samples, "EUR")]
        assert cf.execute("SELECT count(*) FROM extra WHERE \"pred.perf.R2\" <= 0.01 OR \"pred.perf.pval\" >= 0.05").fetchone()[0] == 0
    with gzip.open(str(out / f"M05_ld/cov_MINI_EUR_EAS/{ct}.txt.gz"), "rt") as f:
        header = f.readline().split()
        covariance = [dict(zip(header, line.split())) for line in f]
    assert len(covariance) >= 1
    assert all(math.isfinite(float(r["VALUE"])) for r in covariance)
    twas = table(out / f"M07_twas/{ct}.csv", ",")
    assert len(twas) >= 1
    twas_rows_total += len(twas)
    for r in twas:
        assert r["gene"] in manifest["target_genes"]
        assert int(float(r["n_snps_used"])) >= 1
        pv = float(r["pvalue"])
        assert math.isclose(pv, math.erfc(abs(float(r["zscore"])) / math.sqrt(2)), rel_tol=1e-10)

all_rows = table(out / "M08_summary/allcelltype_result.csv", ",")
assert len(all_rows) == twas_rows_total
assert {r["cell_type"] for r in all_rows} == set(cell_types)
for r in all_rows:
    assert r["gene_id"] in manifest["target_genes"]
    assert r["gene_name"] == manifest["target_genes"][r["gene_id"]]
    assert math.isfinite(float(r["fdr_p"])) and math.isfinite(float(r["bonferroni_p"]))
for name in ["Bonfsig_allcelltype_result.csv", "FDRsig_allcelltype_result.csv", "BonFDRsig_allcelltype_result.csv"]:
    assert (out / "M08_summary" / name).is_file()

print(f"PASS: all 12 steps; {n_samples} samples; chr10; cell_types={cell_types}; {twas_rows_total} TWAS results; {n_snps} SNPs; 2 GWAS studies")
print("Output:", out)
