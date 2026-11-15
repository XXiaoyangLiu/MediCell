#!/usr/bin/env python3
"""Check the actual outputs of the complete private smoke test."""
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
data = root / 'test/data'
manifest = json.loads((data / 'manifest.json').read_text())
for name, info in manifest['files'].items():
    assert hashlib.sha256((data / name).read_bytes()).hexdigest() == info['sha256'], name


def table(path, delimiter='\t'):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter=delimiter))


log = (out / 'run.log').read_text()
assert log.count('[run]') == 12, 'Not all twelve child steps ran'
assert 'MediCell run complete.' in log and 'Traceback' not in log
summary = table(out / 'M03_models/summary/Model_training_CD4T_chr22_summary.txt')
assert len(summary) == 1 and int(summary[0]['n_samples']) == 40
models = table(out / 'M03_models/summary/Model_training_CD4T_chr22_model_summaries.txt')
assert len(models) == 1 and len(models[0]) == 24
weights = table(out / 'M03_models/weights/Model_training_CD4T_chr22_weights.txt')
assert len(weights) == 8 and all(math.isfinite(float(r['beta'])) for r in weights)
with sqlite3.connect(str(out / 'M06_databases/CD4T_models_unfiltered.db')) as c:
    assert c.execute('SELECT n_samples,population FROM sample_info').fetchall() == [(40, 'EUR')]
    assert c.execute('SELECT count(*) FROM extra').fetchone()[0] == 1
    assert c.execute('SELECT count(*) FROM weights').fetchone()[0] == 8
counts = table(out / 'M04_gwas/gwas_normalization_counts.tsv')
assert {r['study_id'] for r in counts} == {'bentham', 'jiang'}
assert all(int(r['matched']) == 8 for r in counts)
with gzip.open(str(out / 'M05_ld/cov_MINI_EUR_EAS/CD4T.txt.gz'), 'rt') as f:
    header = f.readline().split()
    covariance = [dict(zip(header, line.split())) for line in f]
assert len(covariance) == 36
assert all(math.isfinite(float(r['VALUE'])) for r in covariance)
twas = table(out / 'M07_twas/CD4T.csv', ',')
assert len(twas) == 1 and twas[0]['gene'] == manifest['target_gene']
r = twas[0]
assert all(int(float(r[k])) == 8 for k in ['n_snps_used', 'n_snps_in_cov', 'n_snps_in_model'])
pvalue = float(r['pvalue'])
assert math.isclose(pvalue, math.erfc(abs(float(r['zscore'])) / math.sqrt(2)), rel_tol=1e-10)
for name in ['allcelltype_result.csv', 'Bonfsig_allcelltype_result.csv',
             'FDRsig_allcelltype_result.csv', 'BonFDRsig_allcelltype_result.csv']:
    rows = table(out / 'M08_summary' / name, ',')
    assert len(rows) == 1
    assert rows[0]['gene_id'] == manifest['target_gene']
    assert rows[0]['gene_name'] == manifest['target_symbol']
    assert math.isclose(float(rows[0]['fdr_p']), pvalue, rel_tol=1e-12)
    assert math.isclose(float(rows[0]['bonferroni_p']), pvalue, rel_tol=1e-12)
print('PASS: all 12 steps; 40 discovery samples; 1 new model; 8 weights; 2 GWAS studies; LD and TWAS; all 4 summaries')
print('Output:', out)
print('TWAS zscore:', r['zscore'], 'pvalue:', r['pvalue'])
