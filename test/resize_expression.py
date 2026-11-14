#!/usr/bin/env python3
"""Resize only the expression panel; preserve samples, genotypes and GWAS."""
import csv
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import sys

root = Path(__file__).resolve().parents[1]
data = root / 'test/data'
manifest = json.loads((data / 'manifest.json').read_text())
count = int(sys.argv[1])
old_count = manifest['bulk_reference_genes']
backup = root / 'test/attempts' / ('genes' + str(old_count))
assert not backup.exists(), 'Refusing to overwrite a previous fixture snapshot'
shutil.copytree(str(data), str(backup))
samples = (data / 'discovery_samples.txt').read_text().splitlines()
with open(manifest['sources']['reference']['path']) as f:
    cell_types = shlex.split(f.readline())
    reference = {a[0]: [float(v) for v in a[1:]] for a in (shlex.split(line) for line in f) if a}
genes = [manifest['target_gene']]
ranks = [sorted(reference, key=lambda g: (-reference[g][i] / (sum(reference[g]) + 1e-12),
                                         -reference[g][i], g)) for i in range(len(cell_types))]
for rank in range(len(reference)):
    for i in range(len(cell_types)):
        gene = ranks[i][rank]
        if gene not in genes and min(reference[gene]) > 0.1:
            genes.append(gene)
        if len(genes) == count:
            break
    if len(genes) == count:
        break
assert len(genes) == count
bulk = {}
with open(manifest['sources']['bulk']['path']) as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    idx = [header.index(s) for s in samples]
    for row in reader:
        if row[0] in genes:
            bulk[row[0]] = [row[i] for i in idx]
assert set(bulk) == set(genes)
for name, header, values in [('bulk_tpm.tsv', ['gene_id'] + samples, bulk),
                              ('reference_tpm.tsv', ['gene_id'] + cell_types, reference)]:
    with (data / name).open('w') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows([[g] + values[g] for g in genes])
manifest['bulk_reference_genes'] = count
manifest['files'] = {p.name: {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
                     for p in data.iterdir() if p.name != 'manifest.json'}
(data / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
print('Expression genes:', count, 'Discovery samples:', len(samples))
