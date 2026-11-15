#!/usr/bin/env python3
"""Extract a private real-data smoke fixture without changing pipeline algorithms."""
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import sqlite3
import subprocess

os.umask(0o077)
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'test/data'
BCF = '/public/home/dongb/miniforge3/envs/imlatools/bin/bcftools'
SOURCES = {
    'bulk': '/node1/liuxy/chplot/combind_plasma_B/enigma/usefile/pbmc_6_CellType_971samples_pseudobulk_TPM.txt',
    'reference': '/node1/liuxy/chplot/combind_plasma_B/enigma/usefile/TPM_ENIGMA_ref_6CT.matrx',
    'fractions': '/node1/liuxy/chplot/cibersort/usefile/pbmc/_7celltype/CIBERSORTx_Job129_Adjusted_pbmc_54genes.txt',
    'covariates': '/node1/liuxy/sc_eQTL_science2022/OneK1K_sex_age_10PCs.txt',
    'training_vcf': '/node1/liuxy/sc_eQTL_science2022/filter_vcf_r08/chr22.dose.filtered.R2_0.8.vcf.gz',
    'ld_vcf': '/node1/liuxy/index/1000G_2/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz',
    'eur_samples': '/node1/liuxy/index/1000G/EUR_samples.txt',
    'eur_eas_samples': '/node1/liuxy/index/1000G/EUR_EAS_samples.txt',
    'gene_annot': '/node1/liuxy/MediCell/pbmc/predictdb/usefile/gene/gene_annot.parsed.txt',
    'gene_names': '/node1/liuxy/index/gene_pos_grch37.82.csv',
    'bentham': '/node1/zhaoq/scTWAS_PBMC/SPrediXcan/input/GWAS/bentham_2015NG/bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv.gz',
    'jiang': '/node1/zhaoq/scTWAS_PBMC/SPrediXcan/input/GWAS/Jiang_2021NG/GCST90044512_buildGRCh37.tsv.gz',
}
GENE = 'ENSG00000100376'
EXPRESSION_GENES = 128
assert DATA.is_dir() and not any(DATA.iterdir()), 'Refusing to overwrite test/data'
DATA.chmod(0o700)


def write_table(name, header, rows, delimiter='\t'):
    with (DATA / name).open('w') as f:
        writer = csv.writer(f, delimiter=delimiter, lineterminator='\n')
        writer.writerow(header)
        writer.writerows(rows)


def query(vcf, sample_file, fmt):
    result = subprocess.check_output([BCF, 'query', '-R', str(DATA / 'candidate_sites.tsv'),
                                      '-S', str(DATA / sample_file), '-f', fmt, vcf])
    return [line.split('\t') for line in result.decode().splitlines()]


with open(SOURCES['bulk']) as f:
    bulk_ids = next(csv.reader(f, delimiter='\t'))[1:]
with open(SOURCES['covariates']) as f:
    cov = list(csv.reader(f, delimiter='\t'))
with open(SOURCES['fractions']) as f:
    fractions = list(csv.DictReader(f, delimiter='\t'))
fraction_map = {row['Mixture']: row for row in fractions}
vcf_ids = subprocess.check_output([BCF, 'query', '-l', SOURCES['training_vcf']]).decode().splitlines()
shared = [s for s in bulk_ids if 'X' + s in fraction_map and s in cov[0][1:] and s in vcf_ids]
samples = shared[:40]
assert len(samples) == len(set(samples)) == 40
(DATA / 'discovery_samples.txt').write_text('\n'.join(samples) + '\n')
write_table('covariates.tsv', ['sample_id'] + [r[0] for r in cov[1:]],
            [[s] + [r[cov[0].index(s)] for r in cov[1:]] for s in samples])

with open(SOURCES['reference']) as f:
    cell_types = shlex.split(f.readline())
    reference = {a[0]: [float(v) for v in a[1:]] for a in (shlex.split(line) for line in f) if a}
assert GENE in reference and 'CD4T' in cell_types
selected_genes = [GENE]
ranks = [sorted(reference, key=lambda g: (-reference[g][i] / (sum(reference[g]) + 1e-12),
                                         -reference[g][i], g)) for i in range(len(cell_types))]
for rank in range(len(reference)):
    for i in range(len(cell_types)):
        gene = ranks[i][rank]
        if gene not in selected_genes and min(reference[gene]) > 0.1:
            selected_genes.append(gene)
        if len(selected_genes) == EXPRESSION_GENES:
            break
    if len(selected_genes) == EXPRESSION_GENES:
        break
assert len(selected_genes) == EXPRESSION_GENES
bulk_rows = {}
with open(SOURCES['bulk']) as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    idx = [header.index(s) for s in samples]
    for row in reader:
        if row[0] in selected_genes:
            values = [row[i] for i in idx]
            assert all(math.isfinite(float(v)) and float(v) >= 0 for v in values)
            bulk_rows[row[0]] = values
assert set(bulk_rows) == set(selected_genes)
write_table('bulk_tpm.tsv', ['gene_id'] + samples, [[g] + bulk_rows[g] for g in selected_genes])
write_table('reference_tpm.tsv', ['gene_id'] + cell_types, [[g] + reference[g] for g in selected_genes])
fraction_rows = []
for s in samples:
    row = fraction_map['X' + s]
    values = [float(row[ct]) + (float(row['Plasma']) if ct == 'B' else 0) for ct in cell_types]
    assert abs(sum(values) - 1) < 1e-6
    fraction_rows.append([s] + values)
write_table('fractions.tsv', ['Mixture'] + cell_types, fraction_rows)

with open(SOURCES['gene_annot']) as f:
    annot = [r for r in csv.DictReader(f, delimiter='\t') if r['gene_name'] == GENE]
assert len(annot) == 1 and annot[0]['chr'] == '22'
a = annot[0]
write_table('gene_annot.tsv', ['chr', 'gene_name', 'start', 'end', 'gene_type', 'gene_id'],
            [[a['chr'], GENE, a['start'], a['end'], a['gene_type'], GENE]])
with open(SOURCES['gene_names']) as f:
    names = [r for r in csv.DictReader(f) if r['gene_id'] == GENE]
assert len(names) == 1
write_table('gene_names.csv', ['gene_id', 'gene_name'], [[GENE, names[0]['gene_name']]], ',')

with sqlite3.connect(str(ROOT / 'data/CD4T_PCs_and_age_gender_models_filtered_signif.db')) as c:
    candidates = c.execute('SELECT varID,weight FROM weights WHERE gene=? ORDER BY abs(weight) DESC,varID', (GENE,)).fetchall()
variant_map = {v.split('_')[1]: v for v, w in candidates if v.split('_')[2] + v.split('_')[3] not in ['AT', 'TA', 'CG', 'GC']}
(DATA / 'candidate_sites.tsv').write_text(''.join('22\t' + p + '\n' for p in sorted(variant_map, key=int)))
train_rows = query(SOURCES['training_vcf'], 'discovery_samples.txt', '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n')
train = {'_'.join(r[:4]): r[4:] for r in train_rows if len(r) == 44 and all(v != '.' for v in r[4:])}
eur = Path(SOURCES['eur_samples']).read_text().splitlines()
eur_eas = Path(SOURCES['eur_eas_samples']).read_text().splitlines()
eas = [s for s in eur_eas if s not in set(eur)]
ld_samples = eur[:10] + eas[:10]
(DATA / 'ld_samples.txt').write_text('\n'.join(ld_samples) + '\n')
ld_rows = query(SOURCES['ld_vcf'], 'ld_samples.txt', '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n')
ld = {}
for r in ld_rows:
    if len(r) != 24 or any(gt not in ['0|0', '0|1', '1|0', '1|1', '0/0', '0/1', '1/0', '1/1'] for gt in r[4:]):
        continue
    doses = [int(gt[0]) + int(gt[2]) for gt in r[4:]]
    if len(set(doses)) > 1:
        ld['_'.join(r[:4])] = r

studies = {}
configs = [('bentham', ['chrom', 'pos', 'effect_allele', 'other_allele', 'beta', 'se', 'p']),
           ('jiang', ['chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'beta', 'standard_error', 'p_value'])]
for study, columns in configs:
    print('Scanning real GWAS:', study, flush=True)
    found = {}
    with gzip.open(SOURCES[study], 'rt') as f:
        header = f.readline().rstrip('\n').split('\t')
        idx = [header.index(k) for k in columns]
        for line in f:
            if not line.startswith('22\t'):
                continue
            row = line.rstrip('\n').split('\t')
            pos = row[idx[1]]
            if pos not in variant_map:
                continue
            r = [row[i] for i in idx]
            variant = variant_map[pos]
            ref, alt = variant.split('_')[2:]
            try:
                valid = math.isfinite(float(r[4])) and math.isfinite(float(r[5])) and float(r[5]) > 0 and 0 <= float(r[6]) <= 1
            except ValueError:
                valid = False
            if valid and {r[2].upper(), r[3].upper()} == {ref, alt}:
                assert variant not in found, 'Duplicate source GWAS locus'
                found[variant] = r
    studies[study] = found
    print(study, 'usable candidate variants:', len(found), flush=True)
usable = [v for v, w in candidates if v in train and v in ld and all(v in studies[s] for s in studies)
          and all(math.isfinite(float(x)) and 0 <= float(x) <= 2 for x in train[v]) and len(set(train[v])) > 1]
variants = usable[:8]
assert len(variants) >= 2, 'Not enough complete real variants; no synthetic replacement allowed'
variants.sort(key=lambda v: int(v.split('_')[1]))
write_table('genotype.chr22.txt', ['varID'] + samples, [[v] + train[v] for v in variants])
snp_rows = [[v.split('_')[0], v.split('_')[1], v, v.split('_')[2], v.split('_')[3], v] for v in variants]
write_table('snp_annot.chr22.txt', ['chromosome', 'pos', 'varID', 'ref_vcf', 'alt_vcf', 'rsid'], snp_rows)
write_table('gwas_panel.tsv', ['chromosome', 'Pos', 'VariantID', 'ref_vcf', 'alt_vcf', 'RSID'], snp_rows)
(DATA / 'sites.tsv').write_text(''.join('22\t' + v.split('_')[1] + '\n' for v in variants))
subprocess.run([BCF, 'view', '-R', str(DATA / 'sites.tsv'), '-S', str(DATA / 'ld_samples.txt'),
                '-m2', '-M2', '-v', 'snps', '-Oz', '-o', str(DATA / 'reference.chr22.vcf.gz'), SOURCES['ld_vcf']], check=True)
subprocess.run([BCF, 'index', '-t', str(DATA / 'reference.chr22.vcf.gz')], check=True)
actual_ld_samples = subprocess.check_output([BCF, 'query', '-l', str(DATA / 'reference.chr22.vcf.gz')]).decode().splitlines()
assert set(actual_ld_samples) == set(ld_samples)
(DATA / 'ld_samples.txt').write_text('\n'.join(actual_ld_samples) + '\n')
manifest_rows = []
for study in studies:
    name = 'gwas_' + study + '.tsv'
    write_table(name, ['chr', 'pos', 'ea', 'oa', 'beta', 'se', 'p'], [studies[study][v] for v in variants])
    manifest_rows.append([study, str(DATA / name), 'GRCh37', 'log_odds', 'chr', 'pos', 'ea', 'oa',
                          'beta_se', 'beta', 'se', 'p', '', 'none', '', ''])
write_table('gwas_manifest.tsv', 'study_id path build effect_scale chr_col pos_col ea_col oa_col effect_mode effect_col se_col p_col eaf_col n_mode n_value n_col'.split(), manifest_rows)
for name in ['candidate_sites.tsv', 'sites.tsv']:
    (DATA / name).unlink()
manifest = {
    'purpose': 'Private real-data execution smoke test, not inference or historical reproduction',
    'discovery_sample_count': 40, 'discovery_selection': 'First 40 shared IDs in source bulk order; no phenotype or GWAS selection',
    'bulk_reference_genes': EXPRESSION_GENES, 'reference_cell_types': cell_types, 'training_cell_type': 'CD4T',
    'target_gene': GENE, 'target_symbol': names[0]['gene_name'], 'training_variant_count': len(variants),
    'variant_selection': 'Up to 8 largest-absolute-weight nonpalindromic SNPs from an existing CD4T model, requiring complete real training, LD and two-study GWAS coverage; not an unbiased validation dataset',
    'ld_sample_count': 20, 'ld_composition': '10 EUR + 10 EAS; deliberately tiny smoke-test panel',
    'fraction_transformation': 'Map X-prefixed IDs to uniquely matching bulk IDs; merge original Plasma fraction into B to match six-type reference',
    'gene_annotation': 'Preserve original TSS point start/end; add gene_id equal to original ENSG gene_name explicitly',
    'sources': {k: {'path': v, 'bytes': Path(v).stat().st_size, 'mtime_ns': Path(v).stat().st_mtime_ns} for k, v in SOURCES.items()},
    'files': {p.name: {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in DATA.iterdir() if p.is_file()},
}
(DATA / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
print('Prepared:', DATA, 'samples=40, genes=' + str(EXPRESSION_GENES) + ', targets=1, variants=' + str(len(variants)), flush=True)
print('Total bytes:', sum(p.stat().st_size for p in DATA.iterdir() if p.is_file()), flush=True)
