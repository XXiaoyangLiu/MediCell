"""Validate the project-owned dosage boundary without changing vendor math."""
import csv
import gzip
import math
from pathlib import Path
import re
import sys


def alleles(ref, alt):
    if ref not in {'A', 'C', 'G', 'T'} or alt not in {'A', 'C', 'G', 'T'} or ref == alt:
        raise ValueError('Only biallelic SNPs are supported')


def finite(value, low, high):
    x = float(value)
    if not math.isfinite(x) or not low <= x <= high:
        raise ValueError(f'Invalid finite value in [{low},{high}]: {value}')
    return x


def dosage_rows(path, chrom, positions=None):
    seen = set()
    n = None
    with gzip.open(path, 'rt') as f:
        for line in f:
            a = line.rstrip('\n').split('\t')
            if len(a) < 8 or (n is not None and len(a) != n):
                raise ValueError(f'Invalid/inconsistent dosage sample count: {path}')
            n = len(a)
            if positions is not None and a[2] not in positions:
                continue
            c, key, pos, ref, alt, af = a[:6]
            alleles(ref, alt)
            if c != chrom or not re.fullmatch('[1-9][0-9]*', pos) or key != f'{c}_{pos}_{ref}_{alt}':
                raise ValueError(f'Invalid chromosome/VarID: {key}')
            if (c, pos) in seen:
                raise ValueError(f'Duplicate/multiallelic dosage position: {key}')
            seen.add((c, pos))
            finite(af, 0, 1)
            for x in a[6:]:
                finite(x, 0, 2)
            yield a
    if not seen:
        raise ValueError(f'Empty dosage: {path}')


def read_reference(annot, chrom):
    reference = {}
    with open(annot) as f:
        reader = csv.DictReader(f, delimiter='\t')
        if reader.fieldnames != ['chromosome', 'pos', 'varID', 'ref_vcf', 'alt_vcf', 'rsid']:
            raise ValueError('Unknown SNP annotation schema')
        for a in reader:
            c, pos, key, ref, alt = (a[k] for k in ['chromosome', 'pos', 'varID', 'ref_vcf', 'alt_vcf'])
            alleles(ref, alt)
            if c != chrom or not re.fullmatch('[1-9][0-9]*', pos) or key != f'{c}_{pos}_{ref}_{alt}':
                raise ValueError(f'Invalid annotation VarID: {key}')
            if pos in reference:
                raise ValueError(f'Duplicate/multiallelic reference position: {c}:{pos}')
            reference[pos] = (key, ref, alt)
    if not reference:
        raise ValueError('Empty annotation')
    return reference


def sites(annot, chrom, output):
    reference = read_reference(annot, chrom)
    with open(output, 'w') as f:
        for pos in sorted(reference, key=int):
            f.write(chrom + '\t' + pos + '\n')


def align(annot, source, output, chrom):
    reference = read_reference(annot, chrom)
    count = 0
    with gzip.open(output, 'wt') as out:
        for a in dosage_rows(source, chrom, reference):
            if a[2] not in reference:
                continue
            key, ref, alt = reference[a[2]]
            if (a[3], a[4]) == (ref, alt):
                pass
            elif (a[3], a[4]) == (alt, ref):
                a[5] = str(1-float(a[5]))
                a[6:] = [format(2-float(x), '.15g') for x in a[6:]]
            else:
                # Historical matching skips unsupported allele pairs, never guesses strand.
                continue
            a[1], a[3], a[4] = key, ref, alt
            out.write('\t'.join(a)+'\n')
            count += 1
    if not count:
        raise ValueError('No aligned dosage records')
    print(f'Aligned {count} SNPs; AF computed before sample subset (not re-estimated)')


def samples(requested, available, output):
    wanted = Path(requested).read_text().splitlines()
    have = Path(available).read_text().splitlines()
    if not wanted or any(not x or len(x.split())!=1 for x in wanted) or len(set(wanted)) != len(wanted):
        raise ValueError('Empty/duplicate/invalid sample IDs')
    if len(set(have)) != len(have) or not set(wanted) <= set(have):
        raise ValueError('Missing/duplicate VCF samples')
    if len(wanted) < 2:
        raise ValueError('At least two LD samples required')
    Path(output).write_text('\n'.join(wanted)+'\n')


def gt(output, chrom, count):
    seen = set()
    with gzip.open(output, 'wt') as out:
        for line in sys.stdin:
            a = line.rstrip('\n').split('\t')
            if len(a) != 5+int(count):
                raise ValueError('GT sample count mismatch')
            c, pos, ref, alt, af = a[:5]
            alleles(ref, alt)
            if c != chrom or not re.fullmatch('[1-9][0-9]*', pos) or pos in seen:
                raise ValueError('Invalid/duplicate/multiallelic VCF position')
            seen.add(pos)
            finite(af, 0, 1)
            doses = []
            for value in a[5:]:
                if not re.fullmatch('[01][/|][01]', value):
                    raise ValueError(f'Missing/non-diploid/non-biallelic GT: {value}')
                doses.append(str(int(value[0]!='0')+int(value[2]!='0')))
            out.write('\t'.join([c, f'{c}_{pos}_{ref}_{alt}', pos, ref, alt, af]+doses)+'\n')
    if not seen:
        raise ValueError('Empty VCF dosage extraction')


def check_many(folder, *chroms):
    sizes = set()
    for chrom in chroms:
        for a in dosage_rows(Path(folder)/f'chr{chrom}.txt.gz', chrom):
            sizes.add(len(a)-6)
    if len(sizes)!=1:
        raise ValueError('Inconsistent sample counts across chromosomes')


def covariance(path):
    count = 0
    with gzip.open(path, 'rt') as f:
        r = csv.DictReader(f, delimiter='\t')
        if r.fieldnames != ['GENE','RSID1','RSID2','VALUE']:
            raise ValueError('Invalid covariance header')
        for a in r:
            if None in a or any(not a[k] for k in r.fieldnames):
                raise ValueError('Invalid covariance row')
            if a['VALUE'] != 'NA':
                finite(a['VALUE'], -float('inf'), float('inf'))
                count += 1
    if not count:
        raise ValueError('Missing/empty/all-NA covariance output')


if __name__ == '__main__':
    commands = {'sites': sites, 'align': align, 'samples': samples, 'gt': gt, 'check-many': check_many, 'covariance': covariance}
    try:
        commands[sys.argv[1]](*sys.argv[2:])
    except (ValueError, OSError, KeyError, TypeError) as e:
        sys.exit(f'LD contract error: {e}')
