#!/bin/bash
# M04: METAL trans-ethnic meta-analysis (pbmc variant)
# Converts the *_aligned.tsv files into METAL input format and runs METAL.
set -euo pipefail
KIDNEY=0
if [[ "$MC_GWAS_PROFILE" == "sle-legacy-kidney" ]]; then KIDNEY=1; fi
GWAS_OUT="${MC_GWAS_OUT}"
METAL_BIN="${MC_METAL_BIN}"

#GCST011096 EUR+EAS
cat "${GWAS_OUT}/GCST011096_aligned.tsv" | awk '
BEGIN {OFS="\t"}
NR == 1 {print "SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"}
NR > 1 {
    print $9, $4, $5, "NA", $6, $7, $8, "25268"
}' > "${GWAS_OUT}/GCST011096.txt"
#GCST90018917 EUR+EAS
cat "${GWAS_OUT}/GCST90018917_aligned.tsv" | awk '
BEGIN {OFS="\t"}
NR == 1 {print "SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"}
NR > 1 {
    print $9, $3, $4, $5, $6, $7, $8, "659165"
}' > "${GWAS_OUT}/GCST90018917.txt"
#GCST90476183 EUR
cat "${GWAS_OUT}/GCST90476183_aligned.tsv" | awk -v kidney="${KIDNEY}" '
BEGIN {OFS="\t"}
NR == 1 {print "SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"}
NR > 1 {
    se_col = kidney ? $8 : $28
    print $21, $5, $6, $9, $26, se_col, $10, $15
}' > "${GWAS_OUT}/GCST90476183.txt"

#metal
cat > "${GWAS_OUT}/metal_script.txt" << EOF
SCHEME STDERR

MARKER SNP
ALLELE A2 A1
FREQ FREQ
EFFECT BETA
STDERR SE
PVALUE P
SAMPLESIZE N

PROCESS ${GWAS_OUT}/GCST011096.txt
PROCESS ${GWAS_OUT}/GCST90018917.txt
PROCESS ${GWAS_OUT}/GCST90476183.txt

OUTFILE SLE_metal_096_917_183 .tbl
OUTFILE SLE_metal_096_917_183 .tbl
ANALYZE HETEROGENEITY

QUIT
EOF
cd "${GWAS_OUT}"

