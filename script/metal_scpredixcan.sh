#!/bin/bash
cd "/node1/liuxy/MediCell/method_compare_2/usefile/GWAS_scpredixcan/"

#GCST011096
cat GCST011096_aligned.tsv | awk '
BEGIN {OFS="\t"}
NR == 1 {print "SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"}
NR > 1 {
    print $11, $6, $7, "NA", $8, $9, $10, "25268"
}' > GCST011096.txt
#GCST90018917
cat GCST90018917_aligned.tsv | awk '
BEGIN {OFS="\t"}
NR == 1 {print "SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"}
NR > 1 {
    print $11, $5, $6, $7, $8, $9, $10, "659165"
}' > GCST90018917.txt
#GCST90476183
cat GCST90476183_aligned.tsv | awk '
BEGIN {OFS="\t"}
NR == 1 {print "SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"}
NR > 1 { 
    print $21, $3, $4, $7, $22, $6, $8, $15
}' > GCST90476183.txt

#metal
cat > metal_script.txt << 'EOF'
SCHEME STDERR

MARKER SNP
ALLELE A2 A1
FREQ FREQ
EFFECT BETA
STDERR SE
PVALUE P
SAMPLESIZE N

PROCESS GCST011096.txt
PROCESS GCST90018917.txt
PROCESS GCST90476183.txt

OUTFILE SLE_metal_096_917_183 .tbl
OUTFILE SLE_metal_096_917_183 .tbl
ANALYZE HETEROGENEITY

QUIT
EOF
"/public/home/liuxy/METAL-2018-08-28/build/bin/metal" metal_script.txt
