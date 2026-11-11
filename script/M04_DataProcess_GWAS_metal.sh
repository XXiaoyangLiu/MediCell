#!/bin/bash
cd "/node1/liuxy/MediCell/pbmc_2/TWAS/usefile/GWAS/"
#GCST011096 EUR+EAS
cat "/node1/liuxy/MediCell/pbmc_2/TWAS/usefile/GWAS/GCST011096_aligned.tsv" | awk '
BEGIN {OFS="\t"}
NR == 1 {print "SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"}
NR > 1 {
    print $9, $4, $5, "NA", $6, $7, $8, "25268"
}' > /node1/liuxy/MediCell/pbmc_2/TWAS/usefile/GWAS/GCST011096.txt
#GCST90018917 EUR+EAS
cat "/node1/liuxy/MediCell/pbmc_2/TWAS/usefile/GWAS/GCST90018917_aligned.tsv" | awk '
BEGIN {OFS="\t"}
NR == 1 {print "SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"}
NR > 1 {
    print $9, $3, $4, $5, $6, $7, $8, "659165"
}' > /node1/liuxy/MediCell/pbmc_2/TWAS/usefile/GWAS/GCST90018917.txt
#GCST90476183 EUR
cat "/node1/liuxy/MediCell/pbmc_2/TWAS/usefile/GWAS/GCST90476183_aligned.tsv" | awk '
BEGIN {OFS="\t"}
NR == 1 {print "SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"}
NR > 1 {
    print $21, $5, $6, $9, $26, $28, $10, $15
}' > /node1/liuxy/MediCell/pbmc_2/TWAS/usefile/GWAS/GCST90476183.txt

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

PROCESS /node1/liuxy/MediCell/pbmc_2/TWAS/usefile/GWAS/GCST011096.txt
PROCESS /node1/liuxy/MediCell/pbmc_2/TWAS/usefile/GWAS/GCST90018917.txt
PROCESS /node1/liuxy/MediCell/pbmc_2/TWAS/usefile/GWAS/GCST90476183.txt

OUTFILE SLE_metal_096_917_183 .tbl
OUTFILE SLE_metal_096_917_183 .tbl
ANALYZE HETEROGENEITY

QUIT
EOF
"/public/home/liuxy/METAL-2018-08-28/build/bin/metal" metal_script.txt
