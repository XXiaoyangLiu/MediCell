## SLE bentham_2015NG
## hg19
## calculate zscore, and change ref alt beta with db
library(data.table)
library(stringr)
a=fread('/node1/zhaoq/scTWAS_PBMC/SPrediXcan/input/GWAS/bentham_2015NG/bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv.gz')
db_an=fread("/node1/liuxy/chplot/combind_plasma_B/genotype/usefile/snp_annotation.txt")
a <- a[a$chrom %in% c(1:22)]
a$chrom=as.numeric(a$chrom)
db_an$Pos=as.numeric(db_an$Pos)
#db_an$chromosome=as.numeric(db_an$chromosome)

merged <- merge(a, db_an,
                #by.x = c("#CHROM", "POS"),   
                by.x = c("chrom", "pos"),
                by.y = c("chromosome", "Pos"),   
                all = FALSE)  #5210246 SNPs

condition <- merged$ref_vcf == merged$effect_allele & merged$alt_vcf == merged$other_allele
str(which(condition))
merged$other_allele[condition]=merged$ref_vcf[condition]
merged$effect_allele[condition]=merged$alt_vcf[condition]
merged$beta[which(condition)]=-merged$beta[which(condition)]
sig_a_alle=merged[,-c('ref_vcf','alt_vcf','snp_id_originalVCF')]
sig_a_alle$zscore=sig_a_alle$beta/sig_a_alle$se
sig_a_alle_clean <- sig_a_alle[!is.na(sig_a_alle$zscore) & sig_a_alle$beta != 0 & sig_a_alle$se != 0, ]


fwrite(sig_a_alle_clean, file = "/node1/liuxy/MediCell/pbmc/TWAS/usefile/bentham_2015_26502338_sle_processed_beta.tsv",sep='\t')


