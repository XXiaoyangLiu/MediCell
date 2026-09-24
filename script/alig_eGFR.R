##eGFR Karczewski
##算zscore，db和dosage等位基因顺序相反的SNP beta取负值
library(data.table)
library(stringr)
a=fread("/node1/liuxy/add_plot/allimmu_TWAS/new/usefile/eGFR_Karczewski/GCST90691925.tsv.gz")
db_an=fread("/node1/liuxy/chplot/combind_plasma_B/genotype/usefile/snp_annotation.txt")
colnames(a)[colnames(a) == "chromosome"] <- "chrom"
colnames(a)[colnames(a) == "base_pair_location"] <- "pos"
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
merged$P.value <- 10 ^ (-merged$neg_log_10_p_value)
sig_a_alle=merged[,-c('ref_vcf','alt_vcf','snp_id_originalVCF')] #提取出sig_a
colnames(sig_a_alle)[colnames(sig_a_alle) == "standard_error"] <- "se"
###求zscore
sig_a_alle$zscore <- sig_a_alle$beta / sig_a_alle$se
sig_a_alle_clean <- sig_a_alle[!is.na(sig_a_alle$zscore) & sig_a_alle$beta != 0 & sig_a_alle$se != 0, ]

fwrite(sig_a_alle_clean, file = "/node1/liuxy/MediCell/trait_compare_2/eGFR/usefile/eGFR_Karczewski_aligned_with_pbmcdb.txt.gz",sep='\t')
  
