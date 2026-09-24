##grch37
library(data.table)
library(stringr)
snp_index <- fread("/node1/liuxy/MediCell/method_compare_2/usefile/GWAS_scpredixcan/chain/snp_index_hg38_sorted.txt", header = TRUE, col.names = c("chr", "pos_hg38", "ref_index", "alt_index", "snpID"))
snp_index <- unique(snp_index)
snp_index[, chr := as.numeric(chr)]
#GCST011096
a=fread("/node1/liuxy/MediCell/method_compare_2/usefile/GWAS_scpredixcan/GCST011096_hg19.tsv.gz")
a <- a[a$chr_hg38 %in% c(1:22)]
a$chr_hg38=as.character(a$chr_hg38)
merged_data <- merge(a, snp_index, 
                by.x = c("chr_hg38","pos_hg38"),
                by.y = c("chr", "pos_hg38"),   
                all = FALSE)
merged_data$need_flip <- FALSE
     consistent <- (toupper(merged_data$effect_allele) == toupper(merged_data$alt_index) & 
               toupper(merged_data$other_allele) == toupper(merged_data$ref_index))
     need_flip <- (toupper(merged_data$effect_allele) == toupper(merged_data$ref_index) & 
              toupper(merged_data$other_allele) == toupper(merged_data$alt_index)) 
     mismatch <- !consistent & !need_flip
     merged_data$need_flip <- need_flip
     merged_data <- merged_data[mismatch == FALSE]
     merged_data[need_flip == TRUE, 
           `:=`(effect_allele_original = effect_allele,
                other_allele_original = other_allele)]
     merged_data[need_flip == TRUE, 
           `:=`(effect_allele = other_allele_original,
                other_allele = effect_allele_original,
                beta = -beta)]
merged_data$effect_allele <-  toupper(merged_data$effect_allele)
merged_data$other_allele <-  toupper(merged_data$other_allele)
sig_a_alle=merged_data[,-c('odds_ratio','ci_lower','ci_upper','effect_allele_frequency','ref_index','alt_index','effect_allele_original','other_allele_original','need_flip')]
sig_a_alle$zscore=sig_a_alle$beta/sig_a_alle$standard_error
sig_a_alle_clean <- sig_a_alle[!is.na(sig_a_alle$zscore) & sig_a_alle$beta != 0 & sig_a_alle$standard_error != 0, ]
fwrite(sig_a_alle_clean, file = "/node1/liuxy/MediCell/method_compare_2/usefile/GWAS_scpredixcan/GCST011096_aligned.tsv",sep='\t')

#GCST90018917
a=fread("/node1/liuxy/MediCell/method_compare_2/usefile/GWAS_scpredixcan/GCST90018917_hg38.tsv.gz")
#colnames(a) <- c("chromosome","base_pair_location","effect_allele","other_allele","effect_allele_frequency","beta","standard_error","p_value","variant_id")
a <- a[a$chr_hg38 %in% c(1:22)]
a$chr_hg38=as.numeric(a$chr_hg38)
merged_data <- merge(a, snp_index, 
                by.x = c("chr_hg38","pos_hg38"),
                by.y = c("chr", "pos_hg38"),   
                all = FALSE)
merged_data$need_flip <- FALSE
     consistent <- (toupper(merged_data$effect_allele) == toupper(merged_data$alt_index) & 
               toupper(merged_data$other_allele) == toupper(merged_data$ref_index))
     need_flip <- (toupper(merged_data$effect_allele) == toupper(merged_data$ref_index) & 
              toupper(merged_data$other_allele) == toupper(merged_data$alt_index)) 
mismatch <- !consistent & !need_flip
merged_data$need_flip <- need_flip
merged_data <- merged_data[mismatch == FALSE]
merged_data[need_flip == TRUE, 
           `:=`(effect_allele_original = effect_allele,
                other_allele_original = other_allele)]
merged_data[need_flip == TRUE, 
           `:=`(effect_allele = other_allele_original,
                other_allele = effect_allele_original,
                beta = -beta,
                effect_allele_frequency = 1 - effect_allele_frequency)]
sig_a_alle=merged_data[,-c('variant_id','ref_index','alt_index','effect_allele_original','other_allele_original','need_flip')]
sig_a_alle$zscore=sig_a_alle$beta/sig_a_alle$standard_error
sig_a_alle_clean <- sig_a_alle[!is.na(sig_a_alle$zscore) & sig_a_alle$beta != 0 & sig_a_alle$standard_error != 0, ]
fwrite(sig_a_alle_clean, file = "/node1/liuxy/MediCell/method_compare_2/usefile/GWAS_scpredixcan/GCST90018917_aligned.tsv",sep='\t')

#GCST90476183
a=fread("/node1/liuxy/MediCell/pbmc_2/TWAS/usefile/GWAS/GCST90476183.tsv.gz")
a <- a[a$chromosome %in% c(1:22)]
a$chromosome=as.numeric(a$chromosome)
merged_data <- merge(a, snp_index, 
                by.x = c("chromosome", "base_pair_location"),
                by.y = c("chr", "pos_hg38"),   
                all = FALSE)
merged_data$beta = log(merged_data$odds_ratio)
merged_data$zscore_abs = qnorm( 1 - merged_data$p_value / 2 ) 
merged_data$standard_error = abs(merged_data$beta) / abs(merged_data$zscore_abs)
merged_data$need_flip <- FALSE
consistent <- (toupper(merged_data$effect_allele) == toupper(merged_data$alt_index) & 
                toupper(merged_data$other_allele) == toupper(merged_data$ref_index))
need_flip <- (toupper(merged_data$effect_allele) == toupper(merged_data$ref_index) & 
              toupper(merged_data$other_allele) == toupper(merged_data$alt_index)) 
mismatch <- !consistent & !need_flip
merged_data$need_flip <- need_flip
merged_data$mismatch <- mismatch
table(merged_data$need_flip,useNA = "always")
merged_data <- merged_data[mismatch == FALSE]
merged_data[, `:=`(
  effect_allele_original = effect_allele,
  other_allele_original = other_allele
)]
merged_data[need_flip == TRUE, 
           `:=`(effect_allele = other_allele_original,
                other_allele = effect_allele_original,
                beta = -beta,
                effect_allele_frequency = 1 - effect_allele_frequency)]
sig_a_alle=merged_data[,-c('r2','i2','direction','snp_id_originalVCF','effect_allele_original','other_allele_original','need_flip')]
sig_a_alle$zscore=sig_a_alle$beta/sig_a_alle$standard_error
sig_a_alle_clean <- sig_a_alle[!is.na(sig_a_alle$zscore) & sig_a_alle$beta != 0 & sig_a_alle$standard_error != 0, ]
fwrite(sig_a_alle_clean, file = "/node1/liuxy/MediCell/method_compare_2/usefile/GWAS_scpredixcan/GCST90476183_aligned.tsv",sep='\t')
 
#metal QC after metal
a <- fread("/node1/liuxy/MediCell/method_compare_2/usefile/GWAS_scpredixcan/SLE_metal_096_917_1831.tbl")
a$MarkerName
a <- a[! (grepl("\\+", Direction) & grepl("\\-", Direction))]
a <- a[!(a$HetISq > 50 & a$HetPVal < 0.05), ]
fwrite(a, file = "/node1/liuxy/MediCell/method_compare_2/usefile/GWAS_scpredixcan/SLE_metal_096_917_183_filt.tsv",sep='\t')
