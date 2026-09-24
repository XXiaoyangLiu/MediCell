#change trait from UKB(Whole-genome sequencing of 490,640 UK Biobank participants.)
#这篇文章用的trait包括Glu(GCST90474364) SBP(GCST90474636) UA(GCST90474433) uAlb(GCST90474446) uCr(GCST90474438) uK(GCST90474443) uNa(GCST90474451)
#全部都是grch38的,在之前已经比对过hg19
library(data.table)
library(dplyr)
#ref alt校正的方向
db_anno <- fread("/node1/liuxy/chplot/combind_plasma_B/genotype/usefile/snp_annotation.txt")
db_anno$Pos <- as.numeric(db_anno$Pos)
db_anno <- db_anno[db_anno$chrom %in% c(1:22)]
db_anno <- as.data.table(db_anno)
setnames(db_anno, c("ref_vcf", "alt_vcf"), c("ref_db", "alt_db"))
#所有来自UKB的GWAS
traits <- c("Glu","UA","uAlb","uCr","uK","uNa","SBP")
for (trait in traits) {
     trait_dir <- paste0("/node1/liuxy/MediCell/trait_compare/metabolic/usefile/",trait,"/",trait,"_hg19.tsv.gz")
     trait_UKB <- fread(trait_dir)
     trait_UKB <- as.data.table(trait_UKB)
     trait_UKB$zscore <- trait_UKB$beta / trait_UKB$standard_error
     trait_UKB$other_allele <- toupper(trait_UKB$other_allele)
     trait_UKB$effect_allele <- toupper(trait_UKB$effect_allele)
     trait_UKB$snp_id_originalVCF <- paste0("snp_",trait_UKB$chr_hg19,"_",trait_UKB$pos_hg19)
     setnames(trait_UKB, c("base_pair_location"), c("Pos_38"))
     merged_data <- merge(trait_UKB, db_anno, by = "snp_id_originalVCF", all.x = TRUE)
     merged_data$need_flip <- FALSE
     consistent <- (toupper(merged_data$effect_allele) == toupper(merged_data$alt_db) & 
               toupper(merged_data$other_allele) == toupper(merged_data$ref_db))
     need_flip <- (toupper(merged_data$effect_allele) == toupper(merged_data$ref_db) & 
              toupper(merged_data$other_allele) == toupper(merged_data$alt_db)) 
     mismatch <- !consistent & !need_flip
     merged_data$need_flip <- need_flip
     merged_data[need_flip == TRUE, 
           `:=`(effect_allele_original = effect_allele,
                other_allele_original = other_allele)]
     merged_data[need_flip == TRUE, 
           `:=`(effect_allele = other_allele_original,
                other_allele = effect_allele_original,
                beta = -beta,
                zscore = -zscore,
                effect_allele_frequency = 1 - effect_allele_frequency)]
     trait_UKB_aligned <- merged_data[!is.na(ref_db) & !mismatch, 
                             .(chr_hg19, pos_hg19, ref_db, alt_db, 
                               RSID,effect_allele_frequency,
                               standard_error, p_value, n, 
                               beta, zscore)] 
     colnames(trait_UKB_aligned) <- c("chr", "pos", "ref", "alt", 
                              "rsID", "Freq1", 
                              "se", "P.value", "n_total_sum", 
                              "beta", "zscore")
     trati_out <- paste0("/node1/liuxy/MediCell/trait_compare_2/other_trait/usefile/",trait,"_aligned_with_db_pbmc.txt.gz")
     fwrite(trait_UKB_aligned, trati_out, sep = "\t")
}


#tub
db_anno <- fread("/node1/liuxy/MediCell/kidney/glom/predictdb/usefile/glom_snp_annotation.txt")
db_anno$Pos <- as.numeric(db_anno$Pos)
db_anno <- db_anno[db_anno$chrom %in% c(1:22)]
db_anno <- as.data.table(db_anno) 
setnames(db_anno, c("ref_vcf", "alt_vcf"), c("ref_db", "alt_db"))
db_anno[, snp_id_originalVCF := paste0("snp_", sub("_[A-Z]+_[A-Z]+$", "", snp_id_originalVCF))]
#所有来自UKB的GWAS
traits <- c("Glu","UA","uAlb","uCr","uK","uNa","SBP")
for (trait in traits) {
     trait_dir <- paste0("/node1/liuxy/MediCell/trait_compare/metabolic/usefile/",trait,"/",trait,"_hg19.tsv.gz")
     trait_UKB <- fread(trait_dir)
     trait_UKB <- as.data.table(trait_UKB)
     trait_UKB$zscore <- trait_UKB$beta / trait_UKB$standard_error
     trait_UKB$other_allele <- toupper(trait_UKB$other_allele)
     trait_UKB$effect_allele <- toupper(trait_UKB$effect_allele)
     trait_UKB$snp_id_originalVCF <- paste0("snp_",trait_UKB$chr_hg19,"_",trait_UKB$pos_hg19)
     setnames(trait_UKB, c("base_pair_location"), c("Pos_38"))
     merged_data <- merge(trait_UKB, db_anno, by = "snp_id_originalVCF", all.x = TRUE)
     merged_data$need_flip <- FALSE
     consistent <- (toupper(merged_data$effect_allele) == toupper(merged_data$alt_db) & 
               toupper(merged_data$other_allele) == toupper(merged_data$ref_db))
     need_flip <- (toupper(merged_data$effect_allele) == toupper(merged_data$ref_db) & 
              toupper(merged_data$other_allele) == toupper(merged_data$alt_db)) 
     mismatch <- !consistent & !need_flip
     merged_data$need_flip <- need_flip
     merged_data[need_flip == TRUE, 
           `:=`(effect_allele_original = effect_allele,
                other_allele_original = other_allele)]
     merged_data[need_flip == TRUE, 
           `:=`(effect_allele = other_allele_original,
                other_allele = effect_allele_original,
                beta = -beta,
                zscore = -zscore,
                effect_allele_frequency = 1 - effect_allele_frequency)]
     trait_UKB_aligned <- merged_data[!is.na(ref_db) & !mismatch, 
                             .(chr_hg19, pos_hg19, ref_db, alt_db, 
                               rsid,effect_allele_frequency,
                               standard_error, p_value, n, 
                               beta, zscore)] 
     colnames(trait_UKB_aligned) <- c("chr", "pos", "ref", "alt", 
                              "RSID", "Freq1", 
                              "se", "P.value", "n_total_sum", 
                              "beta", "zscore")
     trati_out <- paste0("/node1/liuxy/MediCell/trait_compare_2/other_trait/usefile/",trait,"_aligned_with_db_kidney.txt.gz")
     fwrite(trait_UKB_aligned, trati_out, sep = "\t")
     print(trati_out)
}
