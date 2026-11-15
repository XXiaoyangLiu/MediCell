gene_annot <- read.table('/node1/zhaoq/scTWAS_PBMC/predictDB/model_training/merged_vcf_rsid/allgene/gene_annot/gene_annot.parsed.txt',header = T)
gene_annot$gene_type <- gsub("\";$", "", gene_annot$gene_type)
str(gene_annot$gene_type)
colnames(gene_annot)[6]='gene_type'

rows_plus <- gene_annot$strand == "+"
gene_annot[rows_plus, c("start", "end")] <- gene_annot[rows_plus, "start"]
gene_annot[!rows_plus, c("start", "end")] <- gene_annot[!rows_plus, "end"]
gene_annot$gene_name <- NULL
colnames(gene_annot)[2] <- "gene_name"

gene_annot <- gene_annot[,-7]

write.table(gene_annot,'/node1/liuxy/MediCell/pbmc/predictdb/usefile/gene/gene_annot.parsed.txt',sep='\t',row.names=F)
###attention!!!!后期在s-predixcan中发现一个错误，溯源与gene_annot.parsed.txt中在chr11和chr13中有两个LSP1，核实后发现ENSG00000269099应该为LSP1P1，已手动在gene_annot.parsed.txt中更正
```

#### 1.3.3 remove 10 samples for CIBERSORTx reference
```r
library(data.table)
library(stringr)
vcf <- fread("/node1/liuxy/chplot/combind_plasma_B/genotype/usefile/step5.vcf",header = TRUE)
dim(vcf) ##1035 samples, only use 971 samples
sample971 <- read.table('/node1/zhaoq/scTWAS_PBMC/CIBERSORTx/pbmc_14_CellType_971samples_pseudobulk_expCount.txt')
sample971 <- gsub("^X", "", colnames(sample971))  
col_list <- c("varID",sample971)
vcf <- vcf[, ..col_list]
dim(vcf)
var_id = vcf$varID
#vcf = vcf[,-1]
boolvec <- rowSums(vcf[,-1]) == 0
#vcf <- vcf[!boolvec,]
dim(vcf) 
var_id <- vcf$varID
chr <- str_split_i(var_id, "_", 1)
pos <- str_split_i(var_id, "_", 2)
ref <- str_split_i(var_id, "_", 3)
alt <- str_split_i(var_id, "_", 4)
ori <- paste("snp", chr, pos, sep  = "_")
annot <- data.table(chr, pos,var_id, ref, alt, ori, var_id)
colnames(annot) <- c("chromosome","Pos","VariantID","ref_vcf","alt_vcf","snp_id_originalVCF","RSID")
fwrite(vcf, "/node1/liuxy/chplot/combind_plasma_B/genotype/usefile/step6.vcf", sep = "\t")
fwrite(annot, "/node1/liuxy/chplot/combind_plasma_B/genotype/usefile/snp_annotation.txt", sep = "\t")
