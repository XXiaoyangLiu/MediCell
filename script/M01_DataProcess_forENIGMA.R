 #pseudobulk:Sum
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(Matrix)
#singlecellexpriment装不上，所以用基础包直接做
data <- readRDS('/node1/liuxy/chplot/combind_plasma_B/cibersort/usefile/pbmc_6_CellType_971samples.rds')
pesudo_counts <- data@assays$RNA@counts
pesudo_metadata <- data@meta.data
pesudo_counts <- data@assays$RNA@counts
pesudo_metadata <- data@meta.data
pesudo_metadata$sample_id <- factor(pesudo_metadata$individual)
sample_ids <- levels(pesudo_metadata$sample_id)
pb <- matrix(0, 
             nrow = nrow(pesudo_counts), 
             ncol = length(sample_ids),
             dimnames = list(rownames(pesudo_counts), sample_ids))

for (i in seq_along(sample_ids)) {
    sid <- sample_ids[i]
    cells_in_sample <- which(pesudo_metadata$sample_id == sid)
    if (length(cells_in_sample) > 0) {
        pb[, i] <- rowSums(pesudo_counts[, cells_in_sample, drop = FALSE])
    }
}
write.table(as.matrix(pb), file = '/node1/liuxy/chplot/combind_plasma_B/cibersort/usefile/pbmc_6_CellType_971samples_pseudobulk_expCount.txt', sep = "\t", row.names = TRUE, col.names = NA)
#tpm:
mycounts<-read.csv("/node1/liuxy/chplot/combind_plasma_B/cibersort/usefile/pbmc_6_CellType_971samples_pseudobulk_expCount.txt",sep='\t',check.names = FALSE,row.name = "")
gene_effective_len = read.csv('/node1/zhaoq/scTWAS_PBMC/CIBERSORTx/geneid_efflen_gtf.txt',sep='\t')
gene_effective_len=na.omit(gene_effective_len)
mycounts$geneLength = gene_effective_len[match(rownames(mycounts), gene_effective_len$geneid),2]
kb <- mycounts$geneLength / 1000
head(kb)
countdata <- mycounts[,1:971]
head(countdata)
rpk <- countdata / kb
head(rpk)
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
#20245 genes, 971 samples
tpm=as.data.frame(tpm)
tpm$GeneID <- rownames(tpm)  
tpm <- tpm[, c("GeneID", setdiff(names(tpm), "GeneID"))]
write.table(tpm, file = '/node1/liuxy/chplot/combind_plasma_B/enigma/usefile/pbmc_6_CellType_971samples_pseudobulk_TPM.txt', sep = "\t",quote = FALSE, row.names = FALSE)

#ref_10sample
extracted_data_10samples= readRDS('/node1/liuxy/chplot/combind_plasma_B/cibersort/usefile/pbmc_6_CellType_10samples.rds')
tmp_scRNA_10sample = GetAssayData(extracted_data_10samples, layer = "counts", assay = "RNA")
colnames(tmp_scRNA_10sample) = extracted_data_10samples$cell_label
tmp_scRNA_10sample = data.frame(tmp_scRNA_10sample, check.names = F)   # 一定要有check.names = F，不然的话后面重新命名列名会有问题，matrix对象不允许列名相同，但data.frame对象允许
tmp_scRNA_10sample = data.frame(gene = rownames(tmp_scRNA_10sample), tmp_scRNA_10sample, check.names = F)
table(colnames(tmp_scRNA_10sample))
write.table(tmp_scRNA_10sample, file = '/node1/liuxy/chplot/combind_plasma_B/cibersort/usefile/pbmc_6_CellType_10samples_expCount.txt', sep = "\t") 

mycounts  <- read.csv('/node1/liuxy/chplot/combind_plasma_B/cibersort/usefile/pbmc_6_CellType_10samples_expCount.txt',sep = '\t') #这个文件来自前一个代码块生成
mycounts$gene <- NULL
gene_effective_len <- read.csv('/node1/zhaoq/scTWAS_PBMC/CIBERSORTx/geneid_efflen_gtf.txt', sep = '\t')
gene_effective_len <- na.omit(gene_effective_len)

mycounts$geneLength <- gene_effective_len[match(rownames(mycounts), gene_effective_len$geneid), 2]
kb <- mycounts$geneLength / 1000
countdata <- mycounts[, 1:(ncol(mycounts) - 1)]
rpk <- countdata / kb
tpm <- t(t(rpk) / colSums(rpk) * 1000000)
tpm <- as.data.frame(tpm)
tpm$GeneID <- rownames(tpm)
tpm <- tpm[, c("GeneID", setdiff(names(tpm), "GeneID"))]
write.table(tpm, file = '/node1/liuxy/chplot/combind_plasma_B/enigma/usefile/pbmc_6_CellType_10samples_pseudobulk_TPM.txt', sep = "\t", quote = FALSE, row.names = FALSE)

forENIGMA_ref <- as.data.frame(fread('/node1/liuxy/chplot/combind_plasma_B/enigma/usefile/pbmc_6_CellType_10samples_pseudobulk_TPM.txt', check.name=F))
rn <- forENIGMA_ref$GeneID
cn <- colnames(forENIGMA_ref)[-1]
cell_types_sorted <- sort(names(table(cn))) 
profile_pbmc <- NULL
cell_type <- gsub(".[0-9]+$", "", cn)
unique_ct <- unique(cell_type)

profile_list <- list()
for (ct in unique_ct) {
    cols <- which(cell_type == ct) + 1
    profile_list[[ct]] <- rowMeans(forENIGMA_ref[, cols, drop = FALSE])
}
profile_pbmc <- as.data.frame(profile_list)
rownames(profile_pbmc) <- rn
write.table(profile_pbmc, file="/node1/liuxy/chplot/combind_plasma_B/enigma/usefile/TPM_ENIGMA_ref_6CT.matrx")

#### 2. kidney
#scRNA-seq
gene_info$gene_id_clean <- sub("\\..*$", "", gene_info$gene_id)
young_gene <- read.table("/node1/liuxy/chplot/kidney_new/usefile/scTWAS_info/NG_kidney_fromSX/Young/genes.tsv", header = FALSE, sep = "\t")
gene_map <- setNames(young_gene$V1, young_gene$V2)
###tub
###reference（tpm），这个profile_tub_tpm其实就是前面的tub.sc，有时间的时候统一变量名
profile_tub_tpm=as.data.frame(fread('/public/home/Shenglab/DataShare/ToZhaoqiong/scTWAS_info/NG_kidney_fromSX/tub.scRNAseq.ref.exp.tpm.txt', check.name=F)) #21996 genes
profile_tub_tpm$GeneID <- gene_map[profile_tub_tpm$GeneSymbol]
###只保留蛋白编码基因,且必须在bulk内
profile_tub_tpm <- profile_tub_tpm[profile_tub_tpm$GeneID %in% bulk_tub_new_tpm$Gene,] #bulk_ref是2.CIBERSORTx.md的line329, bulk_tub_new_tpm是2.CIBERSORTx.md的line343
#去掉列名中的点和点后面的数字  
colnames(profile_tub_tpm) <- gsub("\\.\\d+$", "", colnames(profile_tub_tpm))  
rn <- profile_tub_tpm$GeneID
cn <- colnames(profile_tub_tpm)[-c(1,5751)]
profile_tub <- NULL
for(i in seq(length(table(cn)))){
  ct <- names(table(cn))[i]
  ordr <- which(ct==cn)+1 #因为第一列是基因名称
  item <- rowMeans(profile_tub_tpm[,ordr])
  profile_tub <- data.table(profile_tub, item)
}
colnames(profile_tub) <- names(table(cn))
profile_tub <- as.data.frame(profile_tub)
rownames(profile_tub) <- rn
profile_tub <- profile_tub[ ,!colnames(profile_tub)%in% c("Neutrophil", "Mast_cell")] #Neutrophil和Mast_cell在356个样本中cell fraction均小于10%，故在ENIGMA这步去除，后续不参与TWAS分析
write.table(profile_tub, file="/node1/liuxy/chplot/kidney_new/usefile/TPM_ENIGMA_ref_withoutNeu_Mast_ProtCoding.matrx")

##glom
profile_glom_tpm = as.data.frame(fread('/public/home/Shenglab/DataShare/ToZhaoqiong/scTWAS_info/NG_kidney_fromSX/glom.scRNAseq.ref.exp.tpm.txt', check.name=F))  
profile_glom_tpm$GeneID <- gene_map[profile_glom_tpm$GeneSymbol]
profile_glom_tpm_2 <- profile_glom_tpm[profile_glom_tpm$GeneID %in% bulk_glom_tpm_new$Gene,]  #bulk_ref是2.CIBERSORTx.md的line329,bulk_glom_tpm_new是2.CIBERSORTx.md的line349
# 去掉列名中的点和点后面的数字  
colnames(profile_glom_tpm) <- gsub("\\.\\d+$", "", colnames(profile_glom_tpm))  

rn <- profile_glom_tpm$GeneID  
cn <- colnames(profile_glom_tpm)[-c(1,6008)]  
profile_glom <- NULL  

for(i in seq(length(table(cn)))){  
  ct <- names(table(cn))[i]  
  ordr <- which(ct == cn) + 1 # 因为第一列是基因名称  
  item <- rowMeans(profile_glom_tpm[, ordr])  
  profile_glom <- data.table(profile_glom, item)  
}  

colnames(profile_glom) <- names(table(cn))  
profile_glom <- as.data.frame(profile_glom)  
rownames(profile_glom) <- rn  

profile_glom <- profile_glom[, !colnames(profile_glom) %in% c("Neutrophil", "Mast_cell","Plasmacytoid_DC")]  #Neutrophil、Mast_cell和Plasmacytoid_DC在303个样本中cell fraction均小于10%，故在ENIGMA这步去除，后续不参与TWAS分析
write.table(profile_glom, file="/node1/liuxy/chplot/kidney_new/usefile/TPM_ENIGMA_ref_withoutNeu_Mast_pDC_ProtCodinng.matrix")  

### Bulk处理：gene symbol unique;保留蛋白编码基因（for ENIGMA）
##参考基因组GeneID与GeneSymbol转换
library(stringr)  
library(dplyr)  
gtf_df <- read.delim("/public/home/Shenglab/DataShare/Reference/human/hg19/gencode.v19.annotation.gtf", comment.char = "#",
                     header = FALSE, sep = "\t", stringsAsFactors = FALSE)

colnames(gtf_df) <- c("chr", "source", "feature", "start", "end",
                      "score", "strand", "frame", "attribute")
gene_df <- gtf_df %>% filter(feature == "gene")

gene_df$gene_id <- str_extract(gene_df$attribute, "gene_id\\s+([^;]+);")
gene_df$gene_id <- str_replace_all(gene_df$gene_id, "gene_id\\s+|;", "")
gene_df$gene_name <- str_extract(gene_df$attribute, "gene_name\\s+([^;]+);")
gene_df$gene_name <- str_replace_all(gene_df$gene_name, "gene_name\\s+|;", "")
gene_df$gene_type <- str_extract(gene_df$attribute, "gene_type\\s+([^;]+);")
gene_df$gene_type <- str_replace_all(gene_df$gene_type, "gene_type\\s+|\"|;", "")
gene_df$gene_status <- str_extract(gene_df$attribute, "gene_status\\s+([^;]+);")
gene_df$gene_status <- str_replace_all(gene_df$gene_status, "gene_status\\s+|\"|;", "")
gene_df$transcript_type <- str_extract(gene_df$attribute, "transcript_type\\s+([^;]+);")
gene_df$transcript_type <- str_replace_all(gene_df$transcript_type, "transcript_type\\s+|\"|;", "")
gene_df$transcript_status <- str_extract(gene_df$attribute, "transcript_status\\s+([^;]+);")
gene_df$transcript_status <- str_replace_all(gene_df$transcript_status, "transcript_status\\s+|\"|;", "")
gene_df$level <- str_extract(gene_df$attribute, "level\\s+([^;]+);")
gene_df$level <- str_replace_all(gene_df$level, "level\\s+|\"|;", "")
gene_sorted <- gene_df %>%
  select(chr, start, end, strand, gene_id, gene_name, gene_type, gene_status, transcript_status) %>%
  arrange(chr, start)
gene_info <-gene_sorted[gene_sorted$gene_type == "protein_coding", ] 
gene_ids <- gene_info$gene_id   
write.csv(gene_sorted, "/node1/liuxy/index/gencode_v19.csv", row.names=FALSE)
bulk_ref <- data.frame(GeneID = gene_sorted$gene_id, GeneSymbol = gene_sorted$gene_name, stringsAsFactors = FALSE)  #20345 protein coding genes
#bulk_ref <- bulk_ref[!duplicated(bulk_ref),]
bulk_ref <- bulk_ref[bulk_ref$GeneID %in% gene_ids,]
bulk_ref_2 <- bulk_tub_new_tpm[,c(1,ncol(bulk_tub_new_tpm))]
####这里bulk_tub_tpm_2/bulk_glom_tpm_2来自GSE173343，仅用来对应GeneSymbol-GeneID，已经检查过了里面的geneid和symbol的对应关系是正确的,顺序也是一致的
bulk_tub_tpm <- read.table("/public/home/Shenglab/DataShare/ToZhaoqiong/scTWAS_info/NG_kidney_fromSX/tub.bulk.EXP.tpm.Name.txt", header = T)
bulk_tub_tpm_2 <- read.table("/node1/liuxy/chplot/kidney_new/usefile/TubSamplesTPM.txt",header=T) #42855genes  
bulk_tub_tpm$Gene <- bulk_tub_tpm_2$Gene
bulk_tub_new_tpm <- subset(bulk_tub_tpm, Gene %in% gene_info$gene_id) #18953genes
bulk_tub_new_tpm$Gene <- sub("\\..*$", "", bulk_tub_new_tpm$Gene)
bulk_glom_tpm_new <- subset(bulk_glom_tpm, Gene %in% gene_info$gene_id)
bulk_glom_tpm_new$Gene <- sub("\\..*$", "", bulk_glom_tpm_new$Gene)

bulk_glom_tpm <- read.table("/public/home/Shenglab/DataShare/ToZhaoqiong/scTWAS_info/NG_kidney_fromSX/glom.bulk.EXP.tpm.Name.txt", header = T)
bulk_glom_tpm_2 <- read.table("/node1/liuxy/chplot/kidney_new/usefile/GlomSamplesTPM.txt",header=T) #42827genes  
bulk_glom_tpm$Gene <- bulk_glom_tpm_2$Gene
bulk_glom_new_tpm <- subset(bulk_glom_tpm, Gene %in% gene_info$gene_id) #18905 genes 
bulk_glom_new_tpm$Gene <- sub("\\..*$", "", bulk_glom_new_tpm$Gene)
write.table(bulk_tub_new_tpm,file='/node1/liuxy/chplot/kidney_new/usefile/tub.bulk.EXP.tpm.Name.Unique.ProtCoding.txt',sep='\t', row.names = FALSE)
write.table(bulk_glom_new_tpm, file='/node1/liuxy/chplot/kidney_new/usefile/glom.bulk.EXP.tpm.Name.Unique.ProtCoding.txt', sep='\t', row.names = FALSE)
