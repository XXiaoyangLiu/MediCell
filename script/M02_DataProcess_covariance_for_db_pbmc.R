#pbmc
#nodes=tc6000
#source ~/.bashrc
#conda activate r3.6.3
library(stringr)
library(parallel)
library(peer)
library(dplyr)
setwd("/node1/liuxy/chplot/combind_plasma_B/enigma/output/peer/")
covs <- read.table("/node1/liuxy/sc_eQTL_science2022/OneK1K_sex_age_10PCs.txt",header = T) #only need 971 samples
sample971 <- read.table("/node1/liuxy/chplot/combind_plasma_B/enigma/usefile/pbmc_6_CellType_971samples_pseudobulk_TPM.txt",header = T)
cbs_frax <- read.table("/node1/zhaoq/scTWAS_PBMC/ENIGMA/new_result_54genes/CIBERSORTx_Job129_Adjusted_pbmc_54genes.txt",header = T)
cbs_frax <- cbs_frax %>% select(-(last_col(2):last_col()))
cbs_frax$B_cell_merged <- cbs_frax[["Plasma"]] + cbs_frax[["B"]]
cbs_frax <- cbs_frax[, !colnames(cbs_frax) %in% c("Plasma", "B")]
colnames(cbs_frax)[colnames(cbs_frax) == 'B_cell_merged'] <- 'B'
rownames(cbs_frax) <- cbs_frax$Mixture
cbs_frax <- cbs_frax[,-1]
cell_fractions <- cbs_frax[, c("Mono", "NK", "CD8T", "CD4T", "B")]
cell_fractions_t <- t(cell_fractions)
rownames(cell_fractions_t) <- c("Mono", "NK", "CD8T", "CD4T", "B")

covs <- covs[,c("id",colnames(sample971[,2:ncol(sample971)]))] 
dim(covs)  # age sex and 6 genotype princepal components
rownames(covs) <- covs$id
##列名为covirates，行名sample
covs <- as.matrix(t(covs[,-1]))
covs_colnames <- colnames(covs)
common_samples <- intersect(colnames(covs), colnames(cell_fractions_t))
cell_fractions_t_aligned <- cell_fractions_t[, common_samples, drop = FALSE]
covs_aligned <- covs[, common_samples, drop = FALSE]
covs_new <- rbind(covs_aligned, cell_fractions_t_aligned)

####for each cell type
celltypes <- c("B","CD4T","CD8T","NK","Mono","DC")
peer_fctr_clcl <- function(i){
  ctname <- prmtr[i,1]
  thiscovs <- covs 
  expr.INT = read.csv(paste0("/node1/liuxy/chplot/combind_plasma_B/enigma/output/count/tpm0.1sample0.2/", ctname, "_gene_exp.csv"), header=T)
  dim(expr.INT) 
  expr.INT <- as.matrix(expr.INT)
  # The data matrix is assumed to have N rows and G columns, where N is the number of samples, and G is the number of genes.Now we can create the model object,
  model = PEER()
  # set the observed data,
  # PEER_setPhenoMean(model,as.matrix(expr)) (NULL response means no error here)
  PEER_setPhenoMean(model,expr.INT)
  dim(PEER_getPhenoMean(model))
  PEER_setNk(model,2) 
  PEER_getNk(model)
  PEER_setCovariates(model, thiscovs)
  PEER_setAdd_mean(model,TRUE)
  PEER_update(model)
  factors = PEER_getX(model)
  dim(factors)
  weights = PEER_getW(model)
  dim(weights)
  precision = PEER_getAlpha(model)
  dim(precision)
  residuals = PEER_getResiduals(model)
  dim(residuals)
  res <- list(factors, weights, precision, residuals)
  names(res) <- c("factors", "weights", "precision", "residuals")
  resfilename = paste0(ctname, "_", "PCs_and_age_gender_peer_res.RData")
  save(res, file = resfilename)
  print(resfilename)
  # return(res)
}

prmtr <- expand.grid(celltypes) 
mclapply(seq(nrow(prmtr)), peer_fctr_clcl, mc.cores = 6) #后台提交程序时，服务器申请7个核
#整理成cov



#glom
# source ~/.bashrc
# conda activate r3.6.3
library(stringr)
library(parallel)
library(peer)
library(dplyr)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
CELLTYPE <- args[1]

setwd("/node1/liuxy/test_peer_glom/exp0.2/peer/output/")
load("/node1/zhaoq/scTWAS_Kidney/predictDB/kidney_KS/glom/QTLData.RData")
covs <- AAglomPhe2[,c("Batch","Pairend","Age","Gender","RIN","PC1","PC2","PC3","PC4","PC5","Sclerosis")]
del <- c("Pairend","Gender","Batch")
covleft <- covs[,!(colnames(covs)%in% del)]
add2 <- model.matrix(~Pairend+Gender+Batch,data=covs)[,-1]
yes <- summary(AAglomPhe2$Batch) < 5
bat <- c(1:25)
delbatch <- paste0("Batch",bat[yes])
add <- add2[,!(colnames(add2) %in% c(delbatch))]
allcov <- cbind(covleft,add)
finalcov <- as.matrix(allcov)
sample <- fread('/node1/liuxy/chplot/kidney_new/output/glom/exp/Bcell_transformed_expression.txt')
all(colnames(finalcov)==sample$v1)
dim(finalcov)  #303 21

#covs_colnames <- rownames(finalcov)
#common_samples <- intersect(rownames(finalcov), colnames(cell_fractions_t))
#cell_fractions_t_aligned <- cell_fractions_t[, common_samples, drop = FALSE]
#covs_aligned <- finalcov[common_samples,, drop = FALSE]
#covs_new <- cbind(covs_aligned, t(cell_fractions_t_aligned))
#finalcov <- t(covs_new)
ctname <- CELLTYPE
thiscovs <- finalcov 

expr.INT = read.csv(paste0("/node1/liuxy/chplot/kidney_new/output/glom/exp/", ctname, "_gene_exp.csv"), header=T)
expr.INT <- as.matrix(expr.INT)
model = PEER()
PEER_setPhenoMean(model,expr.INT)
PEER_setNk(model,60) 
PEER_setCovariates(model, thiscovs)
PEER_setAdd_mean(model,TRUE)
PEER_update(model)
factors = PEER_getX(model)
weights = PEER_getW(model)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)

res <- list(factors, weights, precision, residuals)
names(res) <- c("factors", "weights", "precision", "residuals")
resfilename = paste0(ctname,"_peer_res.RData")
save(res, file = resfilename)

#整理成cov
CELLTYPES <- c("Bcell","Endo","MC","MNP","NK","NKT","Podo","Th_cell","CD8T")
for(ct in CELLTYPES) {
load(paste0("/node1/liuxy/test_peer_glom/exp0.2/peer/output/",ct,"_peer_res.RData"))
peer_fctr <- res$factors
peer_fctr <- peer_fctr[,-18]
peer_fctr_rmcol <- t(peer_fctr)
tmp <- read.table("/node1/liuxy/chplot/kidney_new/output/glom/exp/Bcell_transformed_expression.txt")
rownames(tmp) <- gsub("^X", "", rownames(tmp))  
samples_ordr <- rownames(tmp)
rm(tmp,res)
gc()
colnames(peer_fctr_rmcol) <- samples_ordr
rownames(peer_fctr_rmcol) <- c("Age","RIN","PC1","PC2","PC3","PC4","PC5","Sclerosis","Pairend1","Gender1","Batch2","Batch3","Batch4","Batch5","Batch6","Batch8","Batch11",paste0("peer", 1:60))
peer_15 <- peer_fctr_rmcol[1:32,]
peer_15 <- as.data.frame(peer_15)
path <- "/node1/liuxy/MediCell/kidney_3/glom/predictdb/usefile/cov/"
write.table(peer_15, paste0(path, ct,"_peer15.txt"), sep="\t", row.names=T) 
}

#tub
# source ~/.bashrc
# conda activate r3.6.3
library(stringr)
library(parallel)
library(peer)
library(dplyr)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
CELLTYPE <- args[1]

setwd("/node1/liuxy/MediCell/kidney/tub/predictdb_2/usefile/cov/")
load("/node1/zhaoq/scTWAS_Kidney/predictDB/kidney_KS/tub/QTLData.RData")
covs <- AAtubPhe2[,c("Pairend","Batch","Age","Gender","Site","RIN","PC1","PC2","PC3","PC4","PC5","Fibrosis")]
del <- c("Pairend","Gender","Site","Batch")
covleft <- covs[,!(colnames(covs)%in% del)]
add2 <- model.matrix(~Pairend+Gender+Batch,data=covs)[,-1]
yes <- summary(AAtubPhe2$Batch) < 5
bat <- c(1:25)
delbatch <- paste0("Batch",bat[yes])
add <- add2[,!(colnames(add2) %in% c(delbatch))]
allcov <- cbind(covleft,add)
finalcov <- as.matrix(allcov)
sample <- fread('/node1/liuxy/chplot/kidney_new/output/tub/exp/Bcell_transformed_expression.txt')
all(colnames(finalcov)==sample$v1)
dim(finalcov)

#covs_colnames <- rownames(finalcov)
#common_samples <- intersect(rownames(finalcov), colnames(cell_fractions_t))
#cell_fractions_t_aligned <- cell_fractions_t[, common_samples, drop = FALSE]
#covs_aligned <- finalcov[common_samples,, drop = FALSE]
#covs_new <- cbind(covs_aligned, t(cell_fractions_t_aligned))
#finalcov <- t(covs_new)
ctname <- CELLTYPE
thiscovs <- finalcov 

expr.INT = read.csv(paste0("/node1/liuxy/chplot/kidney_new/output/tub/exp/", ctname, "_gene_exp.csv"), header=T)
expr.INT <- as.matrix(expr.INT)
model = PEER()
PEER_setPhenoMean(model,expr.INT)
PEER_setNk(model,60) 
PEER_setCovariates(model, thiscovs)
PEER_setAdd_mean(model,TRUE)
PEER_update(model)
factors = PEER_getX(model)
weights = PEER_getW(model)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)

res <- list(factors, weights, precision, residuals)
names(res) <- c("factors", "weights", "precision", "residuals")
resfilename = paste0(ctname,"_peer_res.RData")
save(res, file = resfilename)
print(paste("finish:", resfilename))

tmp <- read.table("/node1/liuxy/chplot/kidney_new/output/tub/exp/CD8T_transformed_expression.txt")
rownames(tmp) <- gsub("^X", "", rownames(tmp))  
samples_ordr <- rownames(tmp)
rm(tmp)
gc()

#整理成cov
library(data.table)
CELLTYPES <- c("Bcell","CD8T","CDICA","CDICB","DCT","LOH","MNP","NK","NKT","PT","Th_cell","Plasmacytoid_DC")
for(ct in CELLTYPES) {
load(paste0("/node1/liuxy/MediCell/kidney/tub/predictdb_2/usefile/cov/",ct,"_peer_res.RData"))
peer_fctr <- res$factors
peer_fctr <- peer_fctr[,-22]
peer_fctr_rmcol <- t(peer_fctr)
tmp <- read.table("/node1/liuxy/chplot/kidney_new/output/tub/exp/Bcell_transformed_expression.txt")
rownames(tmp) <- gsub("^X", "", rownames(tmp))  
samples_ordr <- rownames(tmp)
rm(tmp,res)
gc()
colnames(peer_fctr_rmcol) <- samples_ordr
rownames(peer_fctr_rmcol) <- c("Age","RIN","PC1","PC2","PC3","PC4","PC5","Fibrosis","Pairend1","Gender1","Batch2","Batch3","Batch4","Batch5","Batch6","Batch7","Batch8","Batch9","Batch10","Batch11","Batch12",paste0("peer", 1:60))
peer_15 <- peer_fctr_rmcol[1:32,]
peer_15 <- as.data.frame(peer_15)
path <- "/node1/liuxy/MediCell/kidney_3/tub/predictdb/usefile/cov/"
write.table(peer_15, paste0(path, ct,"_peer15.txt"), sep="\t", row.names=T) 
}
