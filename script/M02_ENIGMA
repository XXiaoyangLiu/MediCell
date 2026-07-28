library(ENIGMA)
library(dplyr)
setwd('/node1/liuxy/chplot/combind_plasma_B/enigma/output/')
profile_pbmc <- read.table('/node1/liuxy/chplot/combind_plasma_B/enigma/usefile/TPM_ENIGMA_ref_6CT.matrx')
profile_pbmc_7 <- read.table('/node1/zhaoq/scTWAS_PBMC/ENIGMA/new_result_54genes/Aggre/TPM_ENIGMA_ref_7CT.matrx')
Bulk <- read.csv('/node1/liuxy/chplot/combind_plasma_B/enigma/usefile/pbmc_6_CellType_971samples_pseudobulk_TPM.txt', sep='\t', row.names=1)
cbs_frax <- read.table("/node1/zhaoq/scTWAS_PBMC/ENIGMA/new_result_54genes/CIBERSORTx_Job129_Adjusted_pbmc_54genes.txt",header=T)

myenigma <- create_ENIGMA(bulk = as.matrix(Bulk), ref = as.matrix(profile_pbmc), ref_type = "aggre")
cbs_frax_for_enigma <- cbs_frax %>% select(-(last_col(2):last_col()))

ordr <- match(colnames(myenigma@bulk), cbs_frax_for_enigma$Mixture)
cbs_frax_for_enigma <- cbs_frax_for_enigma[ordr,]
all(colnames(myenigma@bulk)==cbs_frax_for_enigma$Mixture)

rownames(cbs_frax_for_enigma) <- cbs_frax_for_enigma$Mixture
cbs_frax_for_enigma <- cbs_frax_for_enigma[,-1]
ordr <- match(colnames(myenigma@ref), colnames(cbs_frax_for_enigma))
cbs_frax_for_enigma <- cbs_frax_for_enigma %>% select(ordr)

cbs_frax_for_enigma <- as.matrix(cbs_frax_for_enigma)
myenigma@result_cell_proportion <- cbs_frax_for_enigma
myenigma.bak <- myenigma
myenigma <- ENIGMA_trace_norm(myenigma.bak, alpha = 0.1, do_cpm=F, preprocess = "log")
save(myenigma, file = paste0('myenigma_alpha_0.1.RData'))
enigma.cse.test <- sce2array(myenigma, norm_output = F)
save(enigma.cse.test, file = paste0("./count/enigma.cse.trace_alpha_0.1.RData"))

all(enigma.cse.test>=0) # TRUE
all(enigma.cse.norm>=0) # FALSE

