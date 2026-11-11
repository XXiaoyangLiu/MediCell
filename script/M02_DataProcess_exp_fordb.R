### QC after ENGMA for calculate PEER factor
## R = 3.6.3
library(stringr)
library(parallel)
library(peer)
library(dplyr)

load("/node1/liuxy/chplot/combind_plasma_B/enigma/output/count/enigma.cse.trace_norm_alpha_0.1.RData")
cell_types <- dimnames(enigma.cse.test)[[3]]
filtered_genes <- list()
for (cell_type in cell_types) {
  expr_matrix <- enigma.cse.test[, , cell_type]
  n_samples_above_threshold <- rowSums(expr_matrix > 0.1)
  min_samples_1 <- ceiling(0.01 * ncol(expr_matrix))
  min_samples <- max(min_samples_1, 5)
  keep_genes <- n_samples_above_threshold >= min_samples
  
  filtered_genes[[cell_type]] <- rownames(expr_matrix)[keep_genes]
  message(sprintf("Cell type: %s, Genes kept: %d/%d (%.1f%%)", 
                  cell_type, 
                  sum(keep_genes), 
                  nrow(expr_matrix),
                  100 * sum(keep_genes)/nrow(expr_matrix)))
}

sapply(filtered_genes, length)
which(filtered_genes[['NK']]=='ENSG00000138378') 

filtered_matrices <- list()
cell_types <- dimnames(enigma.cse.test)[[3]]

for (cell_type in cell_types) {
  all_genes <- dimnames(enigma.cse.test)[[1]]
  keep_genes <- filtered_genes[[cell_type]]
  keep_indices <- which(all_genes %in% keep_genes)
  filtered_mat <- enigma.cse.test[keep_indices, , cell_type]
  filtered_matrices[[cell_type]] <- filtered_mat
}

cell_types <- names(filtered_matrices)
for (ctname in cell_types) {
  ct_matrix <- filtered_matrices[[ctname]]
  tmp <- as.data.frame(t(ct_matrix))
  tmp.INT <- matrix(0, nrow(tmp), ncol(tmp))
  for (i in 1:ncol(tmp)) {
    tmp.INT[, i] <- qqnorm(tmp[, i], plot.it = FALSE)$x
  }
  rownames(tmp.INT) <- rownames(tmp)
  colnames(tmp.INT) <- colnames(tmp)
  tmp.INT <- as.data.frame(tmp.INT)
  
  write.table(
    tmp.INT,
    file = paste0("/node1/liuxy/chplot/combind_plasma_B/enigma/output/count/tpm0.1sample0.2/", ctname, "_gene_exp.csv"),
    sep = ",",
    col.names = TRUE,
    row.names = FALSE
  )
  write.table(
    tmp.INT,
    file = paste0("/node1/liuxy/chplot/combind_plasma_B/enigma/output/count/tpm0.1sample0.2/", ctname, "_transformed_expression.txt"),
    sep = "\t",
    row.names = TRUE
  )
}

#glom
load("/node1/liuxy/chplot/kidney_new/output/glom/count/enigma.cse.trace_alpha_0.1.RData")
cell_types <- dimnames(enigma.cse.test)[[3]]
filtered_genes <- list()
for (cell_type in cell_types) {
  expr_matrix <- enigma.cse.test[, , cell_type]
  #每个基因在多少样本中表达 > 0.1
  n_samples_above_threshold <- rowSums(expr_matrix > 0.1)
  #至少需要的样本数 (1%),并且至少要有5个样本表达
  min_samples_1 <- ceiling(0.01 * ncol(expr_matrix))
  min_samples <- max(min_samples_1, 5)
  keep_genes <- n_samples_above_threshold >= min_samples
  
  filtered_genes[[cell_type]] <- rownames(expr_matrix)[keep_genes]
  message(sprintf("Cell type: %s, Genes kept: %d/%d (%.1f%%)", 
                  cell_type, 
                  sum(keep_genes), 
                  nrow(expr_matrix),
                  100 * sum(keep_genes)/nrow(expr_matrix)))
}

sapply(filtered_genes, length)
which(filtered_genes[['NK']]=='ENSG00000138378') 

##在enigma.cse.test中保留过滤后的基因
filtered_matrices <- list()
cell_types <- dimnames(enigma.cse.test)[[3]]

for (cell_type in cell_types) {
  all_genes <- dimnames(enigma.cse.test)[[1]]
  keep_genes <- filtered_genes[[cell_type]]
  keep_indices <- which(all_genes %in% keep_genes)
  filtered_mat <- enigma.cse.test[keep_indices, , cell_type]
  filtered_matrices[[cell_type]] <- filtered_mat
}

str(filtered_matrices)

##将过滤后的基因表达矩阵进行INT矩阵，并输出为predictDB和s-predixcan需要的格式
cell_types <- names(filtered_matrices)

for (ctname in cell_types) {
  ct_matrix <- filtered_matrices[[ctname]]
  #转置（样本×基因）
  tmp <- as.data.frame(t(ct_matrix))
  #INT转换
  tmp.INT <- matrix(0, nrow(tmp), ncol(tmp))
  for (i in 1:ncol(tmp)) {
    tmp.INT[, i] <- qqnorm(tmp[, i], plot.it = FALSE)$x
  }
  
  rownames(tmp.INT) <- rownames(tmp)
  colnames(tmp.INT) <- colnames(tmp)
  tmp.INT <- as.data.frame(tmp.INT)
  
  write.table(
    tmp.INT,
    file = paste0("/node1/liuxy/chplot/kidney_new/output/glom/count/tpm0.1sample0.01/", ctname, "_gene_exp.csv"),
    sep = ",",
    col.names = TRUE,
    row.names = FALSE
  )
  write.table(
    tmp.INT,
    file = paste0("/node1/liuxy/chplot/kidney_new/output/glom/count/tpm0.1sample0.01/", ctname, "_transformed_expression.txt"),
    sep = "\t",
    row.names = TRUE
  )
}

#tub
load("/node1/liuxy/chplot/kidney_new/output/tub/count/enigma.cse.trace_alpha_0.1.RData")
cell_types <- dimnames(enigma.cse.test)[[3]]
filtered_genes <- list()
for (cell_type in cell_types) {
  expr_matrix <- enigma.cse.test[, , cell_type]
  #每个基因在多少样本中表达 > 0.1
  n_samples_above_threshold <- rowSums(expr_matrix > 0.1)
  #至少需要的样本数 (1%),并且至少要有5个样本表达
  min_samples_1 <- ceiling(0.01 * ncol(expr_matrix))
  min_samples <- max(min_samples_1, 5)
  keep_genes <- n_samples_above_threshold >= min_samples
  
  filtered_genes[[cell_type]] <- rownames(expr_matrix)[keep_genes]
  message(sprintf("Cell type: %s, Genes kept: %d/%d (%.1f%%)", 
                  cell_type, 
                  sum(keep_genes), 
                  nrow(expr_matrix),
                  100 * sum(keep_genes)/nrow(expr_matrix)))
}

sapply(filtered_genes, length)
which(filtered_genes[['NK']]=='ENSG00000138378') 

##在enigma.cse.test中保留过滤后的基因
filtered_matrices <- list()
cell_types <- dimnames(enigma.cse.test)[[3]]

for (cell_type in cell_types) {
  all_genes <- dimnames(enigma.cse.test)[[1]]
  keep_genes <- filtered_genes[[cell_type]]
  keep_indices <- which(all_genes %in% keep_genes)
  filtered_mat <- enigma.cse.test[keep_indices, , cell_type]
  filtered_matrices[[cell_type]] <- filtered_mat
}
str(filtered_matrices)
cell_types <- names(filtered_matrices)
for (ctname in cell_types) {
  ct_matrix <- filtered_matrices[[ctname]]
  tmp <- as.data.frame(t(ct_matrix))
  tmp.INT <- matrix(0, nrow(tmp), ncol(tmp))
  for (i in 1:ncol(tmp)) {
    tmp.INT[, i] <- qqnorm(tmp[, i], plot.it = FALSE)$x
  }
  
  rownames(tmp.INT) <- rownames(tmp)
  colnames(tmp.INT) <- colnames(tmp)
  tmp.INT <- as.data.frame(tmp.INT)
  
  write.table(
    tmp.INT,
    file = paste0("/node1/liuxy/chplot/kidney_new/output/tub/count/tpm0.1sample0.01/", ctname, "_gene_exp.csv"),
    sep = ",",
    col.names = TRUE,
    row.names = FALSE
  )
  write.table(
    tmp.INT,
    file = paste0("/node1/liuxy/chplot/kidney_new/output/tub/count/tpm0.1sample0.01/", ctname, "_transformed_expression.txt"),
    sep = "\t",
    row.names = TRUE
  )
}

