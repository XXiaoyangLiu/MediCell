#r4.4.3 971个样本转换成cell fraction
library(data.table)
library(tidyr)
library(dplyr)
truth = readRDS("/node1/zhaoq/scTWAS_PBMC/CIBERSORTx/pbmc_7_CellType_971samples.rds")
frax = data.frame(s = truth$individual, ct = truth$cell_label)
frax = count(frax, s, ct)
frax = pivot_wider(frax, names_from = ct, values_from = n, values_fill = list(n = 0))
frax$s <- paste0("X", frax$s)
identical(frax$s, cibersort$Mixture)

rs <- rowSums(frax[,-1])
tmp <- lapply(frax[,-1], function(ct){ct/rs})
frax <- data.frame(frax[,1], tmp)
rowSums(frax[,-1])
write.table(frax, "/node1/liuxy/add_plot/cibersort_compare/usefile/frax_ground_truth.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#只用ENIGMA
load("/node1/liuxy/add_plot/enigma_frac_exp/output/myenigma_alpha0.1.RData")
enigma_direct_fractions <- as.data.frame(myenigma@result_cell_proportion)
enigma_direct_fractions$Mixture <- rownames(enigma_direct_fractions)
enigma_direct_fractions$Mixture <- gsub("^X", "", enigma_direct_fractions$Mixture)

#比较三种cell fraction(barplot)
library(data.table)
library(tidyr)
library(dplyr)
library(BayesPrism)
library(ggplot2)
setwd("/node1/liuxy/add_plot/cibersort_compare/output/pbmc/")
load("/node1/liuxy/add_plot/cibersort_compare/output/pbmc/pbmc_kidney_bayesprism_full_results.RData")
pbmc.theta <- get.fraction(bp = onek1k.res, which.theta = "final", state.or.type = "type")
pbmc.theta <- as.data.frame(pbmc.theta)
cibersort <- read.table("/node1/zhaoq/scTWAS_PBMC/CIBERSORTx/7CellType/result/CIBERSORTx_Job129_Adjusted_pbmc_54genes.txt",header = T)
truth <- read.table("/node1/liuxy/add_plot/cibersort_compare/usefile/frax_ground_truth.txt",header = T)

truth$Mixture <- truth$s
cibersort$Mixture <- cibersort$Mixture
pbmc.theta$Mixture <- rownames(pbmc.theta)
truth_long <- truth %>%
  select(Mixture, B, CD4T, CD8T, DC, Mono, NK, Plasma) %>%
  pivot_longer(cols = -Mixture, names_to = "celltype", values_to = "frac") %>%
  mutate(method = "Ground Truth")
cbs_long <- cibersort %>%
  select(Mixture, B, CD4T, CD8T, DC, Mono, NK, Plasma) %>%
  pivot_longer(cols = -Mixture, names_to = "celltype", values_to = "frac") %>%
  mutate(method = "CIBERSORTx")
bp_long <- pbmc.theta %>%
  select(Mixture, B, CD4T, CD8T, DC, Mono, NK, Plasma) %>%
  pivot_longer(cols = -Mixture, names_to = "celltype", values_to = "frac") %>%
  mutate(method = "BayesPrism")
all_data <- bind_rows(truth_long, cbs_long, bp_long)
mean_frac <- all_data %>%
  group_by(celltype, method) %>%
  summarise(mean_frac = mean(frac), .groups = "drop")
celltypes <- c("CD4T","CD8T","NK","Plasma","B","Mono","DC")

# CIBERSORTx vs Ground Truth
cor_cbs <- sapply(celltypes, function(ct) {
  cor(truth[[ct]], cibersort[[ct]], use = "complete.obs")
})
# BayesPrism vs Ground Truth  
cor_bp <- sapply(celltypes, function(ct) {
  cor(truth[[ct]], pbmc.theta[[ct]], use = "complete.obs")
})

cor_df <- data.frame(
  celltype = rep(celltypes, 2),
  method = rep(c("CIBERSORTx", "BayesPrism"), each = length(celltypes)),
  correlation = c(cor_cbs, cor_bp)
)
plot_data <- mean_frac %>%
  left_join(cor_df, by = c("celltype", "method")) %>%
  mutate(
    label = ifelse(!is.na(correlation), 
                   sprintf("r = %.3f", correlation), 
                   ""),
    correlation = ifelse(method == "Ground Truth", NA, correlation),
    label = ifelse(method == "Ground Truth", "Ground Truth", label)
  )
plot_data$celltype <- factor(plot_data$celltype, levels = celltypes)
plot_data$method <- factor(plot_data$method, 
                          levels = c("Ground Truth", "CIBERSORTx", "BayesPrism"))

pdf("three_methods_comparison_horizontal.pdf", width = 12, height = 8)
ggplot(plot_data, aes(x = mean_frac, y = celltype, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = label, x = mean_frac + 0.02), 
            position = position_dodge(width = 0.8),
            hjust = 0, size = 4, fontface = "bold") +
  scale_fill_manual(values = c(
    "Ground Truth" = "#4DAF4A",
    "CIBERSORTx" = "#E41A1C", 
    "BayesPrism" = "#377EB8"
  )) +
  labs(
    x = "Mean Cell Fraction",
    y = "Cell Type",
    title = "Comparison of Cell Fraction Estimation Methods",
    fill = "Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 11)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)))
dev.off()

# 比较三种cell fraction
library(data.table)
library(tidyr)
library(dplyr)
library(BayesPrism)
library(ggplot2)
setwd("/node1/liuxy/add_plot/cibersort_compare/output/pbmc/")
load("/node1/liuxy/add_plot/cibersort_compare/output/pbmc/pbmc_kidney_bayesprism_full_results.RData")
pbmc.theta <- get.fraction(bp = onek1k.res, which.theta = "final", state.or.type = "type")
pbmc.theta <- as.data.frame(pbmc.theta)
cibersort <- read.table("/node1/zhaoq/scTWAS_PBMC/CIBERSORTx/7CellType/result/CIBERSORTx_Job129_Adjusted_pbmc_54genes.txt",header = T)
truth <- read.table("/node1/liuxy/add_plot/cibersort_compare/usefile/frax_ground_truth.txt",header = T)
enigma <- read.table("/node1/liuxy/add_plot/enigma_frac_exp/output/cell_fractions_alpha_0.1.txt", header = T)

truth$Mixture <- truth$s
cibersort$Mixture <- cibersort$Mixture
pbmc.theta$Mixture <- rownames(pbmc.theta)
enigma$Mixture <- rownames(enigma)

truth_long <- truth %>%
  select(Mixture, B, CD4T, CD8T, DC, Mono, NK, Plasma) %>%
  pivot_longer(cols = -Mixture, names_to = "celltype", values_to = "frac") %>%
  mutate(method = "Ground Truth")

cbs_long <- cibersort %>%
  select(Mixture, B, CD4T, CD8T, DC, Mono, NK, Plasma) %>%
  pivot_longer(cols = -Mixture, names_to = "celltype", values_to = "frac") %>%
  mutate(method = "CIBERSORTx")

bp_long <- pbmc.theta %>%
  select(Mixture, B, CD4T, CD8T, DC, Mono, NK, Plasma) %>%
  pivot_longer(cols = -Mixture, names_to = "celltype", values_to = "frac") %>%
  mutate(method = "BayesPrism")

enigma_long <- enigma %>%
  select(Mixture, B, CD4T, CD8T, DC, Mono, NK, Plasma) %>%
  pivot_longer(cols = -Mixture, names_to = "celltype", values_to = "frac") %>%
  mutate(method = "ENIGMA")

all_data <- bind_rows(truth_long, cbs_long, enigma_long, bp_long)
all_data$celltype <- factor(all_data$celltype, levels = celltypes)
all_data$method <- factor(all_data$method, 
                         levels = c("Ground Truth", "CIBERSORTx", "ENIGMA", "BayesPrism"))

cor_cbs <- sapply(celltypes, function(ct) {
  cor(truth[[ct]], cibersort[[ct]], use = "complete.obs")
})

cor_bp <- sapply(celltypes, function(ct) {
  cor(truth[[ct]], pbmc.theta[[ct]], use = "complete.obs")
})

cor_enigma <- sapply(celltypes, function(ct) {
  cor(truth[[ct]], enigma[[ct]], use = "complete.obs")
})

cor_df <- data.frame(
  celltype = rep(celltypes, 3),
  method = rep(c("CIBERSORTx", "ENIGMA", "BayesPrism"), each = length(celltypes)),
  correlation = c(cor_cbs, cor_enigma, cor_bp)
)

#箱线图
pdf("four_methods_comparison_boxplot_2.pdf", width = 13, height = 6)
ggplot(all_data, aes(x = celltype, y = frac, fill = method)) +
  geom_boxplot(position = position_dodge(0.9), alpha = 0.7, outlier.shape = NA) +
  geom_text(data = cor_df, 
            aes(x = celltype, y = max(all_data$frac, na.rm = TRUE) * 0.95, 
                label = sprintf("r=%.3f", correlation), group = method),
            position = position_dodge(0.9),
            vjust = 0, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c(
    "Ground Truth" = "#916042",
    "CIBERSORTx" = "#f67d25", 
    "ENIGMA" = "#ffe8d1",
    "BayesPrism" = "#005246"
  )) +
  labs(
    x = "Cell Type",
    y = "Cell Fraction",
    title = "Comparison of Cell Fraction Estimation Methods",
    fill = "Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1)))
dev.off()
write.csv(all_data, "cell_fraction_boxplot_data.csv", row.names = FALSE)

#表达谱
library(data.table)
library(tidyr)
library(dplyr)
library(BayesPrism)
library(ggplot2)
setwd("/node1/liuxy/MediCell/method_compare/deconvolution/bayesprism/")
#保存exp
load("/node1/liuxy/MediCell/method_compare/deconvolution/bayesprism/output/pbmc_bayesprism_full_results.RData")
pbmc.exp <- onek1k.res@posterior.initial.cellState@Z
gene_ref <- read.csv("/node1/liuxy/index/gene_pos_grch37.82.csv")
pbmc_genes <- dimnames(pbmc.exp)[[2]]
coding_gene <- intersect(pbmc_genes,gene_ref$gene_id[gene_ref$gene_type %in% "protein_coding"])

celltypes <- dimnames(pbmc.exp)[[3]]
expression_list <- list()
for(ct in celltypes) {
  expr_matrix <- t(pbmc.exp[, coding_gene, ct])
  expression_list[[ct]] <- expr_matrix
  output_df <- data.frame(Gene = rownames(expr_matrix), expr_matrix)
  filename <- paste0("/node1/liuxy/MediCell/method_compare/deconvolution/bayesprism/output/",ct, "_expression.txt")
  write.table(output_df, filename, sep = "\t", row.names = FALSE, quote = FALSE)
}
#groundtruth
library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(pheatmap)
library(ggplot2)
gene_ref <- read.csv("/node1/liuxy/index/gene_pos_grch37.82.csv")
load("/node1/liuxy/chplot/combind_plasma_B/enigma/output/scTPM971_direct.RData")
pbmc_genes <- rownames(combined_data)
coding_gene <- intersect(pbmc_genes,gene_ref$gene_id[gene_ref$gene_type %in% "protein_coding"])
sample_names_full <- colnames(combined_data)
sample_names <- unique(sapply(strsplit(sample_names_full, "_"), function(x) {
  paste(x[1:(length(x)-1)], collapse = "_")
}))
celltypes <- unique(sapply(strsplit(sample_names_full, "_"), function(x) {
  tail(x, 1)
}))

expression_list <- list()
for(ct in celltypes) {
    ct_columns <- grep(paste0("_", ct, "$"), sample_names_full, value = TRUE)
    expr_matrix <- combined_data[coding_gene, ct_columns, drop = FALSE]
    new_colnames <- sapply(strsplit(colnames(expr_matrix), "_"), function(x) {
      paste(head(x, -1), collapse = "_")
    })
    colnames(expr_matrix) <- new_colnames
    ct_df <- data.frame(Gene = rownames(expr_matrix), expr_matrix)
    expression_list[[ct]] <- ct_df
    filename <- paste0("/node1/liuxy/MediCell/method_compare/deconvolution/gdt/", ct, "_ground_truth_expression.txt")
    write.table(ct_df, filename, sep = "\t", row.names = FALSE, quote = FALSE)
}

#cibersort+ENIGMA
celltypes <- c("CD4T","CD8T","NK","B","Mono","DC")
load("/node1/liuxy/chplot/combind_plasma_B/enigma/output/count/enigma.cse.trace_norm_alpha_0.1.RData")
boolvec = str_detect(dimnames(enigma.cse.test)[[2]], "^X1_")
tmp <- ifelse(boolvec, gsub("^.", "", dimnames(enigma.cse.test)[[2]]), dimnames(enigma.cse.test)[[2]])
dimnames(enigma.cse.test)[[2]] <- tmp
expression_enigma_corrected <- list()
celltypes <- c("CD4T","CD8T","NK","B","Mono","DC")
for(ct in celltypes) {
  enigma_infered <- enigma.cse.test[,,ct]
  colnames(enigma_infered) <- str_remove(colnames(enigma_infered), "^X")
  expression_enigma_corrected[[ct]] <- enigma_infered
  filename <- paste0("/node1/liuxy/add_plot/cibersort_compare/usefile/exp/cibersort_enigma/", ct, "_cb_en_expression.txt")
  write.table(enigma_infered, filename, sep = "\t", row.names = TRUE, quote = FALSE)
}
expression_enigma <- expression_enigma_corrected



#画图
setwd("/node1/liuxy/MediCell/method_compare/deconvolution/plot/")
#读cibersort+ENIGMA
cben_path <- "/node1/liuxy/MediCell/method_compare/deconvolution/cibersort_enigma/"
expression_cben <- list()
for(ct in celltypes) {
  filename <- paste0(cben_path, ct, "_cb_en_expression.txt")
    cben_data <- read.delim(filename, sep = "\t", check.names = FALSE)
    expression_cben[[ct]] <- as.matrix(cben_data)
}
#读ENIGMA
enigma_path <- "/node1/liuxy/MediCell/method_compare/deconvolution/one_step_enigma/output/"
expression_enigma <- list()
for(ct in celltypes) {
  filename <- paste0(enigma_path, "enigma_expression_", ct,".txt")
    enigma_data <- read.delim(filename, sep = "\t", check.names = FALSE)
    rownames(enigma_data) <- enigma_data[, 1]
    enigma_data <- enigma_data[, -1, drop = FALSE]
    expression_enigma[[ct]] <- as.matrix(enigma_data)
}
#读BayesPrism
bp_path <- "/node1/liuxy/MediCell/method_compare/deconvolution/bayesprism/output/"
expression_bayesprism <- list()
for(ct in celltypes) {
  filename <- paste0(bp_path, ct, "_expression.txt")
    bp_data <- read.delim(filename, sep = "\t", check.names = FALSE)
    rownames(bp_data) <- bp_data[, 1]
    bp_data <- bp_data[, -1, drop = FALSE]
    expression_bayesprism[[ct]] <- as.matrix(bp_data)
}
#读ground truth
gt_path <- "/node1/liuxy/MediCell/method_compare/deconvolution/gdt/"
expression_truth <- list()
for(ct in celltypes) {
  filename <- paste0(gt_path, ct, "_ground_truth_expression.txt")
    gt_data <- read.delim(filename, sep = "\t", check.names = FALSE)
    rownames(gt_data) <- gt_data[, 1]
    gt_data <- gt_data[, -1, drop = FALSE]
    expression_truth[[ct]] <- as.matrix(gt_data)
}

calculate_sample_level_correlation_corrected <- function(truth_mat, method_mat, celltype, method_name) {
  if(nrow(method_mat) < ncol(method_mat)) {
    method_mat <- t(method_mat)
  }
  common_genes <- intersect(rownames(truth_mat), rownames(method_mat))

  common_samples <- intersect(colnames(truth_mat), colnames(method_mat))
  if(length(common_genes) == 0 || length(common_samples) == 0) {
    return(NULL)
  }
  truth_common <- truth_mat[common_genes, common_samples, drop = FALSE]
  method_common <- method_mat[common_genes, common_samples, drop = FALSE]
  
  truth_samples <- t(truth_common)
  method_samples <- t(method_common)
  sample_cors <- numeric(nrow(truth_samples))
  for(i in 1:nrow(truth_samples)) {
    sample_cors[i] <- cor(truth_samples[i, ], method_samples[i, ], 
                         method = "pearson", use = "complete.obs")
  }
  names(sample_cors) <- rownames(truth_samples)
  list(
    celltype = celltype,
    method = method_name,
    sample_correlations = sample_cors,
    summary = data.frame(
      CellType = celltype,
      Method = method_name,
      Median_Correlation = median(sample_cors, na.rm = TRUE),
      IQR = IQR(sample_cors, na.rm = TRUE),
      Mean_Correlation = mean(sample_cors, na.rm = TRUE),
      N_Samples = length(common_samples),
      N_Genes = length(common_genes)
    )
  )
}

all_results_corrected <- list()
performance_summary_corrected <- data.frame()
for(ct in celltypes) {
  if(ct %in% names(expression_truth) && 
     ct %in% names(expression_bayesprism) && 
     ct %in% names(expression_enigma)) {
    truth_mat <- expression_truth[[ct]]
    bp_mat <- expression_bayesprism[[ct]]
    enigma_mat <- expression_enigma[[ct]]
    truth_samples <- str_remove(colnames(truth_mat), "^X")
    bp_samples <- str_remove(colnames(bp_mat), "^X")
    enigma_samples <- str_remove(colnames(enigma_mat), "^X")
    colnames(truth_mat) <- truth_samples
    colnames(bp_mat) <- bp_samples
    colnames(enigma_mat) <- enigma_samples
    common_samples <- intersect(truth_samples, enigma_samples)
    common_samples <- intersect(common_samples, bp_samples)
    
    if(length(common_samples) > 0) {
      truth_common <- truth_mat[, common_samples, drop = FALSE]
      bp_common <- bp_mat[, common_samples, drop = FALSE]
      enigma_common <- enigma_mat[, common_samples, drop = FALSE]
      bp_sample_cors <- numeric(length(common_samples))
      enigma_sample_cors <- numeric(length(common_samples))
      
      for(i in 1:length(common_samples)) {
        sample <- common_samples[i]
        common_genes_bp <- intersect(rownames(truth_common), rownames(bp_common))
        common_genes_enigma <- intersect(rownames(truth_common), rownames(enigma_common))
        if(length(common_genes_bp) > 10) {
          bp_sample_cors[i] <- cor(truth_common[common_genes_bp, sample], 
                                  bp_common[common_genes_bp, sample], 
                                  use = "complete.obs")
        }
        if(length(common_genes_enigma) > 10) {
          enigma_sample_cors[i] <- cor(truth_common[common_genes_enigma, sample], 
                                      enigma_common[common_genes_enigma, sample], 
                                      use = "complete.obs")
        }
      }
      performance_summary_corrected <- rbind(performance_summary_corrected,
        data.frame(
          CellType = ct,
          Method = "BayesPrism",
          Mean_Correlation = mean(bp_sample_cors, na.rm = TRUE),
          Median_Correlation = median(bp_sample_cors, na.rm = TRUE),
          N_Samples = sum(!is.na(bp_sample_cors))
        ),
        data.frame(
          CellType = ct,
          Method = "ENIGMA",
          Mean_Correlation = mean(enigma_sample_cors, na.rm = TRUE),
          Median_Correlation = median(enigma_sample_cors, na.rm = TRUE),
          N_Samples = sum(!is.na(enigma_sample_cors))
        )
      )
    }
  }
}
print(performance_summary_corrected)

heatmap_data <- performance_summary_corrected %>%
  select(CellType, Method, Mean_Correlation) %>%
  pivot_wider(names_from = Method, values_from = Mean_Correlation) %>%
  as.data.frame()

rownames(heatmap_data) <- heatmap_data$CellType
heatmap_matrix <- as.matrix(heatmap_data[, -1])
celltype_order <- c("CD4T","CD8T","NK","B","Mono","DC")
heatmap_matrix <- heatmap_matrix[celltype_order, ]

write.csv(performance_summary_corrected, "method_performance_summary_final.csv", row.names = FALSE)

boxplot_data <- data.frame()
for(ct in celltypes) {
  if(ct %in% names(expression_truth) && 
     ct %in% names(expression_bayesprism) && 
     ct %in% names(expression_cben) &&
     ct %in% names(expression_enigma)) {
    
    truth_mat <- expression_truth[[ct]]
    bp_mat <- expression_bayesprism[[ct]]
    cben_mat <- expression_cben[[ct]]
    enigma_mat <- expression_enigma[[ct]]
    
    truth_samples <- str_remove(colnames(truth_mat), "^X")
    bp_samples <- str_remove(colnames(bp_mat), "^X")
    enigma_samples <- str_remove(colnames(enigma_mat), "^X")
    cben_samples <- colnames(cben_mat)
    
    colnames(truth_mat) <- truth_samples
    colnames(bp_mat) <- bp_samples
    colnames(enigma_mat) <- enigma_samples
    
    common_samples <- intersect(truth_samples, enigma_samples)
    common_samples <- intersect(common_samples, bp_samples)
    common_samples <- intersect(common_samples, cben_samples)
    
    if(length(common_samples) > 0) {
      truth_common <- truth_mat[, common_samples, drop = FALSE]
      bp_common <- bp_mat[, common_samples, drop = FALSE]
      enigma_common <- enigma_mat[, common_samples, drop = FALSE]
      cben_common <- cben_mat[, common_samples, drop = FALSE]
      
      # 为每个样本计算相关系数
      for(sample in common_samples) {
        common_genes_bp <- intersect(rownames(truth_common), rownames(bp_common))
        common_genes_enigma <- intersect(rownames(truth_common), rownames(enigma_common))
        common_genes_cben <- intersect(rownames(truth_common), rownames(cben_common))

        if(length(common_genes_bp) > 10) {
          bp_cor <- cor(truth_common[common_genes_bp, sample], 
                       bp_common[common_genes_bp, sample], 
                       use = "complete.obs")
          boxplot_data <- rbind(boxplot_data,
                               data.frame(CellType = ct,
                                         Method = "BayesPrism",
                                         Correlation = bp_cor,
                                         Sample = sample))
        }
        
        if(length(common_genes_enigma) > 10) {
          enigma_cor <- cor(truth_common[common_genes_enigma, sample], 
                           enigma_common[common_genes_enigma, sample], 
                           use = "complete.obs")
          boxplot_data <- rbind(boxplot_data,
                               data.frame(CellType = ct,
                                         Method = "ENIGMA",
                                         Correlation = enigma_cor,
                                         Sample = sample))
        }

        if(length(common_genes_cben) > 10) {
          cben_cor <- cor(truth_common[common_genes_cben, sample], 
                       cben_common[common_genes_cben, sample], 
                       use = "complete.obs")
          boxplot_data <- rbind(boxplot_data,
                               data.frame(CellType = ct,
                                         Method = "Cibersort + ENIGMA",
                                         Correlation = cben_cor,
                                         Sample = sample))
        }
      }
    }
  }
}
celltype_order <- c("CD4T","CD8T","NK","B","Mono","DC")
boxplot_data$CellType <- factor(boxplot_data$CellType, levels = celltype_order)
pdf("final_correlation_boxplot_3.pdf", width = 6, height = 4)
ggplot(boxplot_data, aes(x = CellType, y = Correlation, fill = Method)) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.7, outlier.shape = NA) +
  scale_fill_manual(values = c("BayesPrism" = "#005246", "ENIGMA" = "#ffe8d1", "Cibersort + ENIGMA" = "#f67d25")) +
  labs(title = "Sample-level Expression Correlation with Ground Truth",
       x = "Cell Type", 
       y = "Pearson Correlation") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "top") +
  ylim(0.7, 1.05)
dev.off()
write.csv(boxplot_data, "correlation_boxplot_data.csv", row.names = FALSE)

#每个基因ground truth vs pred表达的相关性
library(ggplot2)
library(reshape2)
library(patchwork)
library(viridis)

target_genes <- c("IRF7", "BTN2A1", "ACAP1", "FAM167A", "KCTD11", 
                  "ORMDL3", "IRF5", "ZGLP1", "TMEM80", "BTN2A2", "BLK")

calc_mean_expression <- function(expr_list, method_name) {
  mean_expr <- data.frame()
  for(ct in names(expr_list)) {
    expr_matrix <- expr_list[[ct]]
    for(gene in target_genes) {
      if(gene %in% rownames(expr_matrix)) {
        values <- as.numeric(expr_matrix[gene, ])
        values <- values[!is.na(values) & is.finite(values)]
        if(length(values) > 0) {
          mean_val <- mean(values)
          if(mean_val == 0) mean_val <- 1e-10
          mean_expr <- rbind(mean_expr, 
                            data.frame(gene = gene,
                                      celltype = ct,
                                      mean_expression = mean_val,
                                      log2_expr = log2(mean_val + 1),
                                      method = method_name))
        }
      }
    }
  }
  return(mean_expr)
}

mean_cben <- calc_mean_expression(expression_cben, "CIBERSORT+ENIGMA")
mean_gt <- calc_mean_expression(expression_truth, "Ground Truth")

# 合并数据
combined_expr <- rbind(mean_cben, mean_gt)

# 重新排序celltype
celltype_order <- c("CD4T","CD8T","NK","B","Plasma","Mono","DC")
combined_expr$celltype <- factor(combined_expr$celltype, levels = celltype_order)

pdf("gene_exp_heatmap_gt_cben.pdf", width = 12, height = 8)
ggplot(combined_expr, aes(x = celltype, y = gene, fill = log2_expr)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(
    low = "#2166AC",
    high = "#B2182B",
    midpoint = median(combined_expr$log2_expr, na.rm = TRUE),
    name = "Log2(Mean Expr)"
  ) +
  facet_grid(. ~ method) +
  labs(title = "Mean Expression Comparison: Target Genes",
       x = "Cell Type", y = "Gene") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 9, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.grid = element_blank(),
    panel.spacing = unit(0.5, "cm")
  )
dev.off()
write.csv(combined_expr, "gene_exp_gt_cben_heatmap_data.csv", row.names = FALSE)

#皮尔森相关系数
calc_celltype_gene_cor <- function() {
  cor_results <- data.frame()
  for(ct in celltypes) {
    cben_mat <- expression_cben[[ct]]
    gt_mat <- expression_truth[[ct]]
    cben_cols <- colnames(cben_mat)
    cben_samples <- sapply(strsplit(cben_cols, "_"), function(x) paste(x, collapse = "_"))
    gt_cols <- colnames(gt_mat)
    gt_samples <- sapply(gt_cols, function(x) {
      if(startsWith(x, "X")) {
        return(substr(x, 2, nchar(x)))
      } else {
        return(x)
      }
    })
    common_samples <- intersect(cben_samples, gt_samples)
    if(length(common_samples) > 3) {
      for(gene in target_genes) {
        if(gene %in% rownames(cben_mat) && gene %in% rownames(gt_mat)) {
          cben_values <- numeric(length(common_samples))
          gt_values <- numeric(length(common_samples))
          for(i in seq_along(common_samples)) {
            sample <- common_samples[i]
            cben_idx <- which(cben_samples == sample)
            gt_idx <- which(gt_samples == sample)
            
            if(length(cben_idx) > 0 && length(gt_idx) > 0) {
              cben_values[i] <- as.numeric(cben_mat[gene, cben_idx[1]])
              gt_values[i] <- as.numeric(gt_mat[gene, gt_idx[1]])
            }
          }
          valid_idx <- !is.na(cben_values) & !is.na(gt_values) & 
                       is.finite(cben_values) & is.finite(gt_values)
          
          if(sum(valid_idx) > 3 && sd(cben_values[valid_idx]) > 0 && sd(gt_values[valid_idx]) > 0) {
            pearson_cor <- cor(cben_values[valid_idx], gt_values[valid_idx], 
                               method = "pearson")
            cor_test <- cor.test(cben_values[valid_idx], gt_values[valid_idx])
            p_value <- cor_test$p.value
            sig_symbol <- ""
            if(p_value < 0.001) {
              sig_symbol <- "***"
            } else if(p_value < 0.01) {
              sig_symbol <- "**"
            } else if(p_value < 0.05) {
              sig_symbol <- "*"
            }
            cor_results <- rbind(cor_results, 
                                data.frame(gene = gene,
                                          celltype = ct,
                                          pearson_cor = pearson_cor,
                                          p_value = p_value,
                                          sig_symbol = sig_symbol,
                                          n_samples = sum(valid_idx)))
          }
        }
      }
    }
  }
  return(cor_results)
}
cor_results <- calc_celltype_gene_cor()
combined_expr_with_cor <- merge(combined_expr, cor_results, 
                                by = c("gene", "celltype"), 
                                all.x = TRUE)

combined_expr_with_cor$cor_label <- ifelse(
  !is.na(combined_expr_with_cor$pearson_cor),
  sprintf("%.2f", combined_expr_with_cor$pearson_cor),
  NA
)
pdf("gene_exp_heatmap_gt_cben_with_correlation.pdf", width = 14, height = 8)
ggplot(combined_expr_with_cor, aes(x = celltype, y = gene, fill = log2_expr)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = cor_label), 
            size = 3, 
            color = "black",
            fontface = "bold",
            na.rm = TRUE) +
  scale_fill_gradient2(
    low = "#2166AC",
    high = "#B2182B",
    midpoint = median(combined_expr_with_cor$log2_expr, na.rm = TRUE),
    name = "Log2(Mean Expr)"
  ) +
  facet_grid(. ~ method) +
  labs(x = "Cell Type", y = "Gene") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 9, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    panel.grid = element_blank(),
    panel.spacing = unit(0.5, "cm")
  )
dev.off()
