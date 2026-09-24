#TWAS filt
library(stringr)
library(dplyr)
library(tidyr)
library(ggVennDiagram)
library(ggplot2)
library(reshape2)
library(patchwork)
library(UpSetR)
library(grid)
library(data.table)
library(dendextend)

gene_grch37 <- read.csv("/node1/liuxy/index/gene_pos_grch37.82.csv")
tissue_dirs <- c(
  pbmc = "/node1/liuxy/MediCell/method_compare_2/output/pseudobulk_TWAS/pbmc/SLE/"
)
summary_df <- data.frame(
  tissue = character(),
  gwas_id = character(),
  total_genes = integer(),
  fdr_sig = integer(),
  bonf_sig = integer(),
  stringsAsFactors = FALSE
)
for (tissue_name in names(tissue_dirs)) {
  TWAS_dir <- tissue_dirs[tissue_name]
  csv_files <- list.files(TWAS_dir, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(csv_files) > 0) {
    gwas_id <- "SLE" 
    output_dir <- file.path(TWAS_dir, "summary")
    cat("Processing tissue:", tissue_name, "\n")
    cat("GWAS:", gwas_id, "\n")
    cat("Directory:", TWAS_dir, "\n")
    
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      cat("Created:", output_dir, "\n")
    } else {
      cat("Directory exists:", output_dir, "\n")
    }
    peer_value <- case_when(
      tissue_name == "pbmc" ~ "peer2"
    )
    df_list <- list()
    for (file in csv_files) {
      tryCatch({
        dat <- read.csv(file)  
        filename <- basename(file) 
        cell_type <- tools::file_path_sans_ext(filename)   
        ngeneindb <- nrow(dat)  
  #      bonferroni_p <- 0.05 / ngeneindb  
        dat$fdr_p <- p.adjust(dat$pvalue, method = "BH")
        dat$bonferroni_p <- p.adjust(dat$pvalue, method = "bonferroni")
        dat$cell_type <- cell_type  
        dat$ngeneindb <- ngeneindb  
 #       dat$bonferroni_p <- bonferroni_p  
        dat$gwas_id <- gwas_id
        dat$tissue <- tissue_name
        dat$peer <- peer_value
        dat$gene_name <- gene_grch37$gene_name[match(dat$gene, gene_grch37$gene_id)]
        df_list[[filename]] <- dat
        cat("Processed:", filename, "\n")
      }, error = function(e) {
        cat("File error:", file, "\nError:", e$message, "\n")
      })
    }
    
    if (length(df_list) == 0) {
      cat("Warning: No valid CSV files processed for", tissue_name, "\n")
    } else {
      final_df <- bind_rows(df_list) 
      colnames(final_df)[1:2] <- c('gene_id', 'gene_name')
      total_genes <- nrow(final_df)
    #  final_df$bonferroni_p <- 0.05 / total_genes
      final_df$fdr_p <- p.adjust(final_df$pvalue, method = "BH")
      final_df$bonferroni_p <- p.adjust(final_df$pvalue, method = "bonferroni")
      bonf_sig_df <- subset(final_df, bonferroni_p < 0.05)
      fdr_sig_df <- subset(final_df, fdr_p < 0.05)
      all_sig_df <- unique(rbind(bonf_sig_df, fdr_sig_df))
      final_df$FDR.sig <- ifelse(final_df$fdr_p < 0.05, 'TRUE', 'FALSE')
      final_df$bonferroni_p.sig <- ifelse(final_df$bonferroni_p < 0.05, 'TRUE', 'FALSE')
      cat("\n统计信息 -", tissue_name, ":\n")
      cat("总基因数量:", nrow(final_df), "\n")
      cat("FDR显著基因数量:", nrow(fdr_sig_df), "\n")
      cat("Bonferroni显著基因数量:", nrow(bonf_sig_df), "\n")
      cat("\n按细胞类型统计:\n")
      write.csv(final_df, file.path(output_dir, "allcelltype_result.csv"), row.names = FALSE) 
      write.csv(bonf_sig_df, file.path(output_dir, "Bonfsig_allcelltype_result.csv"), row.names = FALSE) 
      write.csv(fdr_sig_df, file.path(output_dir, "FDRsig_allcelltype_result.csv"), row.names = FALSE)
      write.csv(all_sig_df, file.path(output_dir, "BonFDRsig_allcelltype_result.csv"), row.names = FALSE)

      cat("\nOutput saved to:", output_dir, "\n")
      summary_df <- rbind(summary_df, data.frame(
        tissue = tissue_name,
        gwas_id = gwas_id,
        total_genes = nrow(final_df),
        fdr_sig = nrow(fdr_sig_df),
        bonf_sig = nrow(bonf_sig_df)
      ))
    }
  } else {
    cat("No CSV files found in:", TWAS_dir, "\n")
  }
}

#scpredixcan
tissue_dirs <- c(pbmc = "/node1/liuxy/MediCell/method_compare_2/output/scpredixcan_TWAS/scpredixcan/")
summary_df <- data.frame(
  tissue = character(),
  gwas_id = character(),
  total_genes = integer(),
  fdr_sig = integer(),
  bonf_sig = integer(),
  stringsAsFactors = FALSE
)
for (tissue_name in names(tissue_dirs)) {
  TWAS_dir <- tissue_dirs[tissue_name]
  csv_files <- list.files(TWAS_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) > 0) {
    gwas_id <- "SLE" 
    output_dir <- file.path(TWAS_dir, "summary")
    cat("Processing tissue:", tissue_name, "\n")
    cat("GWAS:", gwas_id, "\n")
    cat("Directory:", TWAS_dir, "\n")
    
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      cat("Created:", output_dir, "\n")
    } else {
      cat("Directory exists:", output_dir, "\n")
    }
    peer_value <- case_when(
      tissue_name == "pbmc" ~ "None"
    )
    df_list <- list()
    for (file in csv_files) {
      tryCatch({
        dat <- read.csv(file)  
        filename <- basename(file) 
        cell_type <- tools::file_path_sans_ext(filename)   
        ngeneindb <- nrow(dat)  
        bonferroni_p <- 0.05 / ngeneindb  
        dat$fdr_p <- p.adjust(dat$pvalue, method = "BH")
        dat$cell_type <- cell_type  
        dat$ngeneindb <- ngeneindb  
        dat$bonferroni_p <- bonferroni_p  
        dat$gwas_id <- gwas_id
        dat$tissue <- tissue_name
        dat$peer <- peer_value
        df_list[[filename]] <- dat
        cat("Processed:", filename, "\n")
      }, error = function(e) {
        cat("File error:", file, "\nError:", e$message, "\n")
      })
    }
    
    if (length(df_list) == 0) {
      cat("Warning: No valid CSV files processed for", tissue_name, "\n")
    } else {
      final_df <- bind_rows(df_list) 
      colnames(final_df)[1:2] <- c('gene_name', 'gene_id')
      
      bonf_sig_df <- subset(final_df, pvalue < bonferroni_p)
      fdr_sig_df <- subset(final_df, fdr_p < 0.05)
      all_sig_df <- unique(rbind(bonf_sig_df, fdr_sig_df))
      final_df$FDR.sig <- ifelse(final_df$fdr_p < 0.05, 'TRUE', 'FALSE')
      cat(tissue_name, ":\n")
      cat(nrow(final_df), "\n")
      cat(nrow(fdr_sig_df), "\n")
      cat(nrow(bonf_sig_df), "\n")
      write.csv(final_df, file.path(output_dir, "allcelltype_result.csv"), row.names = FALSE) 
      write.csv(bonf_sig_df, file.path(output_dir, "Bonfsig_allcelltype_result.csv"), row.names = FALSE)
      write.csv(fdr_sig_df, file.path(output_dir, "FDRsig_allcelltype_result.csv"), row.names = FALSE)
      write.csv(all_sig_df, file.path(output_dir, "BonFDRsig_allcelltype_result.csv"), row.names = FALSE)
      
      cat("\nOutput saved to:", output_dir, "\n")
      summary_df <- rbind(summary_df, data.frame(
        tissue = tissue_name,
        gwas_id = gwas_id,
        total_genes = nrow(final_df),
        fdr_sig = nrow(fdr_sig_df),
        bonf_sig = nrow(bonf_sig_df)
      ))
    }
  } else {
    cat("No CSV files found in:", TWAS_dir, "\n")
  }
}

#cpmpare plot
setwd("/node1/liuxy/MediCell/method_compare_2/output/plot/")
celltypes <- c("CD4T","CD8T","NK","B","Mono","DC")
gene_pos_df <- read.csv("/node1/liuxy/index/gene_pos_grch37.82.csv")
mhc_genes <- gene_pos_df$gene_name[gene_pos_df$chr %in% "6" & as.numeric(gene_pos_df$start) >= 28477797 & as.numeric(gene_pos_df$end) <= 33448354]
TWAS_Pred <- read.csv("/node1/liuxy/MediCell/pbmc_2/TWAS/output/pbmc/SLE/summary/Bonfsig_allcelltype_result.csv")
twas_pred_filtered <- TWAS_Pred[!TWAS_Pred$gene_name %in% mhc_genes, ]
#twas_pred_filtered <- twas_pred_filtered[twas_pred_filtered$cell_type %in% celltypes, ]
TWAS_ct <- read.csv("/node1/liuxy/MediCell/method_compare_2/output/pseudobulk_TWAS/pbmc/SLE/summary/Bonfsig_allcelltype_result.csv")
TWAS_ct_filtered <- TWAS_ct[!TWAS_ct$gene_name %in% mhc_genes, ]
celltype_mapping <- list(
  "DC" = c("DC"),
  "Mono" = c("MonoC","MonoNC"),
  "NK" = c("NKmat","NKact"),
  "CD4T" = c("CD4all","CD4effCM","CD4TGFbStim"),
  "CD8T" = c("CD8all","CD8eff","CD8unknown"),
  "B" = c("BimmNaive", "Bmem", "Plasma")
)
celltype_mapping_df <- data.frame(
  detailed_celltype = unlist(celltype_mapping),
  main_celltype = rep(names(celltype_mapping), 
                      sapply(celltype_mapping, length)),
  stringsAsFactors = FALSE
)
TWAS_ct_filtered$cell_type_14 <- TWAS_ct_filtered$cell_type
TWAS_ct_filtered$cell_type <- celltype_mapping_df$main_celltype[
  match(TWAS_ct_filtered$cell_type_14, celltype_mapping_df$detailed_celltype)
]

scpredixcan <- read.csv("/node1/liuxy/MediCell/method_compare_2/output/scpredixcan_TWAS/scpredixcan/summary/Bonfsig_allcelltype_result.csv")
colnames(scpredixcan)[colnames(scpredixcan) %in% c("gene_name","gene_id","cell_type")] <- c("ENSG","gene_name","cell_type_21")
scpredixcan_filtered <- scpredixcan[!scpredixcan$gene_name %in% mhc_genes, ]
scpredixcan_filtered <- scpredixcan_filtered[!grepl("^ENSG", scpredixcan_filtered$gene_name), ]
celltype_mapping <- list(
  "DC" = c("conventional_dendritic cell", "dendritic_cell", "plasmacytoid_dendritic_cell"),
  "Mono" = c("CD14-low_CD16-positive_monocyte", "CD14-positive_monocyte"),
  "NK" = c("CD16-negative_CD56-bright_natural_killer_cell_human", 
           "natural_killer_cell"),
  "CD4T" = c("CD4-positive_alpha-beta_cytotoxic_T_cell",
             "CD4-positive_alpha-beta_T_cell",
             "central_memory_CD4-positive_alpha-beta_T_cell",
             "effector_memory_CD4-positive_alpha-beta_T_cell",
             "naive_thymus-derived_CD4-positive_alpha-beta_T_cell",
             "regulatory_T_cell"),
  "CD8T" = c("CD8-positive_alpha-beta_T_cell",
             "central_memory_CD8-positive_alpha-beta_T_cell",
             "effector_memory_CD8-positive_alpha-beta_T_cell",
             "naive_thymus-derived_CD8-positive_alpha-beta_T_cell"),
  "B" = c("memory_B_cell", "naive_B_cell", "transitional_stage_B_cell")
)
celltype_mapping_df <- data.frame(
  detailed_celltype = unlist(celltype_mapping),
  main_celltype = rep(names(celltype_mapping), 
                      sapply(celltype_mapping, length)),
  stringsAsFactors = FALSE
)
scpredixcan_filtered$cell_type <- celltype_mapping_df$main_celltype[
  match(scpredixcan_filtered$cell_type_21, celltype_mapping_df$detailed_celltype)
]
#scpredixcan_filtered <- scpredixcan_filtered[scpredixcan_filtered$cell_type %in% celltypes, ]

set_pred <- unique(twas_pred_filtered$gene_name)
set_ct <- unique(TWAS_ct_filtered$gene_name) 
set_scpred <- unique(scpredixcan_filtered$gene_name)
all_genes <- unique(c(set_pred, set_ct, set_scpred))
comb_matrix <- data.frame(
  gene = all_genes,
  in_pred = as.integer(all_genes %in% set_pred),
  in_ct = as.integer(all_genes %in% set_ct),
  in_scpred = as.integer(all_genes %in% set_scpred)
)
comb_matrix$combination <- paste0(
  comb_matrix$in_pred,
  comb_matrix$in_ct,
  comb_matrix$in_scpred
)
comb_counts <- table(comb_matrix$combination)
max_length <- max(comb_counts) + 1
dt <- as.data.table(comb_matrix)

long_data <- comb_matrix %>%
  group_by(combination) %>%
  dplyr::mutate(row_num = row_number() + 1) %>%
  ungroup() %>%
  select(combination, gene, row_num)
count_rows <- data.frame(
  combination = names(comb_counts),
  gene = as.character(comb_counts),
  row_num = 1
)
all_data <- bind_rows(count_rows, long_data)
wide_data <- all_data %>%
  pivot_wider(
    id_cols = row_num,
    names_from = combination,
    values_from = gene,
    values_fill = NA
  ) %>%
  arrange(row_num) %>%
  select(-row_num)
wide_data <- cbind(
  row_type = c("count", rep("gene", nrow(wide_data)-1)),
  wide_data
)
colnames(wide_data) <- c("Method","scPrediXcan","gdTruthct","gdTruthct_scPrediXcan",
                         "MediCell","MediCell_scPrediXcan","MediCell_gdTruth","all")
write.csv(wide_data, "gene_comparisonmedicell_gdtruth_scpredixcan.csv", row.names = FALSE)

gene_lists <- list(
  "MediCell" = unique(twas_pred_filtered$gene_name),
  "gdTruth_ct" = unique(TWAS_ct_filtered$gene_name),
  "scPredXcan" = unique(scpredixcan_filtered$gene_name)
)
pdf("upset_gdTruth_MediCell_scPredXcan.pdf", width = 8, height = 6)
upset(fromList(gene_lists), 
      order.by = "freq",
      mb.ratio = c(0.6, 0.4),
      text.scale = c(2.0, 2.0, 1.5, 1.5, 2.0, 2.0))
dev.off()

#heatmap
common <- na.omit(wide_data$all[-1])
gdTruthct_scPrediXcan <- na.omit(wide_data$gdTruthct_scPrediXcan[-1])
MediCell_scPrediXcan <- na.omit(wide_data$MediCell_scPrediXcan[-1])
MediCell_gdTruth <- na.omit(wide_data$MediCell_gdTruth[-1])
target_genes <- unique(c(common,MediCell_gdTruth,MediCell_scPrediXcan,gdTruthct_scPrediXcan))
target_genes <- factor(target_genes, levels = target_genes)
cell_order <- c("CD4T", "CD8T", "NK", "B", "Mono", "DC")

get_sig_matrix <- function(dt, method_name){
  dt <- as.data.table(dt)
  dt[, bonferroni_p := pvalue * ngeneindb]
  dt[, sig := as.integer(bonferroni_p < 0.05)]
  dt <- dt[gene_name %in% target_genes & cell_type %in% cell_order]
  dt <- dt[, .(gene_name, cell_type, sig)]
  setnames(dt, "sig", paste0("sig_", method_name))
  return(dt)
}
sig_MediCell <- get_sig_matrix(twas_pred_filtered, "MediCell")
sig_gdTruthct <- get_sig_matrix(TWAS_ct_filtered, "gdTruthct")
sig_scPred <- get_sig_matrix(scpredixcan_filtered, "scPrediXcan")

comb <- merge(sig_MediCell, sig_gdTruthct, by = c("gene_name","cell_type"), all = TRUE)
comb <- merge(comb, sig_scPred, by = c("gene_name","cell_type"), all = TRUE)
comb[is.na(comb)] <- 0
comb[, total_sig := sig_MediCell + sig_gdTruthct + sig_scPrediXcan]
full_grid <- expand.grid(gene_name = target_genes, cell_type = cell_order, stringsAsFactors = FALSE)
setDT(full_grid)
plot_data <- merge(full_grid, comb[, .(gene_name, cell_type, total_sig)], 
                   by = c("gene_name","cell_type"), all.x = TRUE)
plot_data[is.na(total_sig), total_sig := 0]
plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_order)
plot_data$gene_name <- factor(plot_data$gene_name, levels = target_genes)
write.csv(comb, file = "cell_sig_heatmap.csv")


#heatmap2
cell_order <- c("CD4T", "CD8T", "NK", "B", "Mono", "DC")
method_colors <- c("scPrediXcan" = "#a774c9", "MediCell" = "#f67d25", "Pseudobulk" = "#89c760")
get_sig_matrix <- function(dt, method_name){
  dt <- as.data.table(dt)
  dt[, bonferroni_p := pvalue * ngeneindb]
  dt[, sig := as.integer(bonferroni_p < 0.05)]
  dt <- dt[gene_name %in% target_genes & cell_type %in% cell_order]
  dt <- dt[, .(gene_name, cell_type, sig)]
  setnames(dt, "sig", paste0("sig_", method_name))
  return(dt)
}
sig_MediCell <- get_sig_matrix(twas_pred_filtered, "MediCell")
sig_gdTruthct <- get_sig_matrix(TWAS_ct_filtered, "gdTruthct")
sig_scPred <- get_sig_matrix(scpredixcan_filtered, "scPrediXcan")

comb <- merge(sig_MediCell, sig_gdTruthct, by = c("gene_name","cell_type"), all = TRUE)
comb <- merge(comb, sig_scPred, by = c("gene_name","cell_type"), all = TRUE)
comb[is.na(comb)] <- 0

write.csv(comb, file = "cell_sig_heatmap.csv", row.names = F)

dt_long <- melt(
  comb,
  id.vars = c("gene_name", "cell_type"),
  measure.vars = c("sig_MediCell", "sig_gdTruthct", "sig_scPrediXcan"),
  variable.name = "method_raw",
  value.name = "is_sig"
)
dt_long[, method := fcase(
  method_raw == "sig_MediCell", "MediCell",
  method_raw == "sig_gdTruthct", "Pseudobulk",
  method_raw == "sig_scPrediXcan", "scPrediXcan"
)]
method_order <- c("Pseudobulk", "MediCell", "scPrediXcan")
dt_long$method <- factor(dt_long$method, levels = method_order)
dt_long[, cell_method := paste0(cell_type, "_", method)]
full_grid_long <- expand.grid(
  gene_name = target_genes,
  cell_type = cell_order,
  method = method_order,
  stringsAsFactors = FALSE
)
setDT(full_grid_long)
full_grid_long[, cell_method := paste0(cell_type, "_", method)]
plot_data <- merge(
  full_grid_long[, .(gene_name, cell_type, method, cell_method)],
  dt_long[, .(gene_name, cell_method, is_sig)],
  by = c("gene_name", "cell_method"),
  all.x = TRUE
)
plot_data[is.na(is_sig), is_sig := 0]
level_vec <- c()
for(ct in cell_order){
  for(m in method_order){
    level_vec <- c(level_vec, paste0(ct, "_", m))
  }
}
plot_data[, cell_method := factor(cell_method, levels = level_vec)]
plot_data$gene_name <- factor(plot_data$gene_name, levels = rev(target_genes))
plot_data$ct_index <- match(plot_data$cell_type, cell_order)

pdf("gene_cell_method_separate_heatmap.pdf", width = 10, height = 6)
ggplot(plot_data, aes(x = cell_method, y = gene_name)) +
  geom_tile(aes(fill = ifelse(is_sig == 1, method, NA)), color = "gray90", size = 0.3) +
  geom_vline(xintercept = seq(3.5, 15.5, by=3), linewidth=0.7, color="black") +
  scale_x_discrete(expand = c(-0.001, -0.001)) +
  scale_y_discrete(expand = c(-0.001, -0.001)) +
  scale_fill_manual(values = method_colors, na.value = "white", name = "Method", drop = FALSE) +
  labs(x = "Cell type & Method", y = "Gene") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 14),
    legend.position = "bottom",
    panel.grid = element_blank()
  )
dev.off()

#zscore棒棒糖图
get_z_matrix <- function(dt, method_name){
  dt <- as.data.table(dt)
  dt <- dt[gene_name %in% target_genes & cell_type %in% cell_order]
  dt <- dt[, .(gene_name, cell_type, zscore)]
  dt <- dt[, .(zscore = mean(zscore)), by = .(gene_name, cell_type)]
  setnames(dt, "zscore", paste0("z_", method_name))
  return(dt)
}
z_MediCell   <- get_z_matrix(twas_pred_filtered, "MediCell")
z_gdTruthct  <- get_z_matrix(TWAS_ct_filtered, "gdTruthct")
z_scPred     <- get_z_matrix(scpredixcan_filtered, "scPrediXcan")
comb_zscore <- expand.grid(
  gene_name = target_genes,
  cell_type = cell_order,
  stringsAsFactors = FALSE
)
setDT(comb_zscore)
comb_zscore <- merge(comb_zscore, z_MediCell,  by=c("gene_name","cell_type"), all.x=T)
comb_zscore <- merge(comb_zscore, z_gdTruthct, by=c("gene_name","cell_type"), all.x=T)
comb_zscore <- merge(comb_zscore, z_scPred,    by=c("gene_name","cell_type"), all.x=T)
plot_data_zscore <- melt(
  comb_zscore,
  id.vars = c("gene_name", "cell_type"),
  measure.vars = c("z_MediCell", "z_gdTruthct", "z_scPrediXcan"),
  variable.name = "Method",
  value.name = "zscore"
)
plot_data_zscore <- plot_data_zscore[!is.na(zscore)]
plot_data_zscore[, Method := gsub("z_", "", Method)]

shape_map <- c("MediCell" = 16, "gdTruthct" = 15, "scPrediXcan" = 17)
method_colors <- c("MediCell" = "#f67d25", "gdTruthct" = "#89c760", "scPrediXcan" = "#a774c9")
plot_data_zscore$cell_type <- factor(plot_data_zscore$cell_type, levels = cell_order)
plot_data_zscore$gene_name <- factor(plot_data_zscore$gene_name, levels = rev(target_genes))
plot_data_zscore$Method <- factor(plot_data_zscore$Method)
write.csv(comb_zscore, file = "cell_sig_zscore_lollipop.csv", row.names = F)

pdf("gene_cell_zscore_lollipop_full.pdf", width=20, height=7)
ggplot(plot_data_zscore, aes(x = zscore, y = gene_name)) +
  geom_vline(xintercept = 0, color="black", linewidth=0.8) +
  geom_linerange(aes(xmin=0, xmax=zscore), linewidth=0.6, alpha=0.6) +
  geom_point(aes(shape = Method, color = Method), size = 4.5, stroke=1, alpha=0.95) +
  facet_wrap(~cell_type, nrow = 1, ncol = 6, scales = "free_y", strip.position = "top") +
  scale_color_manual(values = method_colors, name="Method") +
  scale_shape_manual(values = shape_map, name="Method") +
  labs(x = "Z-score", y = "Gene") +
  theme_bw() +
  theme(
    strip.text = element_text(size=16),
    strip.background = element_rect(fill="#f7f7f7", color="black"),
    axis.text.y = element_text(size=16),
    axis.text.x = element_text(size=16),
    axis.title = element_text(size=16),
    legend.position = "bottom",
    legend.box = "vertical",
    panel.spacing = unit(0.3, "lines")
  )
dev.off()
