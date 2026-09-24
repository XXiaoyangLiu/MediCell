library(stringr)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(data.table)
library(tidyr)

gene_pos_df <- fread("/node1/liuxy/index/gene_pos_grch37.82.csv", data.table = F)
gene_pos_df$start <- as.numeric(gene_pos_df$start)
gene_pos_df$end   <- as.numeric(gene_pos_df$end) 
mhc_genes <- gene_pos_df %>%
  filter(chr == "6", start >= 28477797, end <= 33448354) %>%
  pull(gene_name)

trait_files <- c(
  "/node1/liuxy/MediCell/kidney_3/glom/TWAS/output/glom/SLE/summary/Bonfsig_allcelltype_result.csv",
  "/node1/liuxy/MediCell/trait_compare_2/eGFR/output/TWAS/glom/eGFR/summary/Bonfsig_allcelltype_result.csv",
  "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/UA/glom/summary/Bonfsig_allcelltype_result.csv",
  "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uAlb/glom/summary/Bonfsig_allcelltype_result.csv",
  "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uCr/glom/summary/Bonfsig_allcelltype_result.csv",
  "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uK/glom/summary/Bonfsig_allcelltype_result.csv",
  "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uNa/glom/summary/Bonfsig_allcelltype_result.csv",
  "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/SBP/glom/summary/Bonfsig_allcelltype_result.csv"
)
dt_list <- lapply(trait_files, function(f) read.csv(f))
all_twas <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
all_twas_nomhc <- all_twas[!all_twas$gene_name %in% mhc_genes,]
gene_trait_count <- all_twas_nomhc[, .(n_trait = uniqueN(gwas_id)), by = gene_name]
unique_gene <- gene_trait_count[n_trait == 1, gene_name]
twas_all_specific <- all_twas_nomhc[gene_name %in% unique_gene]
gwas_SLE <- fread("/node1/liuxy/check/other_SLE_GWAS/metal/20260830/SLE_metal_096_917_183_filt_3.tsv")
gwas_SLE$Allele1 <- toupper(gwas_SLE$Allele1)
gwas_SLE$Allele2 <- toupper(gwas_SLE$Allele2)
gwas_SLE[, chr := as.integer(gsub("^(\\d+)_.*", "\\1", MarkerName))]
gwas_SLE[, pos := as.integer(gsub("^\\d+_(\\d+)_.*", "\\1", MarkerName))]
gwas_SLE$logP <- -log10(gwas_SLE$`P-value`)
prep_manhattan_coord <- function(dat){
  dat %>%
    arrange(chr, pos) %>%
    group_by(chr) %>%
    mutate(
      chr_min = min(pos),
      chr_max = max(pos),
      chr_span = chr_max - chr_min
    ) %>%
    ungroup() %>%
    group_by(chr) %>%
    mutate(rel_pos = (pos - chr_min) / (chr_span + 1)) %>%
    ungroup() %>%
    mutate(
      max_allowed_width = 1e5,
      actual_width = ifelse(chr_span > max_allowed_width, max_allowed_width, chr_span),
      compressed_rel_pos = ifelse(chr_span > max_allowed_width, rel_pos * (max_allowed_width / chr_span), rel_pos),
      chr_start = cumsum(lag(actual_width, default = 0) + 1e5),
      cumulative_pos = chr_start + compressed_rel_pos * actual_width
    )
}
gwas_coord <- prep_manhattan_coord(gwas_SLE)
gwas_coord$group <- "background"
axis_df <- gwas_coord %>%
  group_by(chr) %>%
  summarise(center = mean(cumulative_pos, na.rm = TRUE), .groups = "drop")

plot_color <- c("background" = "gray70", 
                "SLE" = "#d62728",
                "eGFR" = "#f67d25",
                "UA" = "#9983bd",
                "uAlb" = "#5b2d7c",
                "uNa" = "#bf6cac",
                "uK" = "#89c760",
                "uCr" = "#27ae60",
                "SBP" = "#fee301")

priority_list <- c("SLE", "eGFR", "UA", "uAlb", "uCr", "uK", "uNa", "SBP")
out_dir <- "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/plot/"
cell_type_list <- c("Endo", "Podo", "MC", "MNP", "NK", "NKT", "CD8T", "Th_cell", "Bcell")
for(ct in cell_type_list){
  weight_dir <- "/node1/liuxy/MediCell/kidney_3/glom/predictdb/output/weights/"
  file_pattern <- paste0(weight_dir, "Model_training_", ct, "_chr*_weights.txt")
  file_list <- Sys.glob(file_pattern)
  chr_num <- as.integer(gsub(".*chr(\\d+).*", "\\1", file_list))
  file_list <- file_list[order(chr_num)]
  dt_list <- lapply(file_list, function(f) fread(f, sep = "\t", header = TRUE))
  cell_weights <- rbindlist(dt_list)
  setnames(cell_weights, "gene_name", "gene_id")
  twas_ct <- twas_all_specific[twas_all_specific$cell_type %in% ct,]
  twas_snp_link <- merge(
    x = twas_ct,
    y = cell_weights,
    by.x = "gene_id",
    by.y = "gene_id",
    all = FALSE
  )
  twas_snp_link[, position := as.numeric(str_extract(varID, "_(\\d+)_", group = 1))]
  twas_snp_link$snpID <- twas_snp_link$varID
  setnames(twas_snp_link, "gwas_id", "group")
  twas_sub <- select(twas_snp_link, snpID, position, group, gene_name)
  gwas_sub <- select(gwas_coord, snpID = MarkerName, position=pos, cumulative_pos, logP)
  twas_snp_gwas <- inner_join(twas_sub, gwas_sub, by = c("snpID", "position"))
  gene_label <- twas_snp_gwas %>%
    group_by(gene_name, group) %>%
    summarise(
      label_x = mean(cumulative_pos),
      max_y = max(logP) + 0.25,
      label_color = unique(group),
      .groups = "drop"
    )
  gene_trait_count <- gene_label %>% count(gene_name, name = "trait_num")
  keep_single_trait_gene <- gene_trait_count %>% filter(trait_num == 1) %>% pull(gene_name)
  gene_label_single <- gene_label %>% filter(gene_name %in% keep_single_trait_gene)
  twas_snp_gwas_single <- twas_snp_gwas %>% filter(gene_name %in% keep_single_trait_gene)
  pdf_name <- paste0(out_dir, ct, "_TWAS_SNP_manhattan_SLE_bg_glom.pdf")
  pdf(pdf_name, width = 40, height = 25)
  p <- ggplot() +
    geom_point(
      data = gwas_coord,
      aes(x = cumulative_pos, y = logP),
      color = "gray70", size = 0.2
    ) +
    geom_point(
      data = twas_snp_gwas_single,
      aes(x = cumulative_pos, y = logP, color = group),
      size = 0.7
    ) +
    geom_text_repel(
      data = gene_label_single,
      aes(x = label_x, y = max_y, label = gene_name, color = label_color),
      size = 2.6,
      point.padding = 0.15,
      max.overlaps = 25, 
      min.segment.length = 0,
      segment.size = 0.2
    ) +
    scale_x_continuous(
      name = "Chromosome",
      breaks = axis_df$center, 
      labels = axis_df$chr,
      expand = expansion(mult = 0.01)
    ) +
    scale_y_continuous(
      name = expression(-log[10](italic(p))),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    scale_color_manual(values = plot_color, limits = names(plot_color)) +
    labs(title = paste0("Cell type: ", ct)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = element_line(color = "black", linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )
  print(p)
  dev.off()
  cat("已保存：", pdf_name, "\n\n")
}

for(ct in cell_type_list){
  weight_dir <- "/node1/liuxy/MediCell/kidney_3/glom/predictdb/output/weights/"
  file_pattern <- paste0(weight_dir, "Model_training_", ct, "_chr*_weights.txt")
  file_list <- Sys.glob(file_pattern)
  chr_num <- as.integer(gsub(".*chr(\\d+).*", "\\1", file_list))
  file_list <- file_list[order(chr_num)]
  dt_list <- lapply(file_list, function(f) fread(f, sep = "\t", header = TRUE))
  cell_weights <- rbindlist(dt_list)
  setnames(cell_weights, "gene_name", "gene_id")
  twas_ct <- twas_all_specific[twas_all_specific$cell_type %in% ct,]
  twas_snp_link <- merge(
    x = twas_ct,
    y = cell_weights,
    by.x = "gene_id",
    by.y = "gene_id",
    all = FALSE
  )
  twas_snp_link[, position := as.numeric(str_extract(varID, "_(\\d+)_", group = 1))]
  twas_snp_link$snpID <- twas_snp_link$varID
  setnames(twas_snp_link, "gwas_id", "group")
  twas_sub <- select(twas_snp_link, snpID, position, group, gene_name)
  gwas_sub <- select(gwas_coord, snpID = MarkerName, position=pos, cumulative_pos, logP)
  twas_snp_gwas <- inner_join(twas_sub, gwas_sub, by = c("snpID", "position"))
  gene_label <- twas_snp_gwas %>%
    group_by(gene_name, group) %>%
    summarise(
      label_x = mean(cumulative_pos),
      max_y = max(logP) + 0.25,
      label_color = unique(group),
      .groups = "drop"
    )
  gene_trait_count <- gene_label %>% count(gene_name, name = "trait_num")
  keep_single_trait_gene <- gene_trait_count %>% filter(trait_num == 1) %>% pull(gene_name)
  gene_label_single <- gene_label %>% filter(gene_name %in% keep_single_trait_gene)
  twas_snp_gwas_single <- twas_snp_gwas %>% filter(gene_name %in% keep_single_trait_gene)
  pdf_name <- paste0(out_dir, ct, "_TWAS_SNP_manhattan_SLE_bg_glom_nolabel_2.pdf")
  pdf(pdf_name, width = 35, height = 15)
  p <- ggplot() +
    geom_point(
      data = gwas_coord,
      aes(x = cumulative_pos, y = logP),
      color = "gray70", size = 0.2
    ) +
    geom_point(
      data = twas_snp_gwas_single,
      aes(x = cumulative_pos, y = logP, color = group),
      size = 0.7
    ) +
    geom_point(
    data = axis_df,
    aes(x = center, y = 0),
    color = "black", size = 1.3, shape = 20
  ) +
    scale_x_continuous(
      name = "Chromosome",
      breaks = axis_df$center,
      labels = axis_df$chr,
      expand = expansion(mult = 0.01)
    ) +
    scale_y_continuous(
      name = expression(-log[10](italic(p))),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    scale_color_manual(values = plot_color, limits = names(plot_color)) +
    labs(title = paste0("Cell type: ", ct)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = element_line(color = "black", linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )
  print(p)
  dev.off()
}

#tub
library(stringr)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(data.table)
library(tidyr)

gene_pos_df <- fread("/node1/liuxy/index/gene_pos_grch37.82.csv", data.table = F)
gene_pos_df$start <- as.numeric(gene_pos_df$start)
gene_pos_df$end   <- as.numeric(gene_pos_df$end) 
mhc_genes <- gene_pos_df %>%
  filter(chr == "6", start >= 28477797, end <= 33448354) %>%
  pull(gene_name)

trait_files <- c(
  "/node1/liuxy/MediCell/kidney_3/tub/TWAS/output/tub/SLE/summary/Bonfsig_allcelltype_result.csv",
  "/node1/liuxy/MediCell/trait_compare_2/eGFR/output/TWAS/tub/eGFR/summary/Bonfsig_allcelltype_result.csv",
  "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/UA/tub/summary/Bonfsig_allcelltype_result.csv",
  "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uAlb/tub/summary/Bonfsig_allcelltype_result.csv",
  "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uCr/tub/summary/Bonfsig_allcelltype_result.csv",
  "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uK/tub/summary/Bonfsig_allcelltype_result.csv",
  "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uNa/tub/summary/Bonfsig_allcelltype_result.csv",
  "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/SBP/tub/summary/Bonfsig_allcelltype_result.csv"
)
dt_list <- lapply(trait_files, function(f) read.csv(f))
all_twas <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
all_twas_nomhc <- all_twas[!all_twas$gene_name %in% mhc_genes,]
gene_trait_count <- all_twas_nomhc[, .(n_trait = uniqueN(gwas_id)), by = gene_name]
unique_gene <- gene_trait_count[n_trait == 1, gene_name]
twas_all_specific <- all_twas_nomhc[gene_name %in% unique_gene]
gwas_SLE <- fread("/node1/liuxy/check/other_SLE_GWAS/metal/20260830/SLE_metal_096_917_183_filt_3.tsv")
gwas_SLE$Allele1 <- toupper(gwas_SLE$Allele1)
gwas_SLE$Allele2 <- toupper(gwas_SLE$Allele2)
gwas_SLE[, chr := as.integer(gsub("^(\\d+)_.*", "\\1", MarkerName))]
gwas_SLE[, pos := as.integer(gsub("^\\d+_(\\d+)_.*", "\\1", MarkerName))]
gwas_SLE$logP <- -log10(gwas_SLE$`P-value`)
prep_manhattan_coord <- function(dat){
  dat %>%
    arrange(chr, pos) %>%
    group_by(chr) %>%
    mutate(
      chr_min = min(pos),
      chr_max = max(pos),
      chr_span = chr_max - chr_min
    ) %>%
    ungroup() %>%
    group_by(chr) %>%
    mutate(rel_pos = (pos - chr_min) / (chr_span + 1)) %>%
    ungroup() %>%
    mutate(
      max_allowed_width = 1e5,
      actual_width = ifelse(chr_span > max_allowed_width, max_allowed_width, chr_span),
      compressed_rel_pos = ifelse(chr_span > max_allowed_width, rel_pos * (max_allowed_width / chr_span), rel_pos),
      chr_start = cumsum(lag(actual_width, default = 0) + 1e5),
      cumulative_pos = chr_start + compressed_rel_pos * actual_width
    )
}
gwas_coord <- prep_manhattan_coord(gwas_SLE)
gwas_coord$group <- "background"
axis_df <- gwas_coord %>%
  group_by(chr) %>%
  summarise(center = mean(cumulative_pos, na.rm = TRUE), .groups = "drop")

plot_color <- c("background" = "gray70", 
                "SLE" = "#d62728",
                "eGFR" = "#f67d25",
                "UA" = "#9983bd",
                "uAlb" = "#5b2d7c",
                "uNa" = "#bf6cac",
                "uK" = "#89c760",
                "uCr" = "#27ae60",
                "SBP" = "#fee301")

priority_list <- c("SLE", "eGFR", "UA", "uAlb", "uCr", "uK", "uNa", "SBP")
out_dir <- "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/plot/tub_manhattan"
cell_type_list <- c("PT","LOH","DCT","CDICA","CDICB","MNP","Plasmacytoid_DC","NK","NKT","CD8T","Bcell","Th_cell")
for(ct in cell_type_list){
  weight_dir <- "/node1/liuxy/MediCell/kidney_3/tub/predictdb/output/weights/"
  file_pattern <- paste0(weight_dir, "Model_training_", ct, "_chr*_weights.txt")
  file_list <- Sys.glob(file_pattern)
  chr_num <- as.integer(gsub(".*chr(\\d+).*", "\\1", file_list))
  file_list <- file_list[order(chr_num)]
  dt_list <- lapply(file_list, function(f) fread(f, sep = "\t", header = TRUE))
  cell_weights <- rbindlist(dt_list)
  setnames(cell_weights, "gene_name", "gene_id")
  twas_ct <- twas_all_specific[twas_all_specific$cell_type %in% ct,]
  twas_snp_link <- merge(
    x = twas_ct,
    y = cell_weights,
    by.x = "gene_id",
    by.y = "gene_id",
    all = FALSE
  )
  twas_snp_link[, position := as.numeric(str_extract(varID, "_(\\d+)_", group = 1))]
  twas_snp_link$snpID <- twas_snp_link$varID
  setnames(twas_snp_link, "gwas_id", "group")
  twas_sub <- select(twas_snp_link, snpID, position, group, gene_name)
  gwas_sub <- select(gwas_coord, snpID = MarkerName, position=pos, cumulative_pos, logP)
  twas_snp_gwas <- inner_join(twas_sub, gwas_sub, by = c("snpID", "position"))
  gene_label <- twas_snp_gwas %>%
    group_by(gene_name, group) %>%
    summarise(
      label_x = mean(cumulative_pos),
      max_y = max(logP) + 0.25,
      label_color = unique(group),
      .groups = "drop"
    )
  gene_trait_count <- gene_label %>% count(gene_name, name = "trait_num")
  keep_single_trait_gene <- gene_trait_count %>% filter(trait_num == 1) %>% pull(gene_name)
  gene_label_single <- gene_label %>% filter(gene_name %in% keep_single_trait_gene)
  twas_snp_gwas_single <- twas_snp_gwas %>% filter(gene_name %in% keep_single_trait_gene)
  pdf_name <- paste0(out_dir, ct, "_TWAS_SNP_manhattan_SLE_bg_glom.pdf")
  pdf(pdf_name, width = 40, height = 25)
  p <- ggplot() +
    geom_point(
      data = gwas_coord,
      aes(x = cumulative_pos, y = logP),
      color = "gray70", size = 0.2
    ) +
    geom_point(
      data = twas_snp_gwas_single,
      aes(x = cumulative_pos, y = logP, color = group),
      size = 0.7
    ) +
    geom_text_repel(
      data = gene_label_single,
      aes(x = label_x, y = max_y, label = gene_name, color = label_color),
      size = 2.6,
      point.padding = 0.15,
      max.overlaps = 25, 
      min.segment.length = 0,
      segment.size = 0.2
    ) +
    scale_x_continuous(
      name = "Chromosome",
      breaks = axis_df$center, 
      labels = axis_df$chr,
      expand = expansion(mult = 0.01)
    ) +
    scale_y_continuous(
      name = expression(-log[10](italic(p))),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    scale_color_manual(values = plot_color, limits = names(plot_color)) +
    labs(title = paste0("Cell type: ", ct)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = element_line(color = "black", linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )
  print(p)
  dev.off()
  cat("已保存：", pdf_name, "\n\n")
}

for(ct in cell_type_list){
  weight_dir <- "/node1/liuxy/MediCell/kidney_3/tub/predictdb/output/weights/"
  file_pattern <- paste0(weight_dir, "Model_training_", ct, "_chr*_weights.txt")
  file_list <- Sys.glob(file_pattern)
  chr_num <- as.integer(gsub(".*chr(\\d+).*", "\\1", file_list))
  file_list <- file_list[order(chr_num)]
  dt_list <- lapply(file_list, function(f) fread(f, sep = "\t", header = TRUE))
  cell_weights <- rbindlist(dt_list)
  setnames(cell_weights, "gene_name", "gene_id")
  twas_ct <- twas_all_specific[twas_all_specific$cell_type %in% ct,]
  twas_snp_link <- merge(
    x = twas_ct,
    y = cell_weights,
    by.x = "gene_id",
    by.y = "gene_id",
    all = FALSE
  )
  twas_snp_link[, position := as.numeric(str_extract(varID, "_(\\d+)_", group = 1))]
  twas_snp_link$snpID <- twas_snp_link$varID
  setnames(twas_snp_link, "gwas_id", "group")
  twas_sub <- select(twas_snp_link, snpID, position, group, gene_name)
  gwas_sub <- select(gwas_coord, snpID = MarkerName, position=pos, cumulative_pos, logP)
  twas_snp_gwas <- inner_join(twas_sub, gwas_sub, by = c("snpID", "position"))
  gene_label <- twas_snp_gwas %>%
    group_by(gene_name, group) %>%
    summarise(
      label_x = mean(cumulative_pos),
      max_y = max(logP) + 0.25,
      label_color = unique(group),
      .groups = "drop"
    )
  gene_trait_count <- gene_label %>% count(gene_name, name = "trait_num")
  keep_single_trait_gene <- gene_trait_count %>% filter(trait_num == 1) %>% pull(gene_name)
  gene_label_single <- gene_label %>% filter(gene_name %in% keep_single_trait_gene)
  twas_snp_gwas_single <- twas_snp_gwas %>% filter(gene_name %in% keep_single_trait_gene)
  pdf_name <- paste0(out_dir, ct, "_TWAS_SNP_manhattan_SLE_bg_glom_nolabel_2.pdf")
  pdf(pdf_name, width = 35, height = 15)
  p <- ggplot() +
    geom_point(
      data = gwas_coord,
      aes(x = cumulative_pos, y = logP),
      color = "gray70", size = 0.2
    ) +
    geom_point(
      data = twas_snp_gwas_single,
      aes(x = cumulative_pos, y = logP, color = group),
      size = 0.7
    ) +
    geom_point(
    data = axis_df,
    aes(x = center, y = 0),
    color = "black", size = 1.3, shape = 20
  ) +
    scale_x_continuous(
      name = "Chromosome",
      breaks = axis_df$center,
      labels = axis_df$chr,
      expand = expansion(mult = 0.01)
    ) +
    scale_y_continuous(
      name = expression(-log[10](italic(p))),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    scale_color_manual(values = plot_color, limits = names(plot_color)) +
    labs(title = paste0("Cell type: ", ct)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = element_line(color = "black", linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )
  print(p)
  dev.off()
}
