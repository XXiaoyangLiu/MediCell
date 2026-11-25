#manhattan
library(stringr)
library(dplyr)
library(ggrepel)
library(ggplot2)
gene_pos_df <- read.csv("/node1/liuxy/index/gene_pos_grch37.82.csv")
mhc_genes <- gene_pos_df$gene_name[gene_pos_df$chr %in% "6" & as.numeric(gene_pos_df$start) >= 28477797 & as.numeric(gene_pos_df$end) <= 33448354]
#allcelltype
base_paths <- "/node1/liuxy/MediCell/pbmc_2/TWAS/output_check/pbmc/SLE/"
celltypes <- c("CD4T", "CD8T", "NK", "B", "Mono", "DC")
out_dir <- "/node1/liuxy/MediCell/pbmc_2/TWAS/output/plot/manhattan/SLE_char_check/"

for (ct in celltypes) {
    file_pattern <- paste0(ct, "\\.csv$")
    files <- list.files(base_paths, pattern = file_pattern, full.names = TRUE)
    if (length(files) == 0) next
    for (f in files) {
      message("Processing: ", f)
      df <- tryCatch(read.csv(f), error = function(e) { message("read failed: ", f); return(NULL) })
      if (is.null(df)) next
      df$bonfp <- p.adjust(df$pvalue, method = "bonferroni")
      df <- merge(df, gene_pos_df, by.x = "gene", by.y = "gene_id")
      plot_df <- df %>%
        filter(chr %in% c(1:22)) %>%
        mutate(position = as.numeric(start)) %>%
        mutate(chr = factor(chr, levels = c(1:22))) %>%
        arrange(chr, position) %>%
        group_by(chr) %>%
        mutate(chr_min = min(position), chr_max = max(position), chr_span = chr_max - chr_min ) %>%
        ungroup() %>%
        group_by(chr) %>%
        mutate(rel_pos = (position - chr_min) / (chr_span + 1)) %>%
        ungroup() %>%
        mutate(max_allowed_width = 1e5,
               actual_width = ifelse(chr_span > max_allowed_width, max_allowed_width, chr_span),
               compressed_rel_pos = ifelse(chr_span > max_allowed_width, rel_pos * (max_allowed_width / chr_span), rel_pos),
               chr_start = cumsum(lag(actual_width, default = 0) + 1e5),
               cumulative_pos = chr_start + compressed_rel_pos * actual_width)
      axis_df <- plot_df %>%
        group_by(chr) %>%
        summarise(center = mean(cumulative_pos, na.rm = TRUE), span = max(position) - min(position), compressed_span = max(cumulative_pos) - min(cumulative_pos))
      bonfp_threshold <- max(plot_df$pvalue[plot_df$bonfp < 0.05], na.rm = TRUE)
      genes_to_label <- plot_df %>%
        filter(bonfp < 0.05) %>%
        filter(!gene_name.y %in% mhc_genes) %>%
        pull(gene_name.y)
      label_data <- plot_df %>%
        filter(gene_name.y %in% genes_to_label)

      pdf(paste0("/node1/liuxy/MediCell/pbmc_2/TWAS/output/plot/manhattan/SLE_char_check/SLE_manhattan_pbmc_", ct, "_char.pdf"), width = 20, height = 10)
      p <- ggplot(plot_df, aes(x = cumulative_pos, y = -log10(pvalue))) +
        geom_point(aes(color = chr), size = 4, shape = 16) +
        labs(title = paste("TWAS Manhattan Plot in PBMC", ct)) +
        scale_x_continuous(name = "Chromosome", breaks = axis_df$center, labels = axis_df$chr, expand = expansion(mult = 0.01)) +
        geom_text_repel(data = label_data, aes(label = gene_name.y), size = 3, point.padding = 0.3) +
        scale_y_continuous(name = expression(-log[10](italic(p))), expand = expansion(mult = c(0.02, 0.05))) + 
        scale_color_manual(values = rep(c("#894297", "#FFC107"), length.out = length(unique(plot_df$chr)))) +
        theme_minimal() +
        theme(legend.position = "none",
              axis.line.y = element_line(color = "black", linewidth = 0.5),
              axis.ticks.y = element_line(color = "black", linewidth = 0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
            #  axis.text.x = element_text(size = 20, angle = 90),
            #  axis.text.y = element_text(size = 20),
            #  axis.title = element_text(size = 20),
            #  plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
              plot.background = element_rect(fill = "white", color = NA))
        print(p)
      dev.off()

      pdf(paste0("/node1/liuxy/MediCell/pbmc_2/TWAS/output/plot/manhattan/SLE_char_check/SLE_manhattan_pbmc_", ct, "_nochar.pdf"), width = 20, height = 10)
      p <- ggplot(plot_df, aes(x = cumulative_pos, y = -log10(pvalue))) +
        geom_point(aes(color = chr), size = 4, shape = 16) +
        labs(title = paste("TWAS Manhattan Plot in PBMC", ct)) +
        scale_x_continuous(name = "Chromosome", breaks = axis_df$center, labels = axis_df$chr, expand = expansion(mult = 0.01)) +
        #geom_text_repel(data = label_data, aes(label = gene_name.y), size = 3, point.padding = 0.3) +
        scale_y_continuous(name = expression(-log[10](italic(p))), expand = expansion(mult = c(0.02, 0.05))) +
        scale_color_manual(values = rep(c("#894297", "#FFC107"), length.out = length(unique(plot_df$chr)))) +
        theme_minimal() +
        theme(legend.position = "none",
              axis.line.y = element_line(color = "black", linewidth = 0.5),
              axis.ticks.y = element_line(color = "black", linewidth = 0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
            #  axis.text.x = element_text(size = 20, angle = 90),
            #  axis.text.y = element_text(size = 20),
            #  axis.title = element_text(size = 20),
            #  plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
              plot.background = element_rect(fill = "white", color = NA))
        print(p)
      dev.off()
    }
  }


#dotplot
library(ggplot2)
library(dplyr)
library(extrafont)
library(patchwork)
library(cluster)
library(Seurat)
library(pheatmap)
library(tidyr)
library(cowplot)
#MHC
gene_pos_df <- read.csv("/node1/liuxy/index/gene_pos_grch37.82.csv")
mhc_genes <- gene_pos_df$gene_name[gene_pos_df$chr %in% "6" & as.numeric(gene_pos_df$start) >= 28477797 & as.numeric(gene_pos_df$end) <= 33448354]

#pbmc
setwd("/node1/liuxy/MediCell/pbmc_2/TWAS/output/plot/")
#pbmc_twas <- read.csv("/node1/liuxy/chplot/scTWAS/usefile/pbmc_peer2/Bonfsig_peer2_allcelltype_result.csv")
pbmc_twas <- read.csv("/node1/liuxy/MediCell/pbmc_2/TWAS/output_check/pbmc/SLE/summary/Bonfsig_allcelltype_result.csv")
ct_order_pbmc <- c("CD4T","CD8T","NK","B","Mono","DC") 
#新的结果除了有少数p有极小值外effectsize没有膨胀
pbmc_twas$cell_type <- factor(pbmc_twas$cell_type, levels = rev(ct_order_pbmc))
pbmc_twas_filtered <- pbmc_twas[!pbmc_twas$gene_name %in% mhc_genes, ]
gene_pos_df$chr <- factor(gene_pos_df$chr, 
                          levels = c(1:22, "X", "Y", "MT", "M"))
gene_pos_filtered <- gene_pos_df %>%
  filter(gene_name %in% pbmc_twas_filtered$gene_name)
gene_pos_ordered <- gene_pos_filtered %>%
  arrange(chr, start)
gene_order <- gene_pos_ordered$gene_name

onek1k <- readRDS("/node1/liuxy/chplot/combind_plasma_B/cibersort/usefile/pbmc_6_CellType_971samples.rds")
cells_to_keep <- colnames(onek1k)[onek1k$cell_label %in% ct_order_pbmc]
genes_to_keep <- rownames(onek1k)[toupper(gsub("_.*", "", rownames(onek1k))) %in% pbmc_twas_filtered$gene_id]
expr_data <- GetAssayData(onek1k, layer = "data")[unique(pbmc_twas_filtered$gene_id), cells_to_keep]
celltype_factor <- factor(onek1k$cell_label[cells_to_keep], levels = ct_order_pbmc)
avg_expr <- sapply(levels(celltype_factor), function(ct) {
  rowMeans(expr_data[, celltype_factor == ct, drop = FALSE])
})
avg_expr_df <- as.data.frame(avg_expr)
rownames(avg_expr_df) <- gene_pos_df$gene_name[match(rownames(avg_expr_df), gene_pos_df$gene_id)]
avg_expr_scaled <- t(scale(t(avg_expr_df)))

twas_wide <- pbmc_twas_filtered %>%
  select(gene_name, cell_type, effect_size) %>%
  pivot_wider(names_from = cell_type, values_from = effect_size, values_fill = 0)
twas_wide <- twas_wide[match(gene_order, twas_wide$gene_name), ]
twas_matrix <- as.matrix(twas_wide[, -1])
rownames(twas_matrix) <- twas_wide$gene_name
twas_matrix <- twas_matrix[, ct_order_pbmc]
cell_types <- colnames(avg_expr_scaled)
genes <- rownames(avg_expr_scaled)

mean_expr_across_celltypes <- rowMeans(avg_expr_df) 
sd_expr_across_celltypes <- apply(avg_expr_df, 1, sd)
specific_genes_matrix <- matrix(FALSE, nrow = nrow(avg_expr_df), ncol = ncol(avg_expr_df),
                                dimnames = list(rownames(avg_expr_df), colnames(avg_expr_df)))
for(i in 1:nrow(avg_expr_df)) {
  for(j in 1:ncol(avg_expr_df)) {
    if(avg_expr_df[i,j] > (mean_expr_across_celltypes[i] + sd_expr_across_celltypes[i])) {
      specific_genes_matrix[i,j] = TRUE
    }
  }
}
specific_genes_long <- as.data.frame(as.table(specific_genes_matrix))
colnames(specific_genes_long) <- c("gene", "cell_type", "is_specific")
specific_genes_long <- specific_genes_long[specific_genes_long$is_specific, ]
specific_genes_long <- specific_genes_long[specific_genes_long$gene %in% gene_order & 
                                             specific_genes_long$cell_type %in% ct_order_pbmc, ]

expr_data_long <- data.frame(
  cell_type = rep(cell_types, each = length(genes)),
  gene = rep(genes, times = length(cell_types)),
  expression = as.vector(avg_expr_scaled))
expr_data_long$cell_type <- factor(expr_data_long$cell_type, levels = ct_order_pbmc)
expr_data_long$gene <- factor(expr_data_long$gene, levels = genes)
expr_data_long$gene <- factor(expr_data_long$gene, levels = unique(pbmc_twas_filtered$gene_name))

heatmap_data <- data.frame(
  cell_type = rep(colnames(avg_expr_scaled), each = nrow(avg_expr_scaled)),
  gene = rep(rownames(avg_expr_scaled), times = ncol(avg_expr_scaled)),
  zscore = as.vector(avg_expr_scaled)
)
heatmap_data$cell_type <- factor(heatmap_data$cell_type, levels = rev(ct_order_pbmc))
heatmap_data$gene <- factor(heatmap_data$gene, levels = gene_order)

write.table(gene_order, "gene_order_bonf_check.txt",col.names = FALSE,row.names = FALSE)

cell_count <- aggregate(gene_name ~ cell_type, data = pbmc_twas_filtered, FUN = function(x) length(unique(x)))
colnames(cell_count)[2] <- "gene_count"
ct_order_with_count <- paste0(ct_order_pbmc, " (n=", cell_count$gene_count[match(ct_order_pbmc, cell_count$cell_type)], ")")
labels_with_count <- paste0(cell_count$cell_type, " (n=", cell_count$count, ")")
names(labels_with_count) <- cell_count$cell_type
#bonfp plotsize
pdf("combined_heatmap_dotplot_bonfp_samesize_check.pdf", width = 35, height = 15)
TWAS_p <- ggplot(pbmc_twas_filtered, aes(
  x = gene_name,
  y = cell_type,
  fill = ifelse(effect_size > 0, "Positive", "Negative"),
  size = 20)) +
  geom_point(shape = 21, color = "#FFFFFF", size = 14, position = position_nudge(x = 0)) +
#  scale_x_discrete(limits = gene_order) +
  scale_x_discrete(limits = gene_order, expand = c(0.01, 0.01), drop = FALSE) +
  #scale_y_discrete(limits = rev(ct_order_pbmc)) +
  scale_y_discrete(limits = rev(ct_order_pbmc), labels = rev(ct_order_with_count)) +
  scale_fill_manual(
    name = "Effect direction",
    values = c("Positive" = "#f67d25", "Negative" = "#894297"),
    guide = guide_legend(
      override.aes = list(size = 6),
      order = 1
    )
  ) +
#  scale_size_continuous(
#    name = "-log10(p-value)",
#    range = c(min_size, max_size),
#    breaks = c(min_size, min_size + (max_size-min_size)*0.25, 
#               min_size + (max_size-min_size)*0.5, 
#               min_size + (max_size-min_size)*0.75, max_size),
#    labels = c("0", "3.75", "7.5", "11.25", "≥15"),
#    guide = guide_legend(
#      override.aes = list(fill = "#000000", shape = 21, color = "#FFFFFF"),
#      order = 2
#    )
#  ) +
  labs(title = "PBMC") +
  theme_bw(base_size = 24) +
  theme(
   # axis.text.x = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1, size=1, color="transparent"),
    axis.text.y = element_text(size = 20, margin = margin(r = 5)),
    #axis.text.y = element_text(size = 20, margin = margin(r = 5), label = function(y) {paste0(y, " (n=", label_map[y], ")")}),
    axis.title = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.spacing.x = unit(1, "cm"),
    plot.margin = margin(l = 1, b = 0)
  )
#heatmap
heatmap_p <- ggplot(heatmap_data, aes(x = gene, y = cell_type, fill = zscore)) +
 # geom_tile(color = "#FFFFFF", linewidth = 0.2, width = 1.2, height = 1.0) +
  geom_tile(color = "#FFFFFF", linewidth = 0.2, width = 1, height = 1.0) +
  geom_text(data = specific_genes_long, 
            aes(x = gene, y = cell_type, label = "*"), 
            color = "#000000", size = 8, inherit.aes = FALSE) +
  scale_fill_gradientn(
    colors = colorRampPalette(c("#894297","#FFFFFF","#f67d25"))(100),
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    name = "Expression Z-score") +
 # scale_x_discrete(limits = gene_order) +
  scale_x_discrete(limits = gene_order, expand = c(0.01, 0.01), drop = FALSE) +
  scale_y_discrete(limits = rev(ct_order_pbmc)) +
  theme_bw(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(t = 0, l = 5, r = 5)
  ) 
final_layout_detailed <- TWAS_p / heatmap_p +
  plot_layout(heights = c(3, 1), guides = "collect") &
  theme(legend.position = "bottom")
print(final_layout_detailed)
dev.off()

