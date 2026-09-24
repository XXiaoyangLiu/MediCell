library(stringr)
library(dplyr)
library(tidyr)
library(ggrepel)
library(ggplot2)
library(patchwork)
#manhattan
#allcelltype
gene_pos_df <- read.csv("/node1/liuxy/index/gene_pos_grch37.82.csv")
mhc_genes <- gene_pos_df$gene_name[gene_pos_df$chr %in% "6" & as.numeric(gene_pos_df$start) >= 28477797 & as.numeric(gene_pos_df$end) <= 33448354]
#pbmc
tissues <- "pbmc"
eGFR_base_paths <- "/node1/liuxy/MediCell/trait_compare_2/eGFR/output/TWAS/pbmc/eGFR/"
SLE_base_paths <- "/node1/liuxy/MediCell/pbmc_2/TWAS/output/pbmc/SLE/"
celltypes <- list("CD4T", "CD8T", "NK", "B", "Mono", "DC")
out_dir <- "/node1/liuxy/MediCell/trait_compare_2/eGFR/output/plot/manhattan/"
SLE_FDRsig <- read.csv("/node1/liuxy/MediCell/pbmc_2/TWAS/output/pbmc/SLE/summary/FDRsig_allcelltype_result.csv")
SLE_FDRsig <- SLE_FDRsig[!SLE_FDRsig$gene_name %in% mhc_genes,]
eGFR_FDRsig <- read.csv("/node1/liuxy/MediCell/trait_compare_2/eGFR/output/TWAS/pbmc/eGFR/summary/FDRsig_allcelltype_result.csv")
eGFR_FDRsig <- eGFR_FDRsig[!eGFR_FDRsig$gene_name %in% mhc_genes,]

for (ct in celltypes) {
  sle_file <- list.files(SLE_base_paths, pattern=paste0(ct,"\\.csv$"), full.names=TRUE)
  df_sle <- read.csv(sle_file)
  egfr_file <- list.files(eGFR_base_paths, pattern=paste0(ct,"\\.csv$"), full.names=TRUE)
  df_egfr <- read.csv(egfr_file)
  df_sle <- merge(df_sle, gene_pos_df, by.x="gene", by.y="gene_id")
  df_egfr <- merge(df_egfr, gene_pos_df, by.x="gene", by.y="gene_id")
  df_sle$gene_name <- df_sle$gene_name.y
  df_sle$gene_name.x <- NULL
  df_egfr$gene_name <- df_egfr$gene_name.y
  df_egfr$gene_name.x <- NULL
  df_sle <- df_sle %>% filter(chr %in% 1:22, !is.na(start))
  df_sle$trait <- "SLE"
  df_egfr <- df_egfr %>% filter(chr %in% 1:22, !is.na(start))
  df_egfr$trait <- "eGFR"
  
  fdr_threshold <- 0.05
  sle_sig_genes <- SLE_FDRsig %>%
    filter(fdr_p < fdr_threshold) %>%
    filter(cell_type %in% ct) %>%
    pull(gene_name) %>% unique()
  egfr_sig_genes <- eGFR_FDRsig %>%
    filter(fdr_p < fdr_threshold) %>%
    filter(cell_type %in% ct) %>%
    pull(gene_name) %>% unique()
  common_sig_genes <- intersect(sle_sig_genes, egfr_sig_genes)
  
  df <- rbind(df_sle, df_egfr)
  df <- df[, c("gene_name", "chr", "start", "end", "pvalue", "trait")]
  df$chr <- factor(as.integer(df$chr), levels = 1:22)
  cumu_df <- df %>%
    select(chr, start, end) %>%
    distinct() %>%
    filter(chr %in% 1:22) %>%
    mutate(position = as.numeric(start)) %>%
    mutate(chr = factor(chr, levels = 1:22)) %>%
    arrange(chr, position) %>%
    group_by(chr) %>%
    mutate(chr_min = min(position),
           chr_max = max(position),
           chr_span = chr_max - chr_min) %>%
    ungroup() %>%
    mutate(rel_pos = (position - chr_min) / (chr_span + 1),
           max_allowed_width = 1e5,
           actual_width = ifelse(chr_span > max_allowed_width, max_allowed_width, chr_span),
           compressed_rel_pos = ifelse(chr_span > max_allowed_width, rel_pos * (max_allowed_width / chr_span), rel_pos),
           chr_start = cumsum(lag(actual_width, default = 0) + 1e5),
           cumulative_pos = chr_start + compressed_rel_pos * actual_width)

  chr_color <- ifelse(1:22 %% 2 == 1, "#894297", "#FFC107")
  plot_df <- df %>%
    left_join(cumu_df, by = c("chr", "start", "end"))
  axis_df <- cumu_df %>%
    group_by(chr) %>%
    summarise(center = mean(cumulative_pos, na.rm=T), .groups="drop")
  plot_egfr <- plot_df %>% filter(trait == "eGFR")
  plot_sle  <- plot_df %>% filter(trait == "SLE")
  label_egfr <- plot_egfr %>%
    filter(gene_name %in% common_sig_genes, !gene_name %in% mhc_genes)
  label_sle  <- plot_sle %>%
    filter(gene_name %in% common_sig_genes, !gene_name %in% mhc_genes)
  
  out_pdf <- file.path(out_dir, paste0("SLE_eGFR_pbmc_",ct,"_FDRintersect.pdf"))
  pdf(out_pdf, width=20, height=10)
  
  p1 <- ggplot(plot_egfr, aes(x=cumulative_pos, y=-log10(pvalue))) +
    geom_point(aes(color=chr), size=2, alpha=1) +
    geom_segment(
      data = label_egfr,
      aes(x = cumulative_pos, xend = cumulative_pos,
          y = -log10(pvalue), yend = -log10(pvalue) + 0.8),
      linewidth = 0.3, color = "black"
    ) +
    scale_x_continuous(breaks=axis_df$center, labels=axis_df$chr, expand=c(0.01,0.01)) +
    scale_color_manual(values = chr_color) +
    labs(title = ct, y=expression(eGFR~~-log[10](p))) +
    theme_minimal() +
    theme(legend.position="none",
          axis.text.x=element_text(angle=90, size=16),
          axis.text.y=element_text(size=16),
          axis.title.y=element_text(size=18),
          plot.title=element_text(hjust=0.5, size=20, face="bold"))
  if(nrow(label_egfr) > 0){
    p1 <- p1 + geom_text_repel(
      data=label_egfr, 
      aes(label=gene_name), 
      size=4, 
      max.overlaps=30,
      nudge_y = 0.9,
      segment = FALSE,
      direction = "x"
    )
  }
  p2 <- ggplot(plot_sle, aes(x=cumulative_pos, y=-log10(pvalue))) +
    geom_point(aes(color=chr), size=2, alpha=1) +
    geom_segment(
      data = label_sle,
      aes(x = cumulative_pos, xend = cumulative_pos,
          y = -log10(pvalue), yend = -log10(pvalue) - 0.8),
      linewidth = 0.3, color = "black"
    ) +
    scale_x_continuous(breaks=axis_df$center, labels=axis_df$chr, expand=c(0.01,0.01)) +
    scale_y_reverse() +
    scale_color_manual(values = chr_color) +
    labs(y=expression(SLE~~-log[10](p))) +
    theme_minimal() +
    theme(legend.position="none",
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=16),
          axis.title=element_text(size=18),
          plot.title=element_text(hjust=0.5, size=20, face="bold"))
  if(nrow(label_sle) > 0){
    p2 <- p2 + geom_text_repel(
      data=label_sle, 
      aes(label=gene_name), 
      size=4, 
      max.overlaps=20,
      nudge_y = -0.9,
      segment = FALSE,
      direction = "x"
    )
  }
  combined <- wrap_plots(p1, p2, ncol=1, heights=c(1,1))
  print(combined)
  dev.off()
}
