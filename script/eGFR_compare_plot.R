library(VennDiagram)
library(grid)
library(gridExtra)
library(dplyr)
library(ggplot2)
#mhc
gene_pos_df <- read.csv("/node1/liuxy/index/gene_pos_grch37.82.csv")
mhc_genes <- gene_pos_df$gene_id[gene_pos_df$chr %in% "6" & as.numeric(gene_pos_df$start) >= 28477797 & as.numeric(gene_pos_df$end) <= 33448354]
ct_order_pbmc <- c("CD4T","CD8T","NK","B","Mono","DC") 
#pbmc
setwd('/node1/liuxy/MediCell/trait_compare_2/eGFR/output/plot/')
pbmc_twas_SLE <- read.csv("/node1/liuxy/MediCell/pbmc_2/TWAS/output/pbmc/SLE/summary/FDRsig_allcelltype_result.csv")
pbmc_twas_SLE$cell_type <- factor(pbmc_twas_SLE$cell_type, levels = rev(ct_order_pbmc))
pbmc_twas_SLE_filtered <- pbmc_twas_SLE[!pbmc_twas_SLE$gene_id %in% mhc_genes, ]
gene_pos_df$chr <- factor(gene_pos_df$chr, levels = c(1:22, "X", "Y", "MT", "M"))
gene_pos_SLE_filtered <- gene_pos_df %>%
  filter(gene_id %in% pbmc_twas_SLE_filtered$gene_id)
gene_pos_SLE_ordered <- gene_pos_SLE_filtered %>%
  arrange(chr, start)
gene_order_SLE <- gene_pos_SLE_ordered$gene_id

pbmc_twas_eGFR <- read.csv("/node1/liuxy/MediCell/trait_compare_2/eGFR/output/TWAS/pbmc/eGFR/summary/FDRsig_allcelltype_result.csv") 
pbmc_twas_eGFR$cell_type <- factor(pbmc_twas_eGFR$cell_type, levels = rev(ct_order_pbmc))
pbmc_twas_eGFR_filtered <- pbmc_twas_eGFR[!pbmc_twas_eGFR$gene_id %in% mhc_genes, ]
gene_pos_df$chr <- factor(gene_pos_df$chr, levels = c(1:22, "X", "Y", "MT", "M"))
gene_pos_eGFR_filtered <- gene_pos_df %>%
  filter(gene_id %in% pbmc_twas_eGFR_filtered$gene_id)
gene_pos_eGFR_ordered <- gene_pos_eGFR_filtered %>%
  arrange(chr, start)
gene_eGFR_order <- gene_pos_eGFR_ordered$gene_id

for(ct in ct_order_pbmc) {
  assign(paste0("set_SLE_", ct), 
         unique(pbmc_twas_SLE_filtered$gene_id[pbmc_twas_SLE_filtered$cell_type == ct]))
  assign(paste0("set_eGFR_", ct), 
         unique(pbmc_twas_eGFR_filtered$gene_id[pbmc_twas_eGFR_filtered$cell_type == ct]))
}
for(ct in ct_order_pbmc) {
  SLE_set <- get(paste0("set_SLE_", ct))
  eGFR_set <- get(paste0("set_eGFR_", ct))
  common <- intersect(SLE_set, eGFR_set)
  assign(paste0("common_", ct), common)
}
pdf("PBMC_SLE_vs_eGFR_celltype_Venn.pdf", width = 12, height = 16)
grid.newpage()
pushViewport(viewport(layout = grid.layout(4, 2)))
plot_positions <- list(CD4T = c(1, 1), CD8T = c(1, 2), NK = c(2, 1), B = c(2, 2), Mono = c(3, 1), DC = c(3, 2))
for(ct in ct_order_pbmc) {
  pos <- plot_positions[[ct]]
  pushViewport(viewport(layout.pos.row = pos[1], layout.pos.col = pos[2]))
  SLE_set <- get(paste0("set_SLE_", ct))
  eGFR_set <- get(paste0("set_eGFR_", ct))
  common_set <- get(paste0("common_", ct))
  if(length(SLE_set) == 0 & length(eGFR_set) == 0) {
    grid.text(paste(ct, "\n(No significant genes)"), 
              gp = gpar(fontsize = 10))
  } else {
    venn_plot <- draw.pairwise.venn(
      area1 = length(SLE_set),
      area2 = length(eGFR_set),
      cross.area = length(common_set),
      category = c("SLE", "eGFR"),
      fill = c("#894297", "#f67d25"),
      lty = "blank",
      cex = 1.5,
      cat.cex = 1,
      cat.pos = c(0, 0),
      cat.dist = 0.08,
      cat.just = list(c(0.5, 1), c(0.5, 1)),
      ext.pos = 0,
      ext.line.lwd = 0,
      euler.d = TRUE,
      scaled = FALSE,
      inverted = FALSE,
      rotation.degree = 0
    )
  }
  grid.text(ct, x = 0.5, y = 0.92, gp = gpar(fontsize = 12, fontface = "bold"))
  popViewport()
}
dev.off()

statistics_table <- data.frame(
  SLE_TWAS_Genes = sapply(ct_order_pbmc, function(ct) length(get(paste0("set_SLE_", ct)))),
  eGFR_TWAS_Genes = sapply(ct_order_pbmc, function(ct) length(get(paste0("set_eGFR_", ct)))),
  Common_Genes = sapply(ct_order_pbmc, function(ct) length(get(paste0("common_", ct)))),
  Overlap_Percentage_SLE = sapply(ct_order_pbmc, function(ct) {
    SLE_n <- length(get(paste0("set_SLE_", ct)))
    common_n <- length(get(paste0("common_", ct)))
    if(SLE_n > 0) round(common_n/SLE_n*100, 1) else 0
  }),
  Overlap_Percentage_eGFR = sapply(ct_order_pbmc, function(ct) {
    eGFR_n <- length(get(paste0("set_eGFR_", ct)))
    common_n <- length(get(paste0("common_", ct)))
    if(eGFR_n > 0) round(common_n/eGFR_n*100, 1) else 0
  })
)
write.csv(statistics_table, "PBMC_SLE_vs_eGFR_statistics.csv", row.names = FALSE)

#enrich
NK_SLE_func <- read.csv("/node1/liuxy/MediCell/trait_compare/eGFR/usefile/unique/DAVIDChartReport_NK_SLE_unique_2026-05-14.csv")
NK_eGFR_func <- read.csv("/node1/liuxy/MediCell/trait_compare/eGFR/usefile/unique/DAVIDChartReport_NK_eGFR_unique_2026-05-14.csv")
DC_SLE_func <- read.csv("/node1/liuxy/MediCell/trait_compare/eGFR/usefile/unique/DAVIDChartReport_DC_SLE_unique_2026-05-14.csv")
DC_eGFR_func <- read.csv("/node1/liuxy/MediCell/trait_compare/eGFR/usefile/unique/DAVIDChartReport_DC_eGFR_unique_2026-05-14.csv")

pdf("NK_SLE_unique_gobp3.pdf", width = 18, height = 8)
ggplot(NK_SLE_func, aes(x = -log(P.Value), y = Term)) +
  geom_bar(stat = "identity", fill = "#894297", alpha = 0.8) +
  geom_text(aes(label = paste0(Count, "/", List.Total)), 
            hjust = -0.1, size = 5) +
  labs(x = "-log(P.Value)", y = "") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "none",
    panel.spacing = unit(0.8, "lines")
  ) +
  expand_limits(x = max(-log(NK_SLE_func$P.Value)) * 1.3)
dev.off()

library(ggplot2)
library(dplyr)

selected_eGFR_terms <- c(
  "platelet activation",
  "blood coagulation",
  "response to wounding",
  "regulation of body fluid levels",
  "response to virus"
)

NK_eGFR_selected <- NK_eGFR_func %>%
  filter(Term %in% selected_eGFR_terms) %>%
  mutate(Group = "eGFR") 

NK_SLE_func <- NK_SLE_func %>%
  mutate(Group = "SLE")
combined_NK <- bind_rows(NK_eGFR_selected, NK_SLE_func)

combined_NK$Term <- factor(
  combined_NK$Term,
  levels = c(
    rev(NK_eGFR_selected$Term),
    rev(NK_SLE_func$Term) 
  )
)

pdf("NK_eGFR_vs_SLE_enrichment.pdf", width = 6, height = 3)
ggplot(combined_NK, aes(x = -log10(P.Value), y = Term)) +
  geom_bar(aes(fill = Group), stat = "identity", alpha = 0.8) +
  geom_text(aes(label = paste0(Count, "/", List.Total)), 
            hjust = -0.1, size = 4) +
  
  scale_fill_manual(values = c(
    "eGFR" = "#f67d25", 
    "SLE"  = "#894297"
  )) +
  
  labs(x = "-log10(P-value)", y = "") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  expand_limits(x = max(-log10(combined_NK$P.Value)) * 1.3)

dev.off()

library(ggplot2)
library(dplyr)

sle_terms <- c(
  "immune response-activating signaling pathway",
  "immune response-regulating signaling pathway",
  "regulation of immune system process",
  "adaptive immune response",
  "defense response"
)
DC_SLE_selected <- DC_SLE_func %>%
  filter(Term %in% sle_terms) %>%
  mutate(Group = "SLE",
         Term_Label = paste0(Term, " (SLE)")) 

egfr_terms <- c(
  "response to wounding",
  "cellular response to stress",
  "platelet activation",
  "blood coagulation",
  "regulation of immune system process"
)
DC_eGFR_selected <- DC_eGFR_func %>%
  filter(Term %in% egfr_terms) %>%
  mutate(Group = "eGFR",
         Term_Label = paste0(Term, " (eGFR)")) 

combined_DC <- bind_rows(DC_eGFR_selected, DC_SLE_selected)

combined_DC$Term_Label <- factor(
  combined_DC$Term_Label,
  levels = c(
    rev(DC_eGFR_selected$Term_Label),
    rev(DC_SLE_selected$Term_Label)
  )
)

pdf("DC_eGFR_vs_SLE_enrichment.pdf", width = 10, height = 6)

ggplot(combined_DC, aes(x = -log10(P.Value), y = Term_Label)) +
  geom_bar(aes(fill = Group), stat = "identity", alpha = 0.8) +
  geom_text(aes(label = paste0(Count, "/", List.Total)), 
            hjust = -0.1, size = 4) +
  scale_fill_manual(values = c(
    "eGFR" = "#f67d25", 
    "SLE"  = "#894297"
  )) +
  labs(x = "-log10(P-value)", y = "") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  expand_limits(x = max(-log10(combined_DC$P.Value)) * 1.3)

dev.off()

for(ct in ct_order_pbmc) {
  SLE_set    <- get(paste0("set_SLE_", ct))
  eGFR_set   <- get(paste0("set_eGFR_", ct))
  common_set <- get(paste0("common_", ct))
  
  only_SLE  <- setdiff(SLE_set, common_set)
  only_eGFR <- setdiff(eGFR_set, common_set)
  get_gene_info <- function(gid_vec, group_name){
    if(length(gid_vec)==0) return(NULL)
    gene_pos_df %>%
      filter(gene_id %in% gid_vec) %>%
      select(gene_id, gene_name) %>%
      mutate(Cell_Type = ct, Group = group_name) %>%
      relocate(Cell_Type, Group, gene_id, gene_name)
  }
  
  df_c <- get_gene_info(common_set, "Both_SLE_eGFR")
  df_s <- get_gene_info(only_SLE, "Only_SLE_FDRsig")
  df_e <- get_gene_info(only_eGFR, "Only_eGFR_FDRsig")
  
  total_class_df <- dplyr::bind_rows(df_c, df_s, df_e)
  
    write.csv(total_class_df, paste0("PBMC_", ct, "_gene_category_withName.csv"), row.names = FALSE)
}
