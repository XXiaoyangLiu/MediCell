library(ggplot2)
library(dplyr)
library(extrafont)
library(patchwork)
library(cluster)
library(pheatmap)
library(tidyr)
library(textshape)
library(ggVennDiagram)

fill_colors <- c(
  common = "#E0A820",
  kidney = "#6088A0",
  glom = "#C8D786",
  tub = "#404828",
  pbmc = "#C84008"
)
gene_pos_df <- read.csv("/node1/liuxy/index/gene_pos_grch37.82.csv")
mhc_genes <- gene_pos_df$gene_name[gene_pos_df$chr %in% "6" & as.numeric(gene_pos_df$start) >= 28477797 & as.numeric(gene_pos_df$end) <= 33448354]
setwd("/node1/liuxy/MediCell/tissue_compare_2/output/")

pbmc_twas <- read.csv("/node1/liuxy/MediCell/pbmc_2/TWAS/output_check/pbmc/SLE/summary/FDRsig_allcelltype_result.csv")
custom_cell_order <- rev(c("DC","Mono","NK","CD4T","CD8T","B")) 
#pbmc_twas <- pbmc_twas[abs(pbmc_twas$effect_size) < 100,]
#pbmc_twas_filtered <- pbmc_twas[!pbmc_twas$gene_name %in% mhc_genes, ]
#glom
glom_twas <- read.csv("/node1/liuxy/MediCell/kidney_3/glom/TWAS/output/glom/SLE/summary/FDRsig_allcelltype_result.csv")
custom_cell_order <- rev(c("Endo","Podo","MC","MNP","NK","NKT","CD8T","Th_cell","Bcell")) 
#glom_twas <- glom_twas[abs(glom_twas$effect_size) < 100,]
#glom_twas_filtered <- glom_twas[!glom_twas$gene_name %in% mhc_genes, ]
#tub
tub_twas <- read.csv("/node1/liuxy/MediCell/kidney_3/tub/TWAS/output/tub/SLE/summary/FDRsig_allcelltype_result.csv")
custom_cell_order <- rev(c("PT","LOH","DCT","CDICA","CDICB","MNP","Plasmacytoid_DC","NK","NKT","CD8T","Th_cell","Bcell")) 
tub_twas[which(tub_twas$cell_type=='Plasmacytoid_DC'),'cell_type']='pDC'
#tub_twas <- tub_twas[abs(tub_twas$effect_size) < 100,]
#tub_twas_filtered <- tub_twas[!tub_twas$gene_name %in% mhc_genes, ]

pbmc_gene <- unique(pbmc_twas$gene_id)
glom_gene <- unique(glom_twas$gene_id)
glom_gene <- sub("\\..*", "", glom_gene)
tub_gene <- unique(tub_twas$gene_id)
tub_gene <- sub("\\..*", "", tub_gene)
kidney_gene <- unique(c(glom_twas$gene_id,tub_twas$gene_id))
kidney_gene <- sub("\\..*", "", kidney_gene)

pbmc_specific <- setdiff(pbmc_gene,kidney_gene)
kidney_specific <- setdiff(kidney_gene,pbmc_gene)
glom_specific <- setdiff(glom_gene,tub_gene)
tub_specific <- setdiff(tub_gene,glom_gene)
tissue_common <- intersect(kidney_gene,pbmc_gene)
max_length <- max(length(pbmc_specific), length(kidney_specific), length(tissue_common))

#venn
gene_list <- list(
  PBMC = unique(pbmc_gene),
  Glom = unique(glom_gene),
  Tub  = unique(tub_gene)
)
A <- gene_list$PBMC
B <- gene_list$Glom
C <- gene_list$Tub

only_A <- setdiff(A, union(B,C))
only_B <- setdiff(B, union(A,C))
only_C <- setdiff(C, union(A,B))
AB_only <- setdiff(intersect(A,B), C)
AC_only <- setdiff(intersect(A,C), B)
BC_only <- setdiff(intersect(B,C), A)
ABC_all <- intersect(intersect(A,B), C)

venn_all_genes <- list(
  pbmc_specific = only_A,
  glom_specific = only_B,
  tub_specific = only_C,
  pbmc_and_glom = AB_only,
  pbmc_and_tub = AC_only,
  glom_and_tub = BC_only,
  pbmc_glom_tub = ABC_all
)
sapply(venn_all_genes, length)
gene_df <- data.frame()
for(region in names(venn_all_genes)){
  temp <- data.frame(
    gene_id = venn_all_genes[[region]],
    venn_region = region
  )
  gene_df <- rbind(gene_df, temp)
}
idx <- match(gene_df$gene_id, gene_pos_df$gene_id)
gene_df$gene_name <- gene_pos_df$gene_name[idx]
write.csv(gene_df, file = "venn_all_genes_check.csv", row.names = FALSE)
pdf("tissue_venn_ENSG_fdr_check.pdf",width = 6,height = 6)
ggVennDiagram(
  gene_list, 
  label_alpha = 0, 
  label = "count",
  label_size = 10,
  category.names = c("PBMC", "Glom", "Tub")) +
  #scale_fill_gradient(low = "#E0A820", high = "#C8D786") +
  scale_fill_gradient(low = "#FFFFFF", high = "#89c760") +
  theme(
    legend.position = "none",
    text = element_text(size = 20))
dev.off()

#挑特定的term
#GO
#enrich
#enrich 结果来自venn
#20261114
pbmc_specific_func <- read.csv("/node1/liuxy/MediCell/tissue_compare_2/usefile/david_20261114/DAVIDFunctAnnotClusterReport_pbmc_specific_check_2026-11-14.csv")
glom_specific_func <- read.csv("/node1/liuxy/MediCell/tissue_compare_2/usefile/david_20261114/DAVIDFunctAnnotClusterReport_glom_specific_check_2026-11-14.csv")
tub_specific_func <- read.csv("/node1/liuxy/MediCell/tissue_compare_2/usefile/david_20261114/DAVIDFunctAnnotClusterReport_tub_specific_check_2026-11-14.csv")
kidney_total_func <- read.csv("/node1/liuxy/MediCell/tissue_compare_2/usefile/david_20261114/DAVIDFunctAnnotClusterReport_kidney_check_2026-11-14.csv")
common_func <- read.csv("/node1/liuxy/MediCell/tissue_compare_2/usefile/david_20261114/DAVIDFunctAnnotClusterReport_common_check_2026-11-14.csv")

pbmc_specific_func$Group <- "pbmc"
glom_specific_func$Group <- "glom" 
tub_specific_func$Group <- "tub"
kidney_total_func$Group <- "kidney"
common_func$Group <- "common"

pbmc_specific_func$Term <- paste0("pbmc_", pbmc_specific_func$Term)
glom_specific_func$Term <- paste0("glom_", glom_specific_func$Term)
tub_specific_func$Term <- paste0("tub_", tub_specific_func$Term)
kidney_total_func$Term <- paste0("kidney_", kidney_total_func$Term)
common_func$Term <- paste0("common_", common_func$Term)
combined_data <- rbind(
  pbmc_specific_func,
  glom_specific_func,
  tub_specific_func,
  kidney_total_func,
  common_func
)
select_term_pbmc <- c("leukocyte activation","leukocyte activation involved in immune response",
                        "response to lipid")
select_term_pbmc <- paste0("pbmc_",select_term_pbmc)
select_term_glom <- c("regulation of body fluid levels","blood coagulation","response to wounding",
                        "lymphocyte activation","leukocyte chemotaxis","leukocyte activation")
select_term_glom <- paste0("glom_",select_term_glom)
select_term_tub <- c("monoatomic ion homeostasis","monoatomic ion transmembrane transport","inorganic ion homeostasis",
                        "tube development","tube morphogenesis","angiogenesis","myeloid leukocyte migration")
select_term_tub <- paste0("tub_",select_term_tub)
select_term_kidney <- c("leukocyte activation involved in inflammatory response",
                        "phosphorus metabolic process","phosphorus metabolic process",
                        "purine-containing compound metabolic process","response to wounding","regulation of body fluid levels")
select_term_kidney <- paste0("kidney_",select_term_kidney)
select_term_common <- c("antigen processing and presentation of exogenous antigen",
                        "antigen processing and presentation of peptide antigen",
                        "antigen processing and presentation of peptide or polysaccharide antigen via MHC class II",
                        "immune response-activating signaling pathway",
                        "immune response-regulating signaling pathway")
select_term_common <- paste0("common_",select_term_common)
all_select <- c(select_term_common,select_term_kidney,select_term_glom,select_term_tub,select_term_pbmc)
select_data <- combined_data[combined_data$Term %in% all_select,]
select_data$log10P <- -log10(select_data$P.Value)

group_levels <- c("common", "kidney", "glom", "tub", "pbmc")
select_data_sorted <- select_data %>%
  mutate(Group = factor(Group, levels = group_levels)) %>%
  group_by(Group) %>%
  arrange(desc(log10P), -log(P.Value), desc(Count), .by_group = TRUE) %>%
  mutate(Term = factor(Term, levels = rev(unique(Term)))) %>%
  ungroup()

write.csv(select_data_sorted,"BP3_select_terms_FDR_sorted_20261114_check.csv")

pdf("tissue_bp3_fdr_select_20261114_check.pdf", width = 18, height = 12)
ggplot(select_data_sorted, aes(x = log10P, y = Term)) +
  geom_bar(aes(fill = Group), stat = "identity", alpha = 0.8) +
  geom_text(aes(label = paste0(Count, "/", List.Total)), 
            hjust = -0.1, size = 5) +
  labs(x = "-log10(P.Value)", fill = "Tissue") +
  scale_fill_manual(values = fill_colors) +
  theme_minimal() +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    plot.caption = element_text(size = 14, color = "gray40"),
    panel.spacing = unit(0.8, "lines")
  ) +
  expand_limits(x = max(select_data_sorted$log10P) * 1.3)
dev.off()
