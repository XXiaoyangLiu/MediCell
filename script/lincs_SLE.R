
library(dplyr)
library(tidyr)
library(textshape)
library(purrr)
library(tibble)
library(ggplot2)    
setwd("/node1/liuxy/MediCell/lincs_2/usefile/")
ct_list_pbmc <- c("B", "CD4T", "CD8T", "DC", "Mono", "NK")
ct_list_glom <- c("Endo","Podo","MC","MNP","NK","NKT","CD8T","Th_cell","Bcell")
ct_list_tub<- c("PT","LOH","DCT","CDICA","CDICB","MNP","Plasmacytoid_DC","NK","NKT","CD8T","Bcell","Th_cell")
#exp
#pbmc
exp_dir <- "/node1/liuxy/chplot/combind_plasma_B/enigma/output/count/tpm0.1sample0.2/"
ct_list <- c("B", "CD4T", "CD8T", "DC", "Mono", "NK")
ct_expZ_list <- list()
for(ct in ct_list){
  fn <- paste0(exp_dir, ct, "_transformed_expression.txt")
  raw <- read.table(fn, header = T)
  raw_t <- raw %>%
    column_to_rownames(colnames(raw)[1]) %>%
    t() %>%
    as.data.frame()
  gene_z <- apply(raw_t, 1, function(gene_vec){
    mu <- mean(gene_vec, na.rm = TRUE)
    sigma <- sd(gene_vec, na.rm = TRUE)
    if(sigma == 0) return(0) 
    (gene_vec[1] - mu)/sigma
  })
  ct_expZ_list[[ct]] <- tibble(
    gene_name = names(gene_z),
    !!paste0(ct,"_expZ") := gene_z
  )
}
expZ_mat <- reduce(ct_expZ_list, full_join, by = "gene_name") %>%
  replace(is.na(.), 0) %>% 
  column_to_rownames("gene_name")
write.csv(expZ_mat, "/node1/liuxy/MediCell/lincs_2/usefile/exp/pbmc_exp.txt")
save(expZ_mat, file = "/node1/liuxy/MediCell/lincs_2/usefile/exp/pbmc_exp.RData")

#glom
exp_dir <- "/node1/liuxy/chplot/kidney_new/output/glom/exp/"
ct_list <- c("Endo","Podo","MC","MNP","NK","NKT","CD8T","Th_cell","Bcell")
ct_expZ_list <- list()
for(ct in ct_list){
  fn <- paste0(exp_dir, ct, "_transformed_expression.txt")
  raw <- read.table(fn, header = T)
  raw_t <- raw %>%
    column_to_rownames(colnames(raw)[1]) %>%
    t() %>%
    as.data.frame()
  gene_z <- apply(raw_t, 1, function(gene_vec){
    mu <- mean(gene_vec, na.rm = TRUE)
    sigma <- sd(gene_vec, na.rm = TRUE)
    if(sigma == 0) return(0) 
    (gene_vec[1] - mu)/sigma
  })
  ct_expZ_list[[ct]] <- tibble(
    gene_name = names(gene_z),
    !!paste0(ct,"_expZ") := gene_z
  )
}
expZ_mat <- reduce(ct_expZ_list, full_join, by = "gene_name") %>%
  replace(is.na(.), 0) %>% 
  column_to_rownames("gene_name")
write_tsv(rownames_to_column(expZ_mat,"gene_name"), "/node1/liuxy/MediCell/lincs_2/usefile/exp/glom_exp.txt")
save(expZ_mat, file = "/node1/liuxy/MediCell/lincs_2/usefile/exp/glom_exp.RData")
#tub
exp_dir <- "/node1/liuxy/chplot/kidney_new/output/tub/exp/"
ct_list <- c("PT","LOH","DCT","CDICA","CDICB","MNP","Plasmacytoid_DC","NK","NKT","CD8T","Bcell","Th_cell")
ct_expZ_list <- list()
for(ct in ct_list){
  fn <- paste0(exp_dir, ct, "_transformed_expression.txt")
  raw <- read.table(fn, header = T)
  raw_t <- raw %>%
    column_to_rownames(colnames(raw)[1]) %>%
    t() %>%
    as.data.frame()
  gene_z <- apply(raw_t, 1, function(gene_vec){
    mu <- mean(gene_vec, na.rm = TRUE)
    sigma <- sd(gene_vec, na.rm = TRUE)
    if(sigma == 0) return(0) 
    (gene_vec[1] - mu)/sigma
  })
  ct_expZ_list[[ct]] <- tibble(
    gene_name = names(gene_z),
    !!paste0(ct,"_expZ") := gene_z
  )
}
expZ_mat <- reduce(ct_expZ_list, full_join, by = "gene_name") %>%
  replace(is.na(.), 0) %>% 
  column_to_rownames("gene_name")
write.csv(expZ_mat, "/node1/liuxy/MediCell/lincs_2/usefile/exp/tub_exp.txt")
save(expZ_mat, file = "/node1/liuxy/MediCell/lincs_2/usefile/exp/tub_exp.RData")

#TWAS
pbmc_SLE_FDR <- read.csv("/node1/liuxy/MediCell/pbmc_2/TWAS/output/pbmc/SLE/summary/FDRsig_allcelltype_result.csv")
pbmc_SLE_FDR <- pbmc_SLE_FDR[pbmc_SLE_FDR$zscore > 0,]
glom_SLE_FDR <- read.csv("/node1/liuxy/MediCell/kidney_3/glom/TWAS/output/glom/SLE/summary/FDRsig_allcelltype_result.csv")
glom_SLE_FDR <- glom_SLE_FDR[glom_SLE_FDR$zscore > 0,]
tub_SLE_FDR <- read.csv("/node1/liuxy/MediCell/kidney_3/tub/TWAS/output/tub/SLE/summary/FDRsig_allcelltype_result.csv")
tub_SLE_FDR <- tub_SLE_FDR[tub_SLE_FDR$zscore > 0,]

#lincs
lincs <- read.csv("/node1/liuxy/MediCell/lincs_2/usefile/LINCS_20260917_sig_results_minsize10_SLEpos.csv",header=TRUE)
lincs <- lincs[lincs$pert_type %in% "trt_cp",]
lincs <- na.omit(lincs)
lincs <- lincs[!lincs$pubchem_cid %in% "-666",]
pbmc_lincs  <- lincs[lincs$TWAS %in% "SLE_pbmc_pos",]
glom_lincs  <- lincs[lincs$TWAS %in% "SLE_glom_pos",]
tub_lincs   <- lincs[lincs$TWAS %in% "SLE_tub_pos",]

pbmcShare <- pbmc_lincs[!duplicated(pbmc_lincs$pert_iname),]
pbmcShare <- pbmcShare[order(pbmcShare$pval),]
pbmcgenes <- pbmcShare[,c('pert_iname','leadingEdge')]
glomShare <- glom_lincs[!duplicated(glom_lincs$pert_iname),]
glomShare <- glomShare[order(glomShare$pval),]
glomgenes <- glomShare[,c('pert_iname','leadingEdge')]
tubShare <- tub_lincs[!duplicated(tub_lincs$pert_iname),]
tubShare <- tubShare[order(tubShare$pval),]
tubgenes <- tubShare[,c('pert_iname','leadingEdge')]

write.table(tubgenes,file="tubSLE_CommonGenes.txt",sep='\t',row.name=F,col.name=F,quote=F);
write.table(glomgenes,file="glomSLE_CommonGenes.txt",sep='\t',row.name=F,col.name=F,quote=F);
write.table(pbmcgenes,file="pbmcSLE_CommonGenes.txt",sep='\t',row.name=F,col.name=F,quote=F);

tub_leadingEdge <- tubgenes %>%
  separate_rows(leadingEdge, sep = ",") %>%
  rename(gene = leadingEdge) %>%
  mutate(gene = str_trim(gene)) %>%
  filter(gene!="")

glom_leadingEdge <- glomgenes %>%
  separate_rows(leadingEdge, sep = ",") %>%
  rename(gene = leadingEdge) %>%
  mutate(gene = str_trim(gene)) %>%
  filter(gene!="")

pbmc_leadingEdge <- pbmcgenes %>%
  separate_rows(leadingEdge, sep = ",") %>%
  rename(gene = leadingEdge) %>%
  mutate(gene = str_trim(gene)) %>%
  filter(gene!="")

le_gene_pbmc <- unique(pbmc_leadingEdge$gene)
le_gene_glom <- unique(glom_leadingEdge$gene)
le_gene_tub  <- unique(tub_leadingEdge$gene)
pbmc_twas_filtered <- pbmc_SLE_FDR %>%
  filter(gene_name %in% le_gene_pbmc) %>%
  select(gene_id, gene_name, cell_type, zscore)

pbmc_twas_wide <- pbmc_twas_filtered %>%
  pivot_wider(
    names_from = cell_type,
    values_from = zscore,
    values_fill = 0
  ) %>%
  rename_with(.fn = ~ paste0(.x, "_twasZ"), .cols = -c(gene_id, gene_name))

load("/node1/liuxy/MediCell/lincs_2/usefile/exp/pbmc_exp.RData")
pbmc_exp <- expZ_mat
pbmc_exp_df <- as.data.frame(pbmc_exp) %>% rownames_to_column("gene_id")

gene_comb_pbmc <- inner_join(pbmc_twas_wide, pbmc_exp_df, by = "gene_id")
for (ct in ct_list_pbmc) {
  col_twas <- paste0(ct, "_twasZ")
  col_exp  <- paste0(ct, "_expZ")
  gene_comb_pbmc[[paste0(ct, "_rank")]] <- gene_comb_pbmc[[col_twas]] * gene_comb_pbmc[[col_exp]]
}
gene_rank_pbmc <- gene_comb_pbmc %>%
  select(gene_id, gene_name, ends_with("_rank"))

pbmc_full <- left_join(
  pbmc_leadingEdge,
  gene_rank_pbmc,
  by = c("gene" = "gene_name")
)
pbmc_full_score <- pbmc_full %>%
  group_by(pert_iname) %>%
  summarise(
    across(ends_with("_rank"), sum, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(pbmc_full_score, "pbmcSLE_leadingedge_heatmap.csv",row.names = F)

#===================== 4. GLOM SLE =====================
glom_twas_filtered <- glom_SLE_FDR %>%
  filter(gene_name %in% le_gene_glom) %>%
  select(gene_id, gene_name, cell_type, zscore)

glom_twas_wide <- glom_twas_filtered %>%
  pivot_wider(
    names_from = cell_type,
    values_from = zscore,
    values_fill = 0
  ) %>%
  rename_with(.fn = ~ paste0(.x, "_twasZ"), .cols = -c(gene_id, gene_name))

load("/node1/liuxy/MediCell/lincs_2/usefile/exp/glom_exp.RData")
glom_exp <- expZ_mat
glom_exp_df <- as.data.frame(glom_exp) %>% rownames_to_column("gene_id")

gene_comb_glom <- inner_join(glom_twas_wide, glom_exp_df, by = "gene_id")
for (ct in ct_list_glom) {
  col_twas <- paste0(ct, "_twasZ")
  col_exp  <- paste0(ct, "_expZ")
  gene_comb_glom[[paste0(ct, "_rank")]] <- gene_comb_glom[[col_twas]] * gene_comb_glom[[col_exp]]
}
gene_rank_glom <- gene_comb_glom %>%
  select(gene_id, gene_name, ends_with("_rank"))

glom_full <- left_join(
  glom_leadingEdge,
  gene_rank_glom,
  by = c("gene" = "gene_name")
)
glom_full_score <- glom_full %>%
  group_by(pert_iname) %>%
  summarise(
    across(ends_with("_rank"), sum, na.rm = TRUE), 
    .groups = "drop"
  )
write.csv(glom_full_score, "glomSLE_leadingedge_heatmap.csv",row.names = F)
#===================== 5. TUB =====================
tub_twas_filtered <- tub_SLE_FDR %>%
  filter(gene_name %in% le_gene_tub) %>%
  select(gene_id, gene_name, cell_type, zscore)

tub_twas_wide <- tub_twas_filtered %>%
  pivot_wider(
    names_from = cell_type,
    values_from = zscore,
    values_fill = 0
  ) %>%
  rename_with(.fn = ~ paste0(.x, "_twasZ"), .cols = -c(gene_id, gene_name))

load("/node1/liuxy/MediCell/lincs_2/usefile/exp/tub_exp.RData")
tub_exp <- expZ_mat
tub_exp_df <- as.data.frame(tub_exp) %>% rownames_to_column("gene_id")

gene_comb_tub <- inner_join(tub_twas_wide, tub_exp_df, by = "gene_id")
for (ct in ct_list_tub) {
  col_twas <- paste0(ct, "_twasZ")
  col_exp  <- paste0(ct, "_expZ")
  gene_comb_tub[[paste0(ct, "_rank")]] <- gene_comb_tub[[col_twas]] * gene_comb_tub[[col_exp]]
}
gene_rank_tub <- gene_comb_tub %>%
  select(gene_id, gene_name, ends_with("_rank"))

tub_full <- left_join(
  tub_leadingEdge,
  gene_rank_tub,
  by = c("gene" = "gene_name")
)
tub_full_score <- tub_full %>%
  group_by(pert_iname) %>%
  summarise(
    across(ends_with("_rank"), sum, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(tub_full_score, "tubSLE_leadingedge_heatmap.csv",row.names = F)

# PBMC SLE heatmap
pbmc_result <- read.csv("pbmcSLE_leadingedge_heatmap.csv")
heatmap_long <- pbmc_result %>%
  pivot_longer(
    cols = -pert_iname,
    names_to = "cell_type",
    values_to = "score"
  ) %>%
  mutate(cell_type = sub("_rank$", "", cell_type))
drug_order <- pbmc_result$pert_iname

pdf("/node1/liuxy/MediCell/lincs_2/output/pbmcSLE_zscore_heatmap.pdf", width = 10, height = 14)
ggplot(heatmap_long, aes(x = cell_type, y = pert_iname, fill = score)) +
  geom_tile(color = "#FFFFFF", linewidth = 0.2, width = 1, height = 1) +
  scale_fill_gradientn(
    colors = colorRampPalette(c("#246906","#ffeb3d","#b2182b"))(100),
    oob = scales::squish,
    name = "Sum(TWASz × ExpZ)"
  ) +
  scale_x_discrete(limits = ct_list_pbmc, drop = FALSE) +
  scale_y_discrete(limits = drug_order, drop = FALSE) +
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
    plot.margin = margin(t = 10, l = 10, r = 10, b = 10)
  )
dev.off()

# Glom SLE heatmap
glom_result <- read.csv("glomSLE_leadingedge_heatmap.csv")
heatmap_long <- glom_result %>%
  pivot_longer(
    cols = -pert_iname,
    names_to = "cell_type",
    values_to = "score"
  ) %>%
  mutate(cell_type = sub("_rank$", "", cell_type))
drug_order <- glom_result$pert_iname

pdf("/node1/liuxy/MediCell/lincs_2/output/glomSLE_zscore_heatmap.pdf", width = 10, height = 14)
ggplot(heatmap_long, aes(x = cell_type, y = pert_iname, fill = score)) +
  geom_tile(color = "#FFFFFF", linewidth = 0.2, width = 1, height = 1) +
  scale_fill_gradientn(
    colors = colorRampPalette(c("#246906","#ffeb3d","#b2182b"))(100),
    oob = scales::squish,
    name = "Sum(TWASz × ExpZ)"
  ) +
  scale_x_discrete(limits = ct_list_glom, drop = FALSE) +
  scale_y_discrete(limits = drug_order, drop = FALSE) +
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
    plot.margin = margin(t = 10, l = 10, r = 10, b = 10)
  )
dev.off()

# Tub SLE heatmap
tub_result <- read.csv("tubSLE_leadingedge_heatmap.csv")
heatmap_long <- tub_result %>%
  pivot_longer(
    cols = -pert_iname,
    names_to = "cell_type",
    values_to = "score"
  ) %>%
  mutate(cell_type = sub("_rank$", "", cell_type))
drug_order <- tub_result$pert_iname

pdf("/node1/liuxy/MediCell/lincs_2/output/tubSLE_zscore_heatmap.pdf", width = 10, height = 18)
ggplot(heatmap_long, aes(x = cell_type, y = pert_iname, fill = score)) +
  geom_tile(color = "#FFFFFF", linewidth = 0.2, width = 1, height = 1) +
  scale_fill_gradientn(
    colors = colorRampPalette(c("#246906","#ffeb3d","#b2182b"))(100),
    oob = scales::squish,
    name = "Sum(TWASz × ExpZ)"
  ) +
  scale_x_discrete(limits = ct_list_tub, drop = FALSE) +
  scale_y_discrete(limits = drug_order, drop = FALSE) +
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
    plot.margin = margin(t = 10, l = 10, r = 10, b = 10)
  )
dev.off()
