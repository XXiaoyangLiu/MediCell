# R4.1.3
library(ashr)              
library(mashr)              
library(dplyr)  
library(plyr) 
library(data.table)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(ggplot2)
library(RhpcBLASctl)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
cell_color_palettes <- list(
  pbmc = c(
    "DC" = "#fee301",
    "Mono" = "#f67d25",
    "NK" = "#89c760",
    "CD4T" = "#5b2d7c",
    "CD8T" = "#9983bd",
    "B" = "#bf6cac",
    "Plasma" = "#ae2f22"
  ),
  glom = c(
    "Endo" = "#4daf4a",
    "Podo" = "#1b9e77",
    "MC" = "#377eb8",
    "MNP" = "#f67d25",
    "NK" = "#89c760",
    "NKT" = "#9e5cb5",
    "CD8T" = "#9983bd",
    "Th_cell" = "#5b2d7c",
    "Bcell" = "#bf6cac"
  ),
  tub = c(
    "PT" = "#66c2a5",
    "LOH" = "#7fc97f",
    "DCT" = "#00c0a3",
    "CDICA" = "#084081",
    "CDICB" = "#0868ac",
    "MNP" = "#f67d25",
    "pDC" = "#fee301",
    "NK" = "#89c760",
    "NKT" = "#9e5cb5",
    "CD8T" = "#9983bd",
    "Bcell" = "#bf6cac",
    "Th_cell" = "#5b2d7c"
  )
)

setwd("/node1/liuxy/MediCell/tissue_compare_2/output/")
gene_pos_df <- read.csv("/node1/liuxy/index/gene_pos_grch37.82.csv")
mhc_genes <- gene_pos_df$gene_name[gene_pos_df$chr %in% "6" & as.numeric(gene_pos_df$start) >= 28477797 & as.numeric(gene_pos_df$end) <= 33448354]
tissue_dirs <- c("/node1/liuxy/MediCell/pbmc_2/TWAS/output/pbmc/SLE/", "/node1/liuxy/MediCell/kidney_3/glom/TWAS/output/glom/SLE/", "/node1/liuxy/MediCell/kidney_3/tub/TWAS/output/tub/SLE/")
#tissue_dirs <- c("/node1/liuxy/MediCell/pbmc/fulldb_TWAS/output/TWAS/pbmc/SLE/")

twas_result <- list()
for (tissue_path in tissue_dirs) {
  tissue <- basename(dirname(tissue_path))
  csv_files <- list.files(path = tissue_path, pattern = "\\.csv$", full.names = TRUE)
  csv_files <- csv_files[!grepl("summary", csv_files)]
  for (csv_file in csv_files) {
    celltype <- gsub("\\.csv$", "", basename(csv_file))
    tissue_celltype <- paste0(tissue, "_", celltype)
    data <- read.csv(csv_file)
    length(data$gene)
 #   data <- data[!data$gene_name %in% mhc_genes,]
#    data <- data[abs(data$effect_size) < 100 & data$pvalue > 1e-20 & data$pvalue <= 1, ]
    data$se <- sqrt(((data$effect_size)^2) / qchisq(data$pvalue, 1, lower.tail = FALSE))
    twas_result[[tissue_celltype]] <- data
  }
}
all_genes <- unique(unlist(lapply(twas_result, function(x) x$gene)))
tissue_celltype <- names(twas_result)

Bhat <- matrix(NA, nrow=length(all_genes), ncol=length(tissue_celltype), dimnames=list(all_genes, tissue_celltype))
Shat <- matrix(NA, nrow=length(all_genes), ncol=length(tissue_celltype), dimnames=list(all_genes, tissue_celltype))

for (i in seq_along(tissue_celltype)) {
  df <- twas_result[[i]]
  idx <- match(df$gene, rownames(Bhat))
  Bhat[idx, i] <- df$zscore 
  Shat[idx, i] <- df$se
}

data = mash_set_data(Bhat, Shat, alpha = 0, zero_check_to = 0)
m.1by1 = mash_1by1(data, alpha = 0)
lfsr_mat = get_lfsr(m.1by1)
keep = rowSums(lfsr_mat < 0.05, na.rm = TRUE) >= 1  

#keep <- rowSums(is.na(Bhat)) == 0
Bhat <- Bhat[keep, ]
Shat <- Shat[keep, ]
Vhat <- diag(ncol(Bhat))

data = mash_set_data(Bhat, Shat,alpha = 0, zero_check_to = 0)
m.1by1 = mash_1by1(data, alpha = 0)
strong.idx = get_significant_results(m.1by1, 0.05)
# identify a random subset of 5000 tests
random.idx = sample(1:nrow(data$Bhat),2000)
data.random <- mash_set_data(Bhat[random.idx, ], Shat[random.idx, ])
data.strong <- mash_set_data(Bhat[strong.idx, ], Shat[strong.idx, ])

U.c <- cov_canonical(data.random)
U.pca <- cov_pca(data.strong, 5)
U.ed <- cov_ed(data.strong, U.pca)
Ulist <- c(U.c, U.ed)
m.fit <- mash(data.random, Ulist = Ulist, outputlevel = 1)
m.final <- mash(mash_set_data(Bhat, Shat), Ulist = m.fit$fitted_g$Ulist, outputlevel = 2)
#m.fit <- mash(data, Ulist = Ulist, outputlevel = 1)
#m.final <- mash(data, g = m.fit$fitted_g, fixg = TRUE, outputlevel = 2)
#save.image("SLE_mashr_fulldb_test2.RData")
save.image("SLE_mashr.RData")
all_conditions <- colnames(get_lfsr(m.final))
shared_matrix <- get_pairwise_sharing(m.final, lfsr = 0.05)
conditions_df <- data.frame(
  Condition = all_conditions,
  Tissue = gsub("_.*", "", all_conditions),
  CellType = gsub(".*_", "", all_conditions)
)

pdf("SLE_clustered_heatmap_filtered_atleast1sig.pdf", width = 12, height = 10)
Heatmap(shared_matrix, name = "Shared %",
        col = colorRamp2(breaks = c(0, 0.5, 1.0), colors = c("#89c760", "#ffffff", "#ae2f22")),
        cluster_rows = TRUE, cluster_columns = TRUE,
        clustering_distance_rows = function(x) as.dist(1 - as.matrix(x)),
        clustering_distance_columns = function(x) as.dist(1 - as.matrix(x)),
        clustering_method_rows = "complete",
        clustering_method_columns = "complete",
        show_row_dend = TRUE, show_column_dend = TRUE,
        row_dend_width = unit(3, "cm"), column_dend_height = unit(3, "cm"),
        row_names_gp = gpar(fontsize = 11),
        column_names_gp = gpar(fontsize = 11),
        column_names_rot = 45,
        rect_gp = gpar(col = "gray90", lwd = 0.3),
        border = FALSE)
dev.off()
