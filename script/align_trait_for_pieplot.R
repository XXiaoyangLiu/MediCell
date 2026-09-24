library(data.table)
trait_order_full <- c("SLE", "eGFR", "UA", "uAlb", "uCr", "uK", "uNa", "SBP")
gene_pos <- fread("/node1/liuxy/index/gene_pos_grch37.82.csv")
mask_mhc <- gene_pos[chr == "6" & start >= 28477797 & end <= 33448354]
mhc_genes <- mask_mhc$gene_name
priority_rank <- c("SLE"=0, "eGFR"=1, "UA"=2, "uAlb"=3, "uCr"=4, "uK"=5, "uNa"=6, "SBP"=7)
MIN_SHARED <- 4
MAX_GENE_DISPLAY <- 50

dataset_cfg <- list(
  PBMC = list(
    name = "PBMC",
    cell_order = c("CD4T", "CD8T", "NK", "B", "Mono", "DC"),
    limit_gene = TRUE,
    files = c(
      "/node1/liuxy/MediCell/pbmc_2/TWAS/output/pbmc/SLE/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/eGFR/output/TWAS/pbmc/eGFR/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/UA/pbmc/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uCr/pbmc/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uK/pbmc/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uNa/pbmc/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/SBP/pbmc/summary/allcelltype_result.csv"
    )
  ),
  Glom = list(
    name = "Glom",
    cell_order = c("Endo", "Podo", "MC", "MNP", "NK", "NKT", "CD8T", "Th_cell", "Bcell"),
    limit_gene = TRUE,
    files = c(
      "/node1/liuxy/MediCell/kidney_3/glom/TWAS/output/glom/SLE/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/eGFR/output/TWAS/glom/eGFR/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/UA/glom/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uAlb/glom/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uCr/glom/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uK/glom/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uNa/glom/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/SBP/glom/summary/allcelltype_result.csv"
    )
  ),
  Tub = list(
    name = "Tub",
    cell_order = c("PT","LOH","DCT","CDICA","CDICB","MNP","Plasmacytoid_DC","NK","NKT","CD8T","Bcell","Th_cell"),
    limit_gene = FALSE,
    files = c(
      "/node1/liuxy/MediCell/kidney_3/tub/TWAS/output/tub/SLE/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/eGFR/output/TWAS/tub/eGFR/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/UA/tub/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uAlb/tub/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uCr/tub/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uK/tub/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uNa/tub/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/SBP/tub/summary/allcelltype_result.csv"
    )
  )
)

dt_plot_wide <- list()
dt_meta <- list()

for(ds in dataset_cfg){
  name <- ds$name
  cell_order <- ds$cell_order
  limit_gene <- ds$limit_gene
  file_list <- ds$files

  df_list <- lapply(file_list, fread)
  all_dt <- rbindlist(df_list, fill=TRUE)
  all_dt <- all_dt[!gene_name %in% mhc_genes]
  sig_dt <- all_dt[FDR.sig == TRUE, .(gene_name, cell_type, gwas_id)]
  sig_dt <- unique(sig_dt)

  present_traits <- intersect(trait_order_full, unique(sig_dt$gwas_id))
  gene_cell_cnt <- sig_dt[, .(cell_trait_n = uniqueN(gwas_id)), by=.(gene_name, cell_type)]
  gene_max_cell <- gene_cell_cnt[, .(max_cell_trait = max(cell_trait_n)), by=gene_name]
  keep_genes <- gene_max_cell[max_cell_trait >= MIN_SHARED, gene_name]
  sig_filtered <- sig_dt[gene_name %in% keep_genes]
  sig_filtered <- merge(sig_filtered, gene_max_cell, by = "gene_name")
  trait_rank <- setNames(seq_along(trait_order_full), trait_order_full)
  sig_filtered[, priority := trait_rank[gwas_id]]

  gene_agg <- sig_filtered[, .(
    max_cell_trait = first(max_cell_trait),
    min_priority   = min(priority)
  ), by = gene_name]
  gene_agg <- gene_agg[order(-max_cell_trait, min_priority)]
  gene_list_all <- gene_agg$gene_name

  if(limit_gene){
    gene_list <- gene_list_all[1:MAX_GENE_DISPLAY]
    sig_filtered <- sig_filtered[gene_name %in% gene_list]
  }else{
    gene_list <- gene_list_all
  }
  wide_dt <- dcast(sig_filtered, gene_name + cell_type ~ gwas_id, fun.aggregate = length, fill=0)
  for(tr in present_traits){
    if(!tr %in% colnames(wide_dt)){
      wide_dt[, (tr) := 0L]
    }
  }
  keep_cols <- c("gene_name","cell_type", present_traits)
  wide_dt <- wide_dt[, ..keep_cols]
  gene_x <- data.table(gene_name = gene_list, x = seq_along(gene_list))
  cell_y <- data.table(cell_type = cell_order, y = seq_along(cell_order))
  wide_dt <- merge(wide_dt, gene_x, by="gene_name")
  wide_dt <- merge(wide_dt, cell_y, by="cell_type")
  wide_dt[, dataset := name]

  dt_plot_wide[[name]] <- wide_dt

  meta_row <- data.table(
    name = name,
    cell_order = paste(cell_order, collapse = ","),
    gene_list = paste(gene_list, collapse = ","),
    all_traits = paste(present_traits, collapse=",")
  )
  dt_meta[[name]] <- meta_row
}
plot_wide_all <- rbindlist(dt_plot_wide, fill = TRUE)
plot_wide_all[is.na(plot_wide_all)] <- 0L
meta_all <- rbindlist(dt_meta)

fwrite(plot_wide_all, "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/plot/plot_wide.csv", sep=",", quote=FALSE)
fwrite(meta_all, "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/plot/meta_info.csv", sep=",", quote=TRUE)
