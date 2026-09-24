library(data.table)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)

gene_pos_df <- fread("/node1/liuxy/index/gene_pos_grch37.82.csv", data.table = FALSE)
gene_pos_df$start <- as.numeric(gene_pos_df$start)
gene_pos_df$end   <- as.numeric(gene_pos_df$end)

mhc_genes <- gene_pos_df %>%
  filter(chr == "6", start >= 28477797, end <= 33448354) %>%
  pull(gene_name) %>%
  unique()

tissue_configs <- list(
  pbmc = list(
    name = "pbmc",
    files = c(
      "/node1/liuxy/MediCell/pbmc_2/TWAS/output/pbmc/SLE/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/eGFR/output/TWAS/pbmc/eGFR/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/UA/pbmc/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uAlb/pbmc/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uCr/pbmc/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uK/pbmc/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/uNa/pbmc/summary/allcelltype_result.csv",
      "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/SBP/pbmc/summary/allcelltype_result.csv"
    )
  ),
  glom = list(
    name = "glom",
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
  tub = list(
    name = "tub",
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

priority_traits <- c("SLE", "eGFR", "UA", "uAlb", "uCr", "uK", "uNa", "SBP")
process_tissue <- function(config) {
  tissue_name <- config$name
  files <- config$files
  df_list <- lapply(files, function(fp) {
    tmp <- fread(fp, data.table = FALSE)
    if ("FDR.sig" %in% colnames(tmp)) {
      tmp$FDR.sig <- as.logical(tmp$FDR.sig)
    }
    return(tmp)
  })
  all_data <- bind_rows(df_list)
  all_data <- all_data %>% filter(!gene_name %in% mhc_genes)
  mask_sig <- !is.na(all_data$FDR.sig) & all_data$FDR.sig
  final_data <- all_data[mask_sig, ]
  shared_info <- final_data %>%
    group_by(gene_id, gene_name, cell_type) %>%
    summarise(
      shared_trait_count = n_distinct(gwas_id),
      .groups = "drop"
    )
  long_stats <- final_data %>%
    select(gene_id, gene_name, cell_type, gwas_id, zscore, pvalue, fdr_p, FDR.sig) %>%
    distinct(gene_id, gene_name, cell_type, gwas_id, .keep_all = TRUE) %>%
    pivot_longer(
      cols = c(zscore, pvalue, fdr_p, FDR.sig),
      names_to = "stat_type",
      values_to = "stat_val"
    ) %>%
    mutate(col_name = paste0(gwas_id, "_", stat_type)) %>%
    select(-gwas_id, -stat_type)
  trait_wide <- long_stats %>%
    pivot_wider(
      names_from = col_name,
      values_from = stat_val,
      values_fill = NA
    )
  for (tr in priority_traits) {
    col_z <- paste0(tr, "_zscore")
    col_p <- paste0(tr, "_pvalue")
    col_fdr_val <- paste0(tr, "_fdr_p")
    col_fdr_sig <- paste0(tr, "_FDR.sig")
    if (!col_z %in% colnames(trait_wide)) trait_wide[[col_z]] <- NA
    if (!col_p %in% colnames(trait_wide)) trait_wide[[col_p]] <- NA
    if (!col_fdr_val %in% colnames(trait_wide)) trait_wide[[col_fdr_val]] <- NA
    if (!col_fdr_sig %in% colnames(trait_wide)) trait_wide[[col_fdr_sig]] <- NA
  }
  stat_col_order <- c()
  for (tr in priority_traits) {
    stat_col_order <- c(stat_col_order, 
                         paste0(tr, "_zscore"), 
                         paste0(tr, "_pvalue"), 
                         paste0(tr, "_fdr_p"),
                         paste0(tr, "_FDR.sig"))
  }
  tissue_result <- shared_info %>%
    left_join(trait_wide, by = c("gene_id", "gene_name", "cell_type")) %>%
    mutate(tissue = tissue_name) %>%
    select(
      gene_id,
      gene_name,
      tissue,
      cell_type,
      shared_trait_count,
      all_of(stat_col_order)
    )
  build_trait_str <- function(...) {
    row_list <- list(...)
    names(row_list) <- names(tissue_result)
    sig_vec <- c()
    detail_vec <- c()
    
    for(tr in priority_traits){
      fdr_sig_col <- paste0(tr, "_FDR.sig")
      z_col <- paste0(tr, "_zscore")
      p_col <- paste0(tr, "_pvalue")
      fdr_val_col <- paste0(tr, "_fdr_p")
      
      fdr_sig_raw <- row_list[[fdr_sig_col]]
      fdr_is_sig <- FALSE
      if(!is.na(fdr_sig_raw)){
        if(is.logical(fdr_sig_raw)){
          fdr_is_sig <- fdr_sig_raw
        }else{
          fdr_is_sig <- as.numeric(fdr_sig_raw) == 1
        }
      }
      
      if(fdr_is_sig){
        z_val <- sprintf("%.3f", row_list[[z_col]])
        p_val <- sprintf("%.2e", row_list[[p_col]])
        fdr_val <- sprintf("%.2e", row_list[[fdr_val_col]])
        
        sig_vec <- c(sig_vec, tr)
        detail_str <- paste0(tr, ": zscore=", z_val, ";pvalue=", p_val, ";FDR=", fdr_val)
        detail_vec <- c(detail_vec, detail_str)
      }
    }
    sig_trait <- paste(sig_vec, collapse = ";")
    detail_str_total <- paste(detail_vec, collapse = ";")
    return(tibble(sig_trait_list = sig_trait, trait_detail_str = detail_str_total))
  }
  
  str_df <- pmap_dfr(tissue_result, build_trait_str)
  tissue_result <- bind_cols(tissue_result, str_df)
  
  orig_col <- colnames(tissue_result)
  new_col_order <- c(
    orig_col[!orig_col %in% c("sig_trait_list","trait_detail_str")],
    "sig_trait_list",
    "trait_detail_str"
  )
  tissue_result <- tissue_result[, new_col_order]
  
  message(tissue_name, " done, total rows: ", nrow(tissue_result))
  return(tissue_result)
}
all_tissue_results <- lapply(tissue_configs, process_tissue)
final_merged_table <- bind_rows(all_tissue_results)

output_dir <- "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/plot/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (tiss in names(all_tissue_results)) {
  write.csv(
    all_tissue_results[[tiss]],
    file.path(output_dir, paste0("TWAS_", tiss, "_shared_detail_table=.csv")),
    row.names = FALSE,
    quote = FALSE
  )
  message("Export single tissue table: ", paste0("TWAS_", tiss, "_shared_detail_table=.csv"))
}
write.csv(
  final_merged_table,
  file.path(output_dir, "TWAS_all_tissue_merged_shared_detail_table.csv"),
  row.names = FALSE,
  quote = FALSE
)
