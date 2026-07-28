### TWAS summary statistic
## R=r4.4.3

library(stringr)
library(dplyr)

# gene annotation for symbol transcript
gene_grch37 <- read.csv("/node1/liuxy/index/gene_pos_grch37.82.csv")
tissue_dirs <- c(pbmc = "/node1/liuxy/MediCell/pbmc/TWAS/output/TWAS/pbmc/SLE/")
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
      tissue_name == "pbmc" ~ "peer2",
      tissue_name == "tub" ~ "peer55",
      tissue_name == "glom" ~ "peer55"
    )
    df_list <- list()
    for (file in csv_files) {
      tryCatch({
        dat <- read.csv(file)  
        filename <- basename(file) 
        cell_type <- tools::file_path_sans_ext(filename)   
        ngeneindb <- nrow(dat)  
        bonferroni_p <- p.adjust(dat$pvalue, method = "bonferroni")
        dat$fdr_p <- p.adjust(dat$pvalue, method = "BH")
        dat$cell_type <- cell_type  
        dat$ngeneindb <- ngeneindb  
        dat$bonferroni_p <- bonferroni_p  
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
      
      bonf_sig_df <- subset(final_df, bonferroni_p < 0.05)
      fdr_sig_df <- subset(final_df, fdr_p < 0.05)
      all_sig_df <- unique(rbind(bonf_sig_df, fdr_sig_df))
      final_df$FDR.sig <- ifelse(final_df$fdr_p < 0.05, 'TRUE', 'FALSE')
      final_df$bonferroni_p.sig <- ifelse(final_df$bonferroni_p < 0.05, 'TRUE', 'FALSE')
      cat("Statisitc"，tissue_name, ":\n")
      cat("total gene", nrow(final_df), "\n")
      cat("FDR significant genes", nrow(fdr_sig_df), "\n")
      cat("Bonferroni significant genes", nrow(bonf_sig_df), "\n")

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


