------------------------2.2 make database ---------------------------------
#R4.1.3
#直接在R交互式运行，很快
library(RSQLite)  
library(dplyr)  
library(stringr)  

"%&%" <- function(a,b) paste(a, b, sep = '')  

driver <- dbDriver('SQLite')  

celltypes <- list.files("/node1/liuxy/MediCell/pbmc/predictdb/output_2/summary/")  
celltypes <- celltypes[str_starts(celltypes, "Model_training")]  
celltypes <- sub("Model_training_(.*?)_chr.*", "\\1", celltypes)
celltypes <-unique(celltypes)  

for(ct in celltypes) {  
  # 读取第1条染色体的模型总结  
  tryCatch({  
    model_summaries <- read.table(  
      paste0('/node1/liuxy/MediCell/pbmc/predictdb/output_2/summary/Model_training_', ct, '_chr1_model_summaries.txt'),  
      header = TRUE, stringsAsFactors = FALSE, fill = T) ##为了解决报错
    original_colnames <- colnames(model_summaries)
    mat <- as.matrix(model_summaries)
    new_mat <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
    new_mat[,1] <- mat[,1]
    new_mat[,2] <- mat[,1]
    new_mat[,3:ncol(new_mat)] <- mat[,2:(ncol(mat)-1)]
    model_summaries_fixed <- as.data.frame(new_mat, stringsAsFactors = FALSE)
    colnames(model_summaries_fixed) <- original_colnames
    model_summaries <- model_summaries_fixed
    }, error = function(e) {  
        cat("Error reading model summaries for cell type:", ct, "\n")  
        cat("Error message:", e$message, "\n")  
        return(NULL)  # 返回NULL以跳过错误的细胞类型  
        })  
  # 读取第1条染色体的总结信息  
  tiss_summary <- read.table(  
    paste0('/node1/liuxy/MediCell/pbmc/predictdb/output_2/summary/Model_training_', ct, '_chr1_summary.txt'),  
    header = TRUE, stringsAsFactors = FALSE  
  ) 
  n_samples <- tiss_summary$n_samples  
  # 读取2到22条染色体的数据  
  for (i in 2:22) {  
    model_summary_file <- paste0('/node1/liuxy/MediCell/pbmc/predictdb/output_2/summary/Model_training_', ct, '_chr', i, '_model_summaries.txt')  
    model_tmp <- read.table(model_summary_file, header = TRUE, stringsAsFactors = FALSE, fill = T)
    tryCatch({
        original_colnames <- colnames(model_tmp)
        mat <- as.matrix(model_tmp)
        new_mat <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
        new_mat[,1] <- mat[,1]
        new_mat[,2] <- mat[,1]
        new_mat[,3:ncol(new_mat)] <- mat[,2:(ncol(mat)-1)]
        model_tmp_fixed <- as.data.frame(new_mat, stringsAsFactors = FALSE)
        colnames(model_tmp_fixed) <- original_colnames
        model_summaries <- rbind(model_summaries, model_tmp_fixed)
    }, error = function(e) {  
      cat("Error reading model summaries from file:", model_summary_file, "\n")  
      cat("Error message:", e$message, "\n")  
    }) 
    summary_file <- paste0('/node1/liuxy/MediCell/pbmc/predictdb/output_2/summary/Model_training_', ct, '_chr', i, '_summary.txt')  
    tiss_summary <- rbind(tiss_summary, read.table(summary_file, header = TRUE, stringsAsFactors = FALSE))  
  }  
  model_summaries <- rename(model_summaries, gene = gene_name)   # 修改列名  
  # 创建数据库连接  
  conn_name <- paste0('/node1/liuxy/MediCell/pbmc/TWAS/usefile/db_2/', ct, '_models.db')  
  conn <- dbConnect(drv = driver, conn_name)  
  # 写入模型摘要到数据库  
  tryCatch({  
    dbWriteTable(conn, 'model_summaries', model_summaries, overwrite = TRUE)  
    dbExecute(conn, "CREATE INDEX gene_model_summary ON model_summaries (gene)")  
  }, error = function(e) {  
    cat("Error writing model summaries to database for cell type:", ct, "\n")  
    cat("Error message:", e$message, "\n")  
  })  
  # 读取权重文件  
  weights <- read.table(paste0('/node1/liuxy/MediCell/pbmc/predictdb/output_2/weights/Model_training_', ct, '_chr1_weights.txt'), header = TRUE, stringsAsFactors = FALSE)  
  for (i in 2:22) {  
    weights_file <- paste0('/node1/liuxy/MediCell/pbmc/predictdb/output_2/weights/Model_training_', ct, '_chr', i, '_weights.txt')  
    tryCatch({  
      weights <- rbind(weights, read.table(weights_file, header = TRUE, stringsAsFactors = FALSE))  
    }, error = function(e) {  
      cat("Error reading weights from file:", weights_file, "\n")  
      cat("Error message:", e$message, "\n")  
    })  
  }  
  weights <- rename(weights, gene = gene_name)  # 修改列名  
  # 写入权重到数据库  
  tryCatch({  
    dbWriteTable(conn, 'weights', weights, overwrite = TRUE)  
    dbExecute(conn, "CREATE INDEX weights_rsid ON weights (rsid)")  
    dbExecute(conn, "CREATE INDEX weights_gene ON weights (gene)")  
    dbExecute(conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")  
  }, error = function(e) {  
    cat("Error writing weights to database for cell type:", ct, "\n")  
    cat("Error message:", e$message, "\n")  
  })  
  # 写入样本信息到数据库  
  sample_info <- data.frame(n_samples = n_samples, population = 'EUR')   
  tryCatch({  
    dbWriteTable(conn, 'sample_info', sample_info, overwrite = TRUE)  
  }, error = function(e) {  
    cat("Error writing sample info to database for cell type:", ct, "\n")  
    cat("Error message:", e$message, "\n")  
  })  
  # 写入构建信息到数据库  
  construction <- tiss_summary %>% select(chrom, cv_seed) %>% rename(chromosome = chrom)  
  tryCatch({  
    dbWriteTable(conn, 'construction', construction, overwrite = TRUE)  
  }, error = function(e) {  
    cat("Error writing construction to database for cell type:", ct, "\n")  
    cat("Error message:", e$message, "\n")  
  })  
  # 断开数据库连接  
  dbDisconnect(conn)  
}  

#------------------------2.3 filter database ---------------------------------
###filter db###
for(ct in celltypes) {
  # for (pc in pc_count) {
  unfiltered_db <- paste0('/node1/liuxy/MediCell/pbmc/TWAS/usefile/db_2/', ct, '_models.db')
  filtered_db <- paste0('/node1/liuxy/MediCell/pbmc/TWAS/usefile/db_2/', ct, '_models_filtered_signif.db')
  driver <- dbDriver("SQLite")
  in_conn <- dbConnect(driver, unfiltered_db)
  out_conn <- dbConnect(driver, filtered_db)
  model_summaries <- dbGetQuery(in_conn, 'select * from model_summaries where CAST(zscore_pval AS REAL) < 0.05 and CAST(rho_avg_squared AS REAL) > 0.01')
  model_summaries <- model_summaries %>% rename(pred.perf.R2 = rho_avg_squared, genename=gene_id, pred.perf.pval = zscore_pval, n.snps.in.model = n_snps_in_model) #还是改了这里，一定要有genename，所以先这样
  model_summaries$pred.perf.qval <- NA
  dbWriteTable(out_conn, 'extra', model_summaries, overwrite = TRUE)
  construction <- dbGetQuery(in_conn, 'select * from construction')
  dbWriteTable(out_conn, 'construction', construction, overwrite = TRUE)
  sample_info <- dbGetQuery(in_conn, 'select * from sample_info')
  dbWriteTable(out_conn, 'sample_info', sample_info, overwrite = TRUE)
  weights <- dbGetQuery(in_conn, 'select * from weights')
  weights <- weights %>% filter(gene %in% model_summaries$gene) %>% rename(eff_allele = alt, ref_allele = ref, weight = beta)
  dbWriteTable(out_conn, 'weights', weights, overwrite = TRUE)
  dbExecute(out_conn, "CREATE INDEX weights_rsid ON weights (rsid)")
  dbExecute(out_conn, "CREATE INDEX weights_gene ON weights (gene)")
  dbExecute(out_conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
  dbExecute(out_conn, "CREATE INDEX gene_model_summary ON extra (gene)")
  dbDisconnect(in_conn)
  dbDisconnect(out_conn)
  # }
}
for(ct in celltypes) {
  # for (pc in pc_count) {
  unfiltered_db <- paste0('/node1/liuxy/MediCell/pbmc/TWAS/usefile/db_2/', ct, '_models.db')
  filtered_db <- paste0('/node1/liuxy/MediCell/pbmc/TWAS/usefile/db_2/', ct, '_models_unfiltered.db')
  driver <- dbDriver("SQLite")
  in_conn <- dbConnect(driver, unfiltered_db)
  out_conn <- dbConnect(driver, filtered_db)
  model_summaries <- dbGetQuery(in_conn, 'select * from model_summaries')
  model_summaries <- model_summaries %>% rename(pred.perf.R2 = rho_avg_squared, genename=gene_id, pred.perf.pval = zscore_pval, n.snps.in.model = n_snps_in_model) #还是改了这里，一定要有genename，所以先这样
  model_summaries$pred.perf.qval <- NA
  dbWriteTable(out_conn, 'extra', model_summaries, overwrite = TRUE)
  construction <- dbGetQuery(in_conn, 'select * from construction')
  dbWriteTable(out_conn, 'construction', construction, overwrite = TRUE)
  sample_info <- dbGetQuery(in_conn, 'select * from sample_info')
  dbWriteTable(out_conn, 'sample_info', sample_info, overwrite = TRUE)
  weights <- dbGetQuery(in_conn, 'select * from weights')
  weights <- weights %>% filter(gene %in% model_summaries$gene) %>% rename(eff_allele = alt, ref_allele = ref, weight = beta)
  dbWriteTable(out_conn, 'weights', weights, overwrite = TRUE)
  dbExecute(out_conn, "CREATE INDEX weights_rsid ON weights (rsid)")
  dbExecute(out_conn, "CREATE INDEX weights_gene ON weights (gene)")
  dbExecute(out_conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
  dbExecute(out_conn, "CREATE INDEX gene_model_summary ON extra (gene)")
  dbDisconnect(in_conn)
  dbDisconnect(out_conn)
  # }
}
