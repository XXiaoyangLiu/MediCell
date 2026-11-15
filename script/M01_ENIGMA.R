library(ENIGMA)
script_file <- sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(script_file)), "helpers", "input_contracts.R"))

bulk_tpm_file <- Sys.getenv("MC_BULK_TPM")
ref_file <- Sys.getenv("MC_ENIGMA_REF")
cbs_frax_file <- Sys.getenv("MC_CBS_FRAX")
out_dir <- Sys.getenv("MC_OUT_DIR")
alpha_label <- Sys.getenv("MC_ALPHA", "0.1")
alpha <- mc_number("MC_ALPHA", "0.1", lower = 0)

inputs <- mc_load_enigma_inputs(bulk_tpm_file, ref_file, cbs_frax_file)
myenigma <- create_ENIGMA(bulk = inputs$bulk, ref = inputs$ref, ref_type = "aggre")
myenigma@result_cell_proportion <- inputs$fractions
myenigma.bak <- myenigma
myenigma <- ENIGMA_trace_norm(myenigma.bak, theta = NULL, R = NULL, alpha = alpha,
                              preprocess = "log", Normalize = FALSE)
# public ENIGMA API (WWXkenmo/ENIGMA): theta/R are required positionally but
# immediately overwritten inside the function by the object's own
# @result_cell_proportion / @ref slots; Normalize=FALSE corresponds to the
# original lab version's do_cpm=FALSE (no normalization of the output CSE).
save(myenigma, file = file.path(out_dir, paste0('myenigma_alpha_', alpha_label, '.RData')))
enigma.cse.test <- sce2array(myenigma, norm_output = F)
save(enigma.cse.test, file = file.path(out_dir, paste0('enigma.cse.trace_alpha_', alpha_label, '.RData')))
