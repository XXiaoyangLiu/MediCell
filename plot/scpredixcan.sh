#!/bin/bash
#scpredixcan
source ~/.bashrc
conda activate imlatools
cd "/node1/liuxy/MediCell/method_compare_2/output/scpredixcan_TWAS/"
declare -A celltypes=(
    ["pbmc"]="CD14-low_CD16-positive_monocyte central_memory_CD4-positive_alpha-beta_T_cell effector_memory_CD8-positive_alpha-beta_T_cell mucosal_invariant_T_cell 
    plasmablast CD14-positive_monocyte central_memory_CD8-positive_alpha-beta_T_cell erythrocyte naive_B_cell plasmacytoid_dendritic_cell CD16-negative_CD56-bright_natural_killer_cell_human
    conventional_dendritic_cell gamma-delta_T_cell naive_thymus-derived_CD4-positive_alpha-beta_T_cell platelet CD4-  positive_alpha-beta_cytotoxic_T_cell dendritic_cell hematopoietic_precursor_cell
    naive_thymus-derived_CD8-positive_alpha-beta_T_cell regulatory_T_cell CD4-positive_alpha-beta_T_cell double_negative_thymocyte innate_lymphoid_cell natural_killer_cell transitional_stage_B_cell
    CD8-positive_alpha-beta_T_cell effector_memory_CD4-positive_alpha-beta_T_cell memory_B_cell peripheral_blood_mononuclear_cell")
declare -A tissue_paths=(
    ["pbmc_db"]="/node1/liuxy/add_plot/pbmc_compare/usefile/onek1k/scPredixcan/db/"
    ["pbmc_cov"]="/node1/liuxy/add_plot/pbmc_compare/usefile/onek1k/scPredixcan/cov/")

disease="SLE"
gwas_files="/node1/liuxy/MediCell/method_compare_2/usefile/GWAS_scpredixcan/SLE_metal_096_917_183_filt.tsv"
output_base="/node1/liuxy/MediCell/method_compare_2/output/scpredixcan_TWAS/"

for gwas_path in "${gwas_files[@]}"; do
    disease_name=$(basename $(dirname "$gwas_path"))4
    for tissue in "${!celltypes[@]}"; do
        tissue_output_dir="${output_base}/${tissue}/${disease_name}"
        mkdir -p "$tissue_output_dir"
        db_path="${tissue_paths[${tissue}_db]}"
        cov_path="${tissue_paths[${tissue}_cov]}"
        
        for ct_value in ${celltypes[$tissue]}; do
            output_file="scpredixcan/${ct_value}.csv"
            log_file="${tissue_output_dir}/${ct_value}.log"
            
            if [[ "$tissue" == "pbmc" ]]; then
                db_file="${db_path}/${ct_value}.db"
            else
                db_file="${db_path}/${ct_value}.db"
            fi
            cov_file="${cov_path}/${ct_value}_covariances.txt.gz"            
            if [[ ! -f "$db_file" ]]; then
                echo "DB file not found: $db_file" >&2 
                continue
            fi
            if [[ ! -f "$cov_file" ]]; then
                echo "Covariance file not found: $cov_file" >&2
                continue
            fi
            echo "Processing ${disease_name} - ${tissue} - ${ct_value}"
            /public/home/liuxy/MetaXcan-master/software/SPrediXcan.py \
                --model_db_path "$db_file" \
                --covariance "$cov_file" \
                --keep_non_rsid \
                --model_db_snp_key VarID \
                --gwas_file "$gwas_files" \
                --snp_column MarkerName \
                --effect_allele_column Allele1 \
                --non_effect_allele_column Allele2 \
                --beta_column Effect \
                --pvalue_column P-value \
                --se_column StdErr \
                --input_pvalue_fix 0 \
                --output_file "$output_file" > "$log_file" 2>&1
        done
    done
done
