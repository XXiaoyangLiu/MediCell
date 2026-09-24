#!/bin/bash
#scpredixcan
source ~/.bashrc
conda activate imlatools
cd "/node1/liuxy/MediCell/method_compare_2/output/pseudobulk_TWAS/"
declare -A celltypes=(
    ["pbmc"]="BimmNaive Bmem CD4all CD4effCM CD4TGFbStim CD8all CD8eff MonoC MonoNC NKact NKmat"
)
declare -A tissue_paths=(
    ["pbmc_db"]="/node1/liuxy/MediCell/method_compare/TWAS/output/pseudobulk_db/"
    ["pbmc_cov"]="/node1/liuxy/MediCell/method_compare/TWAS/output/pseudobulk_TWAS_2/"
)
disease_name="SLE"
gwas_files="/node1/liuxy/MediCell/pbmc_2/TWAS/usefile/GWAS/SLE_metal_096_917_183_filt.tsv"
output_base="/node1/liuxy/MediCell/method_compare_2/output/pseudobulk_TWAS/"
for gwas_path in "${gwas_files[@]}"; do
    disease_name=$(basename $(dirname "$gwas_path"))
    disease_name="SLE"
    for tissue in "${!celltypes[@]}"; do
        tissue_output_dir="${output_base}/${tissue}/${disease_name}"
        mkdir -p "$tissue_output_dir"
        db_path="${tissue_paths[${tissue}_db]}"
        cov_path="${tissue_paths[${tissue}_cov]}"
        
        for ct_value in ${celltypes[$tissue]}; do
            output_file="${tissue_output_dir}/${ct_value}.csv"
            log_file="${tissue_output_dir}/${ct_value}.log"
            
            if [[ "$tissue" == "pbmc" ]]; then
                db_file="${db_path}/${ct_value}_models_filtered_signif.db"
            else
                db_file="${db_path}/${ct_value}.db"
            fi
            cov_file="${cov_path}/${ct_value}.txt.gz"            
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
