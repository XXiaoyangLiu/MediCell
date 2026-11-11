#!/bin/bash
source ~/.bashrc
conda activate imlatools
declare -A celltypes=(
    ["pbmc"]="B CD4T CD8T DC Mono NK"
)
disease_name="SLE"
gwas_files="/node1/liuxy/MediCell/pbmc_2/TWAS/usefile/GWAS/SLE_metal_096_917_183_filt.tsv"

declare -A tissue_paths=(
    ["pbmc_db"]="/node1/liuxy/MediCell/pbmc/TWAS/usefile/db_2/"
    ["pbmc_cov"]="/node1/liuxy/MediCell/pbmc_2/TWAS/usefile/cov_EUR_EAS/"
)
output_base="/node1/liuxy/MediCell/pbmc_2/TWAS/output/"
for tissue in "${!celltypes[@]}"; do
    tissue_output_dir="${output_base}/${tissue}/${disease_name}"
    mkdir -p "$tissue_output_dir"
    db_path="${tissue_paths[${tissue}_db]}"
    cov_path="${tissue_paths[${tissue}_cov]}"
    for ct_value in ${celltypes[$tissue]}; do
        output_file="${tissue_output_dir}/${ct_value}.csv"
        log_file="${tissue_output_dir}/${ct_value}.log"
        db_file="${db_path}/${ct_value}_PCs_and_age_gender_models_unfiltered.db"
        cov_file="${cov_path}/${ct_value}.txt.gz"
        if [[ ! -f "$db_file" ]]; then
            echo "db not exist: $db_file" >&2
            continue
        fi
        if [[ ! -f "$cov_file" ]]; then
            echo "cov not exist: $cov_file" >&2
            continue
        fi
        echo "${disease_name} | ${tissue} | ${ct_value}"
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



#glom
declare -A celltypes=(
    ["glom"]="Endo Podo MC NK NKT Th_cell CD8T Bcell MNP"
)
disease_name="SLE"
gwas_files="/node1/liuxy/check/other_SLE_GWAS/metal/20260830/SLE_metal_096_917_183_filt_3.tsv"

declare -A tissue_paths=(
    ["glom_db"]="/node1/liuxy/MediCell/kidney_3/glom/TWAS/usefile/db/"
    ["glom_cov"]="/node1/liuxy/MediCell/kidney_3/glom/TWAS/usefile/cov_EUR_EAS/"
)
output_base="/node1/liuxy/MediCell/kidney_3/glom/TWAS/output/"
for tissue in "${!celltypes[@]}"; do
    tissue_output_dir="${output_base}/${tissue}/${disease_name}"
    mkdir -p "$tissue_output_dir"
    db_path="${tissue_paths[${tissue}_db]}"
    cov_path="${tissue_paths[${tissue}_cov]}"
    for ct_value in ${celltypes[$tissue]}; do
        output_file="${tissue_output_dir}/${ct_value}.csv"
        log_file="${tissue_output_dir}/${ct_value}.log"
        db_file="${db_path}/${ct_value}_models_filtered_signif.db"
        cov_file="${cov_path}/${ct_value}.txt.gz"
        if [[ ! -f "$db_file" ]]; then
            echo "db not exist: $db_file" >&2
            continue
        fi
        if [[ ! -f "$cov_file" ]]; then
            echo "cov not exist: $cov_file" >&2
            continue
        fi
        echo "${disease_name} | ${tissue} | ${ct_value}"
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

#tub
#!/bin/bash
source ~/.bashrc
conda activate imlatools
declare -A celltypes=(
    ["tub"]="Bcell CD8T CDICA CDICB DCT LOH MNP NK NKT PT Th_cell Plasmacytoid_DC"
)
disease_name="SLE"
gwas_files="/node1/liuxy/check/other_SLE_GWAS/metal/20260830/SLE_metal_096_917_183_filt_3.tsv"

declare -A tissue_paths=(
    ["tub_db"]="/node1/liuxy/MediCell/kidney_3/tub/TWAS/usefile/db/"
    ["tub_cov"]="/node1/liuxy/MediCell/kidney_3/tub/TWAS/usefile/cov_EUR_EAS/"
)
output_base="/node1/liuxy/MediCell/kidney_3/tub/TWAS/output/"
for tissue in "${!celltypes[@]}"; do
    tissue_output_dir="${output_base}/${tissue}/${disease_name}"
    mkdir -p "$tissue_output_dir"
    db_path="${tissue_paths[${tissue}_db]}"
    cov_path="${tissue_paths[${tissue}_cov]}"
    for ct_value in ${celltypes[$tissue]}; do
        output_file="${tissue_output_dir}/${ct_value}.csv"
        log_file="${tissue_output_dir}/${ct_value}.log"
        db_file="${db_path}/${ct_value}_models_filtered_signif.db"
        cov_file="${cov_path}/${ct_value}.txt.gz"
        if [[ ! -f "$db_file" ]]; then
            echo "db not exist: $db_file" >&2
            continue
        fi
        if [[ ! -f "$cov_file" ]]; then
            echo "cov not exist: $cov_file" >&2
            continue
        fi
        echo "${disease_name} | ${tissue} | ${ct_value}"
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

