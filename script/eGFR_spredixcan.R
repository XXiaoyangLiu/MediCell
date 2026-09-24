#!/bin/bash
#scpredixcan
source ~/.bashrc
conda activate imlatools
cd "/node1/liuxy/MediCell/trait_compare_2/eGFR/output/TWAS/"
declare -A celltypes=(
#    ["pbmc"]="CD4T CD8T NK B Mono DC"
    ["tub"]="PT LOH DCT CDICA CDICB MNP Plasmacytoid_DC NK NKT CD8T Bcell Th_cell"
    ["glom"]="Endo Podo MC NK NKT Th_cell CD8T Bcell MNP"
)
declare -A tissue_paths=(
#    ["pbmc_db"]="/node1/liuxy/MediCell/pbmc/TWAS/usefile/db_2/"
#    ["pbmc_cov"]="/node1/liuxy/MediCell/pbmc_2/TWAS/usefile/cov_EUR/"
    ["glom_db"]="/node1/liuxy/MediCell/kidney_3/glom/TWAS/usefile/db/"
    ["glom_cov"]="/node1/liuxy/MediCell/kidney_3/glom/TWAS/usefile/cov_EUR/"
    ["tub_db"]="/node1/liuxy/MediCell/kidney_3/tub/TWAS/usefile/db/"
    ["tub_cov"]="/node1/liuxy/MediCell/kidney_3/tub/TWAS/usefile/cov_EUR/"
)
disease_name="eGFR"
gwas_files="/node1/liuxy/MediCell/trait_compare/eGFR/usefile/eGFR_GWAS/GCST90474052_glom_withdb.tsv.gz" 
#gwas_files="/node1/liuxy/MediCell/trait_compare_2/eGFR/usefile/eGFR_Karczewski_aligned_with_pbmcdb.txt.gz"
output_base="/node1/liuxy/MediCell/trait_compare_2/eGFR/output/TWAS/"
for gwas_path in "${gwas_files[@]}"; do
    disease_name=$(basename $(dirname "$gwas_path"))
    disease_name="eGFR"
    for tissue in "${!celltypes[@]}"; do
        tissue_output_dir="${output_base}/${tissue}/${disease_name}"
        mkdir -p "$tissue_output_dir"
        db_path="${tissue_paths[${tissue}_db]}"
        cov_path="${tissue_paths[${tissue}_cov]}"
        
        for ct_value in ${celltypes[$tissue]}; do
            output_file="${tissue_output_dir}/${ct_value}.csv"
            log_file="${tissue_output_dir}/${ct_value}.log"
            
            if [[ "$tissue" == "pbmc" ]]; then
                db_file="${db_path}/${ct_value}_PCs_and_age_gender_models_filtered_signif.db"
            else
                db_file="${db_path}/${ct_value}_models_filtered_signif.db"
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
                --snp_column VariantID \
                --effect_allele_column ALT \
                --non_effect_allele_column REF \
                --zscore_column zscore \
                --beta_column BETA \
                --pvalue_column P.value \
                --input_pvalue_fix 0 \
                --output_file "$output_file" > "$log_file" 2>&1
        done
    done
done 
