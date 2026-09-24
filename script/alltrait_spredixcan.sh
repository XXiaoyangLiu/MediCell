 #!/bin/bash
conda activate imlatools
cd "/node1/liuxy/MediCell/trait_compare/metabolic/output_2/"
declare -A celltypes=(
#    ["pbmc"]="CD4T CD8T NK B Mono DC"
    ["glom"]="Endo Podo MC NK NKT Th_cell CD8T Bcell MNP"
    ["tub"]="PT LOH DCT CDICA CDICB MNP Plasmacytoid_DC NK NKT CD8T Bcell Th_cell"
)
declare -A tissue_paths=(
#    ["pbmc_db"]="/node1/liuxy/MediCell/pbmc/TWAS/usefile/db_2/"
#    ["pbmc_cov"]="/node1/liuxy/MediCell/pbmc_2/TWAS/usefile/cov_EUR/"
    ["glom_db"]="/node1/liuxy/MediCell/kidney_3/glom/TWAS/usefile/db/"
    ["glom_cov"]="/node1/liuxy/MediCell/kidney_3/glom/TWAS/usefile/cov_EUR/"
    ["tub_db"]="/node1/liuxy/MediCell/kidney_3/tub/TWAS/usefile/db/"
    ["tub_cov"]="/node1/liuxy/MediCell/kidney_3/tub/TWAS/usefile/cov_EUR/"
)
declare -A gwas_map=(
    ["SBP"]="/node1/liuxy/MediCell/trait_compare_2/other_trait/usefile/SBP_aligned_with_db_kidney.txt.gz"
    ["Glu"]="/node1/liuxy/MediCell/trait_compare_2/other_trait/usefile/Glu_aligned_with_db_kidney.txt.gz"
    ["UA"]="/node1/liuxy/MediCell/trait_compare_2/other_trait/usefile/UA_aligned_with_db_kidney.txt.gz"
    ["uAlb"]="/node1/liuxy/MediCell/trait_compare_2/other_trait/usefile/uAlb_aligned_with_db_kidney.txt.gz"
    ["uCr"]="/node1/liuxy/MediCell/trait_compare_2/other_trait/usefile/uCr_aligned_with_db_kidney.txt.gz"
    ["uK"]="/node1/liuxy/MediCell/trait_compare_2/other_trait/usefile/uK_aligned_with_db_kidney.txt.gz"
    ["uNa"]="/node1/liuxy/MediCell/trait_compare_2/other_trait/usefile/uNa_aligned_with_db_kidney.txt.gz"
)
output_root="/node1/liuxy/MediCell/trait_compare_2/other_trait/output/"

for disease_name in "${!gwas_map[@]}"; do
    echo -e "$disease_name"
    gwas_path="${gwas_map[$disease_name]}"
    for tissue in "${!celltypes[@]}"; do
        echo -e "$tissue"
        tissue_output_dir="${output_root}/${disease_name}/${tissue}/"
        mkdir -p "$tissue_output_dir"
        db_path="${tissue_paths[${tissue}_db]}"
        cov_path="${tissue_paths[${tissue}_cov]}"
        for ct_value in ${celltypes[$tissue]}; do
            echo -e "$ct_value"
            output_file="${tissue_output_dir}/${ct_value}.csv"
            log_file="${tissue_output_dir}/${ct_value}.log"

            db_file="${db_path}/${ct_value}_models_filtered_signif.db"
            cov_file="${cov_path}/${ct_value}.txt.gz"
            if [[ ! -f "$db_file" ]]; then
                echo "DB file not found: $db_file" >&2
                continue
            fi
            if [[ ! -f "$cov_file" ]]; then
                echo "Covariance file not found: $cov_file" >&2
                continue
            fi
            /public/home/liuxy/MetaXcan-master/software/SPrediXcan.py \
                --model_db_path "$db_file" \
                --covariance "$cov_file" \
                --keep_non_rsid \
                --model_db_snp_key VarID \
                --gwas_file "$gwas_path" \
                --snp_column RSID \
                --effect_allele_column alt \
                --non_effect_allele_column ref \
                --zscore_column zscore \
                --beta_column beta \
                --pvalue_column P.value \
                --input_pvalue_fix 0 \
                --output_file "$output_file" > "$log_file" 2>&1
        done
    done
done
