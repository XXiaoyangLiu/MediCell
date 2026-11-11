#!/bin/bash
CELLS=("B" "CD4T" "CD8T" "DC" "Mono" "NK")
for ct in "${CELLS[@]}"; do
    echo "${ct}"
    qsub -N "pbmc_exp0.2_cov_${ct}" \
    -l nodes=tc6000:ppn=1 \
    -l walltime=96:00:00 \
    -o "/node1/liuxy/MediCell/pbmc_2/TWAS/usefile/cov_EUR_EAS/cov_${ct}.o" \
    -e "/node1/liuxy/MediCell/pbmc_2/TWAS/usefile/cov_EUR_EAS/cov_${ct}.e" \
    << END
#!/bin/bash
source ~/.bashrc
conda activate imlatools
python3 "/public/home/liuxy/MetaXcan-master/software/M01_covariances_correlations.py" \
--weight_db "/node1/liuxy/MediCell/pbmc/TWAS/usefile/db_2/${ct}_PCs_and_age_gender_models_unfiltered.db" \
--input_folder /node1/liuxy/MediCell/pbmc_2/TWAS/usefile/cov_EUR_EAS/ \
--delimiter \$'\\t' \
--covariance_output /node1/liuxy/MediCell/pbmc_2/TWAS/usefile/cov_EUR_EAS/${ct}.txt.gz
END
done
