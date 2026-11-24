#!/bin/bash
source ~/.bashrc
conda activate imlatools
cd "/node1/liuxy/MediCell/method_compare_2/usefile/cov_EUR_EAS/"
CELLS=(BimmNaive Bmem CD4all CD4effCM CD4TGFbStim CD8all CD8eff CD8unknown MonoC MonoNC NKact NKmat Plasma)
for ct in "${CELLS[@]}"; do
qsub -N "cov_${ct}" \
-l nodes=tc6000:ppn=1 \
-l mem=100G \
-l walltime=2400:00:00 \
-o "/node1/liuxy/MediCell/method_compare_2/usefile/cov_EUR_EAS/cov_${ct}.o" \
-e "/node1/liuxy/MediCell/method_compare_2/usefile/cov_EUR_EAS/cov_${ct}.e" \
<< END
#!/bin/bash
source ~/.bashrc
conda activate imlatools
python3 "/public/home/liuxy/MetaXcan-master/software/M01_covariances_correlations.py" \
--weight_db "/node1/liuxy/MediCell/method_compare/TWAS/output/pseudobulk_db/${ct}_models_filtered_signif.db" \
--input_folder /node1/liuxy/MediCell/kidney_20260811/TWAS/usefile/cov_EUR_EAS/ \
--delimiter \$'\t' \
--covariance_output /node1/liuxy/MediCell/method_compare_2/usefile/cov_EUR_EAS/${ct}.txt.gz
END
done
