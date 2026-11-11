#!/bin/bash
#glom
CELLS=("Bcell" "Endo" "MC" "MNP" "NK" "NKT" "Podo" "Th_cell" "CD8T")
for ct in "${CELLS[@]}"; do
    echo "${ct}"
    qsub -N "glom_exp0.2_cov_${ct}" \
    -l nodes=tc6000:ppn=1 \
    -l walltime=96:00:00 \
    -o "/node1/liuxy/MediCell/kidney_3/glom/TWAS/usefile/cov_EUR_EAS/cov_${ct}.o" \
    -e "/node1/liuxy/MediCell/kidney_3/glom/TWAS/usefile/cov_EUR_EAS/cov_${ct}.e" \
    << END
#!/bin/bash
source ~/.bashrc
conda activate imlatools
python3 "/public/home/liuxy/MetaXcan-master/software/M01_covariances_correlations.py" \
--weight_db "/node1/liuxy/MediCell/kidney_3/glom/TWAS/usefile/db/${ct}_models_unfiltered.db" \
--input_folder /node1/liuxy/MediCell/kidney_20260811/TWAS/usefile/cov_EUR_EAS/ \
--delimiter \$'\\t' \
--covariance_output /node1/liuxy/MediCell/kidney_3/glom/TWAS/usefile/cov_EUR_EAS/${ct}.txt.gz
END
done

#tub
#!/bin/bash
CELLS=("Bcell" "CD8T" "CDICA" "CDICB" "DCT" "LOH" "MNP" "NK" "NKT" "PT" "Th_cell" "Plasmacytoid_DC")
for ct in "${CELLS[@]}"; do
    echo "${ct}"
    qsub -N "tub_exp0.2_cov_${ct}" \
    -l nodes=tc6000:ppn=1 \
    -l walltime=96:00:00 \
    -o "/node1/liuxy/MediCell/kidney_3/tub/TWAS/usefile/cov_EUR_EAS/cov_${ct}.o" \
    -e "/node1/liuxy/MediCell/kidney_3/tub/TWAS/usefile/cov_EUR_EAS/cov_${ct}.e" \
    << END
#!/bin/bash
source ~/.bashrc
conda activate imlatools
python3 "/public/home/liuxy/MetaXcan-master/software/M01_covariances_correlations.py" \
--weight_db "/node1/liuxy/MediCell/kidney_3/tub/TWAS/usefile/db/${ct}_models_unfiltered.db" \
--input_folder /node1/liuxy/MediCell/kidney_20260811/TWAS/usefile/cov_EUR_EAS/ \
--delimiter \$'\\t' \
--covariance_output /node1/liuxy/MediCell/kidney_3/tub/TWAS/usefile/cov_EUR_EAS/${ct}.txt.gz
END
done
