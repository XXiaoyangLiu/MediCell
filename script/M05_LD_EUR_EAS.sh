#!/bin/bash
BCF_PLUGIN="/public/home/liuxy/bcftools-1.18/plugins"
BCF_TOOL="/public/home/liuxy/bcftools-1.18/bcftools"
SAMPLE_LIST="/node1/liuxy/index/1000G/EUR_EAS_samples.txt"
OUT_DIR="/node1/liuxy/index/1000G_2/1000G_EUR_EAS"
mkdir -p ${OUT_DIR}
for chr in {1..22}; do
    INPUT="/node1/liuxy/index/1000G_2/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    OUTPUT="${OUT_DIR}/chr${chr}_EUR_EAS.txt.gz"

    qsub -N "1000G_EUR_EAS_chr${chr}" \
    -l nodes=node1:ppn=1 \
    -l walltime=240:00:00 \
    -o "${OUT_DIR}/chr${chr}_EUR_EAS.o" \
    -e "${OUT_DIR}/chr${chr}_EUR_EAS.e" \
<< END
#!/bin/bash
BCF_PLUGIN="${BCF_PLUGIN}"
BCF_TOOL="${BCF_TOOL}"
SAMPLE_LIST="${SAMPLE_LIST}"
INPUT="${INPUT}"
OUTPUT="${OUTPUT}"

mkdir -p "\$(dirname "\$OUTPUT")"

BCFTOOLS_PLUGINS=\$BCF_PLUGIN \$BCF_TOOL +fill-tags "\$INPUT" -Ou -- -t AF | \
\$BCF_TOOL view -S \$SAMPLE_LIST --force-samples -Ou | \
\$BCF_TOOL query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF[\t%GT]\n' | awk '
{
    chr = \$1
    pos = \$2
    ref = \$3
    alt = \$4
    af_val = \$5
    id = chr "_" pos "_" ref "_" alt

    printf "%s\t%s\t%s\t%s\t%s\t%s", chr, id, pos, ref, alt, af_val

    for (i=6; i<=NF; i++) {
        gt = \$i
        if (gt == "." || gt == ".|.") {
            printf "\tNA"
        } else {
            split(gt, a, /[|\/]/)
            dos = (a[1] != "0") + (a[2] != "0")
            printf "\t%d", dos
        }
    }
    printf "\n"
}
' | gzip > "\$OUTPUT"

echo "Complete: \$OUTPUT"
END

echo "Submitted chr${chr}"
done

#pbmc


#kidney
#!/bin/bash
for chr in {1..22}; do
qsub -N "final_fix_chr${chr}" \
-v chr="${chr}" \
-l nodes=tc6000:ppn=1,walltime=48:00:00 \
-o "/node1/liuxy/MediCell/kidney_20260811/TWAS/usefile/cov_EUR_EAS/fix_${chr}.o" \
-e "/node1/liuxy/MediCell/kidney_20260811/TWAS/usefile/cov_EUR_EAS/fix_${chr}.e" \
<< 'SCRIPTEND'
#!/bin/bash
ANNOT="/node1/liuxy/test_peer_glom/exp0.2/predictDB/usefile/new_gt/snp_annot/snp_annot.chr${chr}.txt"
IN_GZ="/node1/liuxy/index/1000G_2/1000G_EUR_EAS/chr${chr}_EUR_EAS.txt.gz"
OUT_GZ="/node1/liuxy/MediCell/kidney_20260811/TWAS/usefile/cov_EUR_EAS/chr${chr}.txt.gz"

zcat "$IN_GZ" | awk '
BEGIN{
    FS = "\t"; OFS = "\t"
}
NR == FNR {   
    if(NR == 1) next
    c = $1; p = $2;
    std_r = $4; std_a = $5;
    std_id = c "_" p "_" std_r "_" std_a
    key = c "," p
    ref_map[key] = std_r
    alt_map[key] = std_a
    id_map[key] = std_id
    next
}
{
    c = $1; pos = $3;
    g_r = $4; g_a = $5;
    key = c "," pos

    if( !(key in ref_map) ) next

    std_r = ref_map[key]
    std_a = alt_map[key]
    new_id = id_map[key]

    if (g_r == std_r && g_a == std_a) {
        final_ref = g_r
        final_alt = g_a
        final_af = $6
        dos_line_start = 7
    }
    else if (g_r == std_a && g_a == std_r) {
        final_ref = std_r
        final_alt = std_a
        final_af = 1.0 - $6

        for(k=7; k<=NF; k++){
            if($k != "NA"){
                $k = 2 - $k
            }
        }
        dos_line_start = 7
    }
    else {
        next
    }

    $2 = new_id
    $4 = final_ref
    $5 = final_alt
    $6 = final_af

    print $0
}
' "$ANNOT" <(zcat "$IN_GZ") | gzip -c > "$OUT_GZ"

echo "chr${chr} correction finished, output: $OUT_GZ"
SCRIPTEND

echo "Submitted chromosome ${chr}"
sleep 1
done 
