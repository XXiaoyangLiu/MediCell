# note: all gene for db#####
# 2025/06/06 zhaoq##########

## PBMC
### 1. genotype processing######

#### MAF0.05
```shell
#!/bin/bash
#PBS -N vcf_filter
#PBS -l nodes=1:ppn=4
#PBS -l walltime=02:00:00
#PBS -l mem=16gb
#PBS -j oe
#PBS -o /node1/zhaoq/Science2022_Data/job_output.log

input_dir="/public/home/Shenglab/DataShare/ToWufan/sc_eQTL_science2022/filter_vcf_r08"

output_dir="/node1/zhaoq/Science2022_Data/filter_vcf_r08_MAF0.05"
mkdir -p "$output_dir"

for vcf_file in "$input_dir"/*.vcf.gz; do
base_name=$(basename "$vcf_file" .vcf.gz)

output_file="$output_dir/${base_name}_MAF0.05.vcf.gz"

#/public/home/Shenglab/Softwares/vcftools/bin/bcftools view -i 'MAF>0.05' "$vcf_file" -Oz -o "$output_file"  

snp_count=$(/public/home/Shenglab/Softwares/vcftools/bin/bcftools view -H "$output_file" | wc -l)
echo -e "${base_name}_MAF0.05.vcf.gz\t$snp_count" >> "$count_file"

/public/home/Shenglab/Softwares/vcftools/bin/bcftools index "$output_file"
echo "Filtered $vcf_file and saved to $output_file"
done
cd /node1/zhaoq/scTWAS_PBMC/predictDB/genotype_preprocess/merge_chr
for chrom in {1..22}; do
/public/home/Shenglab/Softwares/vcftools/bin/bcftools query -l /public/home/Shenglab/DataShare/ToWufan/sc_eQTL_science2022/filter_vcf_r08_MAF0.05/chr${chrom}.dose.filtered.R2_0.8_MAF0.05.vcf.gz > samples_chr${chrom}.txt
done

for chrom in {2..22}; do
echo "Comparing chr1 and chr${chrom}..."
diff samples_chr1.txt samples_chr${chrom}.txt
done

ls -1 /public/home/Shenglab/DataShare/ToWufan/sc_eQTL_science2022/filter_vcf_r08_MAF0.05/*.vcf.gz > vcf_files.txt

/public/home/Shenglab/Softwares/vcftools/bin/bcftools  concat -f vcf_files.txt -o merged.vcf.gz -O z
/public/home/Shenglab/Softwares/vcftools/bin/bcftools index merged.vcf.gz

cd /public/home/Shenglab/DataShare/Reference/db_snp/hg19

/public/home/zhaoq/java/jdk-21.0.6/bin/java -jar /public/home/Shenglab/Softwares/snpEff/SnpSift.jar annotate /public/home/Shenglab/DataShare/Reference/db_snp/hg19/00-All.vcf.gz /node1/zhaoq/scTWAS_PBMC/predictDB/genotype_preprocess/merge_chr/raw_merge_chr/merged.vcf.gz > /node1/zhaoq/scTWAS_PBMC/predictDB/genotype_preprocess/merge_chr/id_tran/merged.dbsnp.vcf
bgzip -f merged.dbsnp.vcf

#从rsID注释前改
"/public/home/liuxy/bcftools-1.18/bcftools" index "/node1/liuxy/chplot/kidney_new/usefile/genotype/merged.vcf.gz"
"/public/home/liuxy/bcftools-1.18/bcftools" query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%MAF\n' "/node1/liuxy/chplot/kidney_new/usefile/genotype/merged.vcf.gz" > /node1/liuxy/chplot/kidney_new/usefile/genotype/onek1k_maf.txt

cd /node1/zhaoq/scTWAS_PBMC/predictDB/genotype_preprocess/merge_chr/1vcf_process/merged_vcf_rsid
gunzip -c  "/node1/liuxy/chplot/kidney_new/usefile/genotype/merged.vcf.gz" | tail -n +91 | gzip > "/node1/liuxy/chplot/kidney_new/usefile/genotype/step1.vcf.gz" 

zcat step1.vcf.gz | awk -F'\t' '  
    /^#/ {  
        print $0  
        next  
    }  
    {  
        sub(/.*;/, "", $3)  
        print $0  
    }  
' | gzip > step1.1.vcf.gz   

gunzip -c step1.1.vcf.gz | awk '{ printf "%s_%s_%s_%s",$1,$2,$4,$5; for (i=10; i<=NF; i++) printf "\t%s",$i; printf "\n" }'  - | gzip > step2.vcf.gz 
gunzip -c step1.1.vcf.gz | awk '{ printf "%s_%s_%s_%s\t%s\n", $1, $2, $4, $5, $3 }' | gzip > varID.rsid.vcf.gz ##output for match varID and rsid

##注释的genotype自身的dosage，是根据概率计算（非012）
gunzip -c step2.vcf.gz | awk 'BEGIN {FS="\t"; OFS="\t"} {printf "%s",$1; for (i=2; i<=NF; i++) {split($i, arr, ":"); printf "\t%s",arr[1]} printf "\n"}' - | gzip > step3.vcf.gz
gunzip -c step3.vcf.gz | awk 'NR==1 {print $0; next} {printf "%s",$1; for(i=2; i<=NF; i++) {split($i, arr, /[\/|]/); sum=0; for(j=1; j in arr; j++) {if(arr[j] ~ /[0-9]+/) sum+=arr[j];else sum+=0} printf "\t%s",sum} printf "\n"}'  - | gzip > step4.vcf.gz
gunzip -c step4.vcf.gz | sed 's/#CHROM_POS_REF_ALT/varID/g' -> step5.vcf  # done with genotype input file
