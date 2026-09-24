#!/bin/bash
change GWAS for scpredixcan
cd "/node1/liuxy/MediCell/method_compare_2/usefile/GWAS_scpredixcan/chain/"
DB_DIR="/node1/liuxy/add_plot/pbmc_compare/usefile/onek1k/scPredixcan/db"
OUTPUT_FILE="all_celltypes_snps.txt"
> "$OUTPUT_FILE"

for db in "$DB_DIR"/*.db
    do
    if [ -f "$db" ]
        then
            celltype=$(basename "$db" .db)
        sqlite3 "$db" "SELECT DISTINCT rsid FROM weights;" >> "$OUTPUT_FILE"
    fi
done
sort -u "$OUTPUT_FILE" > "${OUTPUT_FILE}.sorted"

awk '
BEGIN {OFS="\t"; print "chr", "pos_hg38", "ref", "alt", "snpID"}
{snpID=$0;
 split(snpID, parts, "_");
 chr=substr(parts[1], 4);
 pos=parts[2];
 ref=parts[3];
 alt=parts[4];
 print chr, pos, ref, alt, snpID;
}' all_celltypes_snps.txt | awk -F"\t" 'NR==1 {print; next}
{if ($1=="X") key="23";
 else if ($1=="Y") key="24";
 else if ($1=="M"||$1=="MT") key="25";
 else key=sprintf("%02d", $1);
 if (NR != 1) printf "%s\t%09d\t%s\n", key, $2, $0;
}' | sort -k1,1n -k2,2n | cut -f3- > snp_index_hg38_sorted.txt

#awk '{if (NR != 1) print "chr"$2"\t"$3"\t"($3+1)"\t"$1}' '/node1/liuxy/MediCell/pbmc_2/TWAS/usefile/GWAS/GCST011096_buildGRCh37.tsv' > GCST011096_for_liftover.bed
#CrossMap bed "/node1/liuxy/index/hg19ToHg38.over.chain.gz" GCST011096_for_liftover.bed  > GCST011096_hg38.bed
#awk '
#BEGIN{OFS="\t"}
#NR==FNR{
#    gsub("chr","",$6);
#    snpid=$4;
#    a[snpid]=$6"\t"$7;
#    next
#}
#{
#    gsub(/\r/,"");
#}
#FNR==1{print $0,"chr_hg38","pos_hg38";next}
#{
#    sid=$1;
#    if (sid in a){
#        print $0, a[sid]
#    }else{
#        print $0,"NA","NA"
#    }
#}
#' "/node1/liuxy/MediCell/method_compare_2/usefile/GWAS_scpredixcan/chain/GCST011096_hg38.bed" <(cat "/node1/liuxy/MediCell/pbmc_2/TWAS/usefile/GWAS/GCST011096_buildGRCh37.tsv") | \
#gzip > "/node1/liuxy/MediCell/method_compare_2/usefile/GWAS_scpredixcan/GCST011096_hg19.tsv.gz"

#zcat /node1/liuxy/MediCell/pbmc_2/TWAS/usefile/GWAS/GCST90018917_buildGRCh37.tsv.gz |awk '{if (NR != 1) print "chr"$1"\t"$2"\t"($2+1)"\t"$1"_"$2}' > GCST90018917_for_liftover.bed
#CrossMap bed "/node1/liuxy/index/hg19ToHg38.over.chain.gz" GCST90018917_for_liftover.bed > GCST90018917_hg38.bed
#awk '
#BEGIN{OFS="\t"}
#NR==FNR{
#    gsub("chr","",$6);
#    snpid=$4;
#    a[snpid]=$6"\t"$7;
#    next
#}
#{
#    gsub(/\r/,"");
#}
#FNR==1{print $0,"chr_hg38","pos_hg38";next}
#{
#    sid=$1"_"$2;
#    if (sid in a){
#        print $0, a[sid]
#    }else{
#        print $0,"NA","NA"
#    }
#}
#' "/node1/liuxy/MediCell/method_compare_2/usefile/GWAS_scpredixcan/chain/GCST90018917_hg38.bed" <(zcat "/node1/liuxy/MediCell/pbmc_2/TWAS/usefile/GWAS/GCST90018917_buildGRCh37.tsv.gz") > "/node1/liuxy/MediCell/method_compare_2/usefile/GWAS_scpredixcan/GCST90018917_hg19.tsv.gz"
