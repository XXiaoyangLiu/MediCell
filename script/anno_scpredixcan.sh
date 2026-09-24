#!/bin/bash
INFILE="/node1/liuxy/MediCell/pbmc_2/TWAS/usefile/GWAS/GCST90018917_buildGRCh37.tsv.gz"
OUTFILE_PREFIX="/node1/liuxy/MediCell/method_compare_2/usefile/GWAS_scpredixcan/GCST90018917"

#zcat ${INFILE} | awk '
#BEGIN{FS="\t";OFS="\t"}
#NR==1{next}
#{
#    chr=$1; pos=$2; EA=$3; OA=$4;
#    print chr, pos, pos, OA, EA, $0
#}' > ${OUTFILE_PREFIX}.avinput

#perl /public/home/Shenglab/Softwares/annovar/annotate_variation.pl \
#-filter \
#-dbtype avsnp151 \
#-build hg19 \
#${OUTFILE_PREFIX}.avinput \
#/public/home/Shenglab/Softwares/annovar/humandb/

gawk '
BEGIN{FS="\t"}
FILENAME==ARGV[1]{
    # $3 chr; $4 pos; $6 ref(OA); $7 alt(EA); $2 rsid
    key = $3 ":" $4 ":" $6 ":" $7
    rs_map[key] = $2
    next
}
FNR==1{print; next}
{
    key = $1 ":" $2 ":" $4 ":" $3
    if(key in rs_map){
        $9 = rs_map[key]
    }
    print
}
' ${OUTFILE_PREFIX}.avinput.hg19_avsnp151_dropped <(zcat ${INFILE}) | bgzip > ${OUTFILE_PREFIX}_rsid.tsv.gz
