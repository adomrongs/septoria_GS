#!/bin/bash

module load BCFtools/1.15.1-GCC-11.3.0
module load VCFtools/0.1.16-GCC-10.2.0
module load PLINK/1.9b_6.21-x86_64
module load bzip2/1.0.8-GCCcore-12.2.0
module load SAMtools/1.9-GCC-8.2.0-2.31.1

if [ "$#" -lt 2 ]; then
    exit 1
fi

lista_archivos=$1
carpeta_salida=$2

mkdir -p "$carpeta_salida"

region="NC_018210.1:43020-47020"

while IFS= read -r archivo; do
    if [ -f "$archivo" ]; then
        nombre_base=$(basename "$archivo")
        
        if [[ "$archivo" == *.bam ]]; then
            samtools view -b "$archivo" "$region" -o "${carpeta_salida}/${nombre_base%.bam}_subset.bam"
            
        elif [[ "$archivo" == *.vcf || "$archivo" == *.vcf.gz" ]]; then
            bcftools view -r "$region" "$archivo" -o "${carpeta_salida}/${nombre_base%.vcf}_subset.vcf"
        fi
    fi
done < "$lista_archivos"