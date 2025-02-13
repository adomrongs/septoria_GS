#!/bin/bash

# Definir la carpeta de entrada y salida
directorio_bam="../../GATK2/Results/aligned/dedup_bam2/"
carpeta_salida="data/modified_data/bam_subset"
region="NC_018210.1:43020-47020"

# Crear la carpeta de salida si no existe
mkdir -p "$carpeta_salida"

# Procesar cada archivo .bam en el directorio de entrada
for archivo in ${directorio_bam}/*.bam; do
    if [ -f "$archivo" ]; then
        nombre_base=$(basename "$archivo")
        samtools view -b "$archivo" "$region" -o "${carpeta_salida}/${nombre_base%.bam}_subset.bam"
    fi
done
