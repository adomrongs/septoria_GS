#!/bin/bash

module load BCFtools/1.15.1-GCC-11.3.0
module load VCFtools/0.1.16-GCC-10.2.0
module load PLINK/1.9b_6.21-x86_64
module load bzip2/1.0.8-GCCcore-12.2.0

# Definir archivos y directorios
INPUT_VCF="data/raw_data/septoria_ale_clean.vcf.gz"  # Archivo VCF original
LIST_DIR="outputs/lists"                              # Carpeta con las listas
OUTPUT_DIR="outputs/ld"                    # Carpeta de salida

# Crear directorio de salida si no existe
mkdir -p $OUTPUT_DIR

# Iterar sobre cada lista en LIST_DIR, excluyendo samples.txt
for LIST_FILE in $LIST_DIR/*.list; do
    # Excluir samples.txt
    if [[ $(basename "$LIST_FILE") == "samples.txt" ]]; then
        echo "Ignorando $LIST_FILE"
        continue
    fi

    REGION=$(basename "$LIST_FILE" .list)  # Extraer nombre de la región

    echo "Procesando región: $REGION"

    # 1. Filtrar el VCF para incluir solo las muestras de la lista
    FILTERED_VCF="$OUTPUT_DIR/${REGION}.vcf.gz"
    bcftools view -S $LIST_FILE -Oz -o $FILTERED_VCF $INPUT_VCF

    echo "Filtrado completado para $REGION."

    # 2. Convertir VCF filtrado a formato PLINK
    PLINK_PREFIX="$OUTPUT_DIR/${REGION}"
    plink --vcf $FILTERED_VCF --make-bed --out $PLINK_PREFIX --double-id --allow-extra-chr

    echo "Conversión a PLINK completada para $REGION."

    # 3. Calcular LD dentro de 20 kb
    LD_RESULTS="$OUTPUT_DIR/${REGION}"
    plink --bfile $PLINK_PREFIX --ld-window-kb 20 --ld-window-r2 0 --r2 --out $LD_RESULTS --allow-extra-chr

    echo "LD calculado para $REGION. Resultados en $LD_RESULTS"
done

echo "Proceso completado para todas las regiones, excluyendo samples.txt."