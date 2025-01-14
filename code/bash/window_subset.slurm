#!/bin/bash
#SBATCH --job-name=subset
#SBATCH --output=subset_%j.log
#SBATCH --error=subset_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=05:00:00
#SBATCH --mem=40G


module load BCFtools/1.15.1-GCC-11.3.0
module load VCFtools/0.1.16-GCC-10.2.0
module load PLINK/1.9b_6.21-x86_64
module load bzip2/1.0.8-GCCcore-12.2.0

regions_file='data/modified_data/window.list'
ldmatrix='outputs/postGWAS_sep/ldmatrix'
mkdir ${ldmatrix}

process_vcf() {
  local region=$1
  
  local output_plink="outputs/postGWAS_sep/${region}"
  mkdir -p "${output_plink}"
  
  bcftools view -r "${region}" data/raw_data/septoria_ale_clean.vcf.gz -o data/raw_data/${region}.vcf
  
  bgzip -c data/raw_data/${region}.vcf > data/raw_data/${region}.vcf.gz
  
  vcftools --gzvcf data/raw_data/${region}.vcf.gz --plink --out "${output_plink}/${region}"
  
  plink --file "${output_plink}/${region}" --make-bed --out "${output_plink}/${region}"
  
  plink --file "${output_plink}/${region}" --r2 --matrix --out "${ldmatrix}/${region}"
}

while IFS= read -r region; do
echo "Procesando región: ${region}"
process_vcf "$region"
done < "$regions_file"