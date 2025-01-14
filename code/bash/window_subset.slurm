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

# Load required modules
module load BCFtools/1.15.1-GCC-11.3.0
module load VCFtools/0.1.16-GCC-10.2.0 
module load PLINK/1.9b_6.21-x86_64
module load bzip2/1.0.8-GCCcore-12.2.0

# Variables
regions_file='data/modified_data/window.list'
ldmatrix='outputs/postGWAS_sep/ldmatrix'

# Create directories if not exist
mkdir -p "${ldmatrix}"
mkdir -p "data/raw_data"

# Function to process VCF files
process_vcf() {
  local region=$1
  
  echo "Processing region: ${region}"
  local output_plink="outputs/postGWAS_sep/${region}"
  mkdir -p "${output_plink}"
  
  # Extract region and compress directly to .vcf.gz
  bcftools view -r "${region}" data/raw_data/septoria_ale_clean.vcf.gz -O z -o data/raw_data/${region}.vcf.gz
  
  # Convert to PLINK format
  vcftools --gzvcf data/raw_data/${region}.vcf.gz --plink --out "${output_plink}/${region}"
  
  # Generate PLINK binary files
  plink --file "${output_plink}/${region}" --make-bed --out "${output_plink}/${region}"
  
  # Calculate LD matrix
  plink --bfile "${output_plink}/${region}" --r2 --matrix --out "${ldmatrix}/${region}"
}

# Read regions and process each
while IFS= read -r region; do
  process_vcf "$region"
done < "${regions_file}"