#!/bin/bash
#SBATCH --job-name=SNPeff
#SBATCH --output=SNPeff_%j.log
#SBATCH --error=SNPeff_%j.err
#SBATCH --partition=fast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=40G

#----------------------------------------------------------------------
# Define the variables and load the required modules
#----------------------------------------------------------------------
dirsnpeff="data/modified_data/"
list="data/raw_data/renamechr.txt"
snpeff="/home/arondon/other/snpEff/snpEff.jar"
vcf="data/raw_data/septoria_ale_clean2ann.vcf"
out="data/raw_data/septoria_ale_clean_ann.vcf.gz"

module load Java/1.8.0_212
module load GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8 
module load picard/2.26.10-Java-15
module load BCFtools/1.9-GCC-8.2.0-2.31.1

#----------------------------------------------------------------------
# Renaming chromosomes
#----------------------------------------------------------------------

echo "Renaming chromosomes..."

bcftools annotate --rename-chrs ${list} data/raw_data/septoria_ale_clean.vcf.gz > ${vcf}
bcftools view -Oz -o ${vcf}.gz ${vcf}


echo "Process finished"

#----------------------------------------------------------------------
# SNP annotation
#----------------------------------------------------------------------

echo "Annotating filtered snps..."

java -Xmx32g -jar ${snpeff} download Zymoseptoria_tritici
java -Xmx32g -jar ${snpeff} -v Zymoseptoria_tritici ${vcf} > $out

echo "Process finished"