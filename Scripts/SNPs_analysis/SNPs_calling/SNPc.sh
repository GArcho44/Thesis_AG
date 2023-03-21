#!/bin/bash -i
#
#SBATCH --job-name=SNP_call
#SBATCH --output=SNPc.txt
#
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=32000
#
# This is a script for calling SNPs using the SNP caller Mutect2

# Index the reference genome
samtools faidx ref_gen.fasta

# Create a dictionary file for the reference genome
gatk CreateSequenceDictionary -R ref_gen.fasta -O ref_gen.dict

# Index the bam file
samtools index mg.reads.sorted.bam

# Call variance using Mutect2
gatk Mutect2 --af-of-alleles-not-in-resource 0.33\
 -R ref_gen.fasta\
 -I mg.reads.sorted.bam\
 -O mutect2.variants.vcf.gz\
 -L contigs.bed\

# Unzip the vcf.gz file
gzip -d mutect2.variants.vcf.gz

# Filter the raw mutect2 vcf file
Rscript SNPs_filtering.R

# Remove the temporary file "mutect2.filtered.tmp"
rm mutect2.filtered.tmp.vcf

# This is the end of script
