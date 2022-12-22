#!/bin/bash -i
#
#SBATCH --job-name=sim
#SBATCH --output=sim.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=48000
#
# Reads simulation, mapping to reference & SNP calling. VALIDADATION of the SNP callers.

# This is a script for simulating and mapping reads to reference genome & Calling SNPs with three different tools: Mutect2, BCFtools & HaplotypeCaller

# Activate the appropriate conda environment (seqtools) for using InSilicoSeq tool
conda activate seqtools

# Simulate reads using the InSilicoSeq tool
iss generate --n_reads 11M --draft contigs.bin3.fasta mutatedg4.fasta --abundance_file abundance.txt\
 --model novaseq --output simulated.reads

# De-activate conda environment
conda deactivate

# Remove the temporary files
rm simulated.reads.iss.tmp.*

# Preprocess the simulated reads with the IMP3 pipeline 
/zfs/omics/projects/metatools/TOOLS/IMP3/runIMP3 -l -t 2 config.preprocessing.yaml

# Copy the preprocessed  output to the current directory
cp IMP3_Preprocessing/Preprocessing/mg.r1.preprocessed.fq .
cp IMP3_Preprocessing/Preprocessing/mg.r2.preprocessed.fq .

# Activate the appropriate conda environment for using BWA
conda activate /zfs/omics/projects/metatools/TOOLS/IMP3/conda/76e02d85877600c41ac84cb7bc27a189

# First the reference genome need to be indexed
bwa index contigs.bin3.fasta

# Run reads mapping
bwa mem -v 1 -t 1 contigs.bin3.fasta mg.r1.preprocessed.fq mg.r2.preprocessed.fq | samtools view --threads 2 -bS - > mapped.reads.bam

# De-activate the conda environment
conda deactivate

# Index the reference genome
samtools faidx contigs.bin3.fasta

# Create a dictionary file for the reference genome
gatk CreateSequenceDictionary -R contigs.bin3.fasta -O contigs.bin3.dict

# Sort the initial bam file
samtools sort mapped.reads.bam -o mapped.reads.sorted.bam

# Index the initial bam file
samtools index mapped.reads.sorted.bam

# Add the appropriate RG
samtools addreplacerg -r '@RG\tID:mg\tSM:val.test' mapped.reads.sorted.bam -o mapped.reads.sorted1.bam

# Remove the old bam file and index
rm mapped.reads.sorted.bam
rm mapped.reads.sorted.bam.bai

# Re-name the final bam file
mv mapped.reads.sorted1.bam mapped.reads.sorted.bam

# Index the final bam file
samtools index mapped.reads.sorted.bam

# Call variance using BCFtools
bcftools mpileup -f contigs.bin3.fasta mapped.reads.sorted.bam | bcftools call --ploidy 1 -mv -o bcftools.variants.val.s18.vcf

# Call variance using HaplotypeCaller
gatk HaplotypeCaller -ploidy 1\
 -R contigs.bin3.fasta\
 -I mapped.reads.sorted.bam\
 -O HCaller.variants.val.vcf.gz\

# Call varince using Mutect2
gatk Mutect2 --af-of-alleles-not-in-resource 0.33 -R contigs.bin3.fasta -I mapped.reads.sorted.bam -O mutect2.variants.val.vcf.gz

# This is the end of script
