# This is a script for filtering the raw SNPs called from Mutect2

# Load appropriate libraries
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)

# Read Mutect2 raw calls
mutect2 <- fread(file='mutect2.variants.vcf', sep='\t', header = TRUE, skip = '#CHROM')
mutect2 <- as_tibble(mutect2)

# Read the whole vcf file
mutect2_full <- fread(file='mutect2.variants.vcf', sep='\t', header = FALSE)
# Read only the calls
mutect2_calls <- fread(file='mutect2.variants.vcf', sep='', header = FALSE, skip = '#CHROM')
# Store the rest part of the file (will be concatenated with the filtered SNPs calls)
mutect2_rest <- mutect2_full %>% filter(!V1 %in% mutect2_calls$V1)

## Filtering ##

# Filter for non-SNP variants
mutect2 <- mutect2 %>% filter(REF=='A'|REF=='T'|REF=='G'|REF=='C') %>% filter(ALT=='A'|ALT=='T'|ALT=='G'|ALT=='C')

# Filter based on alt allele frequency

# Get the correct information from the AS_SB_Table (INFO column)
mutect2$frequencies <- sub("DP.*", "", mutect2$INFO) 
mutect2$frequencies <- sub(".*AS_SB_TABLE=", "", mutect2$frequencies) 
mutect2$frequencies <- str_replace_all(mutect2$frequencies, "[[:punct:]]", " ")
mutect2$frequencies <- regmatches(mutect2$frequencies, gregexpr("[[:digit:]]+", mutect2$frequencies))

# Use this information to calculate the frequency of the alt allele per position
for (i in 1:nrow(mutect2)) {
  ref <- as.numeric(unlist(mutect2$frequencies[i])[1]) + as.numeric(unlist(mutect2$frequencies[i])[2])
  alt <- as.numeric(unlist(mutect2$frequencies[i])[3]) + as.numeric(unlist(mutect2$frequencies[i])[4])
  
  mutect2$perc[i] <- alt / (ref + alt)
  
}

# Here it is possible to set the maximum alternative allele frequency. Also remove the perc and frequencies columns
mutect2_filtered_tmp <- mutect2 %>% filter(perc < 0.9) %>% select(-perc, -frequencies)

# Write the filtered vcf file

# Write the filtered called snps in a temporary file (tmp)
fwrite(mutect2_filtered_tmp, file="mutect2.filtered.tmp.vcf", sep='\t')
# Read the temporary file 
mutect2_tmp <- fread(file='mutect2.filtered.tmp.vcf', sep='', header = FALSE)
# Join the rest with the temporary file to get the complete (final) vcf file
mutect2_filtered <- full_join(mutect2_rest, mutect2_tmp)
# Write the final file
fwrite(mutect2_filtered, file="mutect2.filtered.vcf", sep='\t')

# This is the end of the script