# This is a script to preprocess the pangenomics data (part 2)

# Load appropriate libraries
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)

## Read the data ##

# Read the curated data set
c.dataset <- read.csv("data.csv")

# Read the tab loci list
loci.list <- read_delim("loci_list.tab",col_names=FALSE) %>% select(X1,X2,X5) 
colnames(loci.list) <- c("gene_region","Gene_id","sample")

# Function to get the complement DNA strand 
get_complement <- function(nucleotide) {
  complement <- switch(nucleotide,
                       "A" = "T",
                       "T" = "A",
                       "C" = "G",
                       "G" = "C",
                       nucleotide)
  return(complement)
}

# Processing 

# Bring the concatenated file in the appropriate format
# Make a column indicating the sample

# Join the two data frames with Gene_id and filter out paralogs 
com.df <- full_join(loci.list,c.dataset,by="Gene_id") %>% 
  filter(!is.na(Site.ID)) %>% distinct() %>% 
  pivot_wider(names_from=sample,values_from = gene_region, values_fn=list) %>%
  filter(rowSums(sapply(.[20:29], function(x) lengths(x) > 1)) == 0) %>%
  select(-Gene_id)

# Un-list the elements of the data frame
df <- data.frame(lapply(com.df, function(x) unlist(x, recursive = FALSE)))

process_samples <- function(sample_id,df) {
  
  
  # Read the coordinate file
  coord_file <- sample_id
  ori.coord <- read_delim(coord_file, col_names = TRUE) %>% 
    select(-Gene, -Type, -Product)
  sample_id <- gsub(".co.*","", sample_id)
  
  # Select data for the sample
  x <- df %>% select(sample_id, Site.ID, pos_in_gene,
                     paste0("Ref_allele_", sample_id, sep=""), paste0("N_", sample_id, sep="") , gene_len)
  ori.coord <- ori.coord %>% filter(Name %in% x[,1])
  colnames(ori.coord)[1] <- sample_id
  
  # Join data frames and calculate new coordinates
  coord.file <- full_join(x, ori.coord)
  colnames(coord.file)[5] <- "N"
  coord.filen <- coord.file %>% mutate(
    POS = ifelse(Strand == "Reverse",
                 End - pos_in_gene + N + 1,
                 Start + pos_in_gene - N)
  )
  
  # Read the SNPs file
  snps <- read.csv(paste0("S_",sample_id,".csv", sep=""), row.names = 1)
  colnames(snps) <- c("contig", "POS", "REF", "ALT", "freq")
  colnames(coord.filen)[11] <- "contig"
  
  # Perform inner join
  joint <- inner_join(coord.filen, snps, by = c("contig", "POS"))
  
  # Apply complement if strand is reverse
  joint$REF <- ifelse(joint$Strand == "Reverse", sapply(joint$REF, get_complement), joint$REF)
  joint$ALT <- ifelse(joint$Strand == "Reverse", sapply(joint$ALT, get_complement), joint$ALT)
  
  
  final.df <- full_join(coord.filen, joint)
  
  return(final.df)
}
  
## Call the function and get the final output ##

# Call the function for each sample

file_names <- list.files("./", pattern = "co-ords.tab")

for (sample_id in file_names) {
  
  data <- process_samples(sample_id, df) %>% select(Site.ID, ALT, freq)
  sample_id <- gsub(".co.*","", sample_id)
  var_name <- paste0("f_", sample_id, sep="")
  colnames(data) <- c("Site.ID", paste0("Minor_allele_", sample_id, sep=""), paste0("Minor_allele_freq_", sample_id, sep = ""))
  assign(var_name, data)
  
}

# Merge the data frames 

# first put all the data into a list
df_list <- list(f_ID1_d016_s109, f_ID1_d062_s75, f_ID1_d096_s52, f_ID1_d152_s7, f_ID1_d212_s69, f_ID1_d243_s46, f_ID1_d243_s46, f_ID1_d302_s1, f_ID1_d332_s12, f_ID1_d355_s36) # (+)

# then merge them into a new dataframe 
merged_data <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list) 

# Create the final data frame by combining the initial information with the above data frame
final.df <- Reduce(function(x, y) merge(com.df[, 1:14], merged_data, all=TRUE), df_list) 

# Sort by position
final.df <- final.df[order(final.df$Position), ]

# Write into a csv file
write.csv(final.df, file="./final.data.csv")
