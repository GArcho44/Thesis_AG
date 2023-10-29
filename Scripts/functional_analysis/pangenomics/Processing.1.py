from Bio import SeqIO
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'

# Provide the path to files
fasta_file = 'core_alignment.fasta'
gff_file = 'core_alignment.gff'

# Function to create the initial data frame #


def create_dataframe(f_file):
    sequences = list(SeqIO.parse(f_file, 'fasta'))
    records = SeqIO.parse(f_file, 'fasta')
    fasta_ids = [record.id for record in records]

    alignment_length = len(sequences[0].seq)

    data = []
    for i in range(alignment_length):
        column = [sequence.seq[i] for sequence in sequences]

        counts = {base: column.count(base) for base in 'ACGT'}
        most_frequent_base = max(counts, key=counts.get)    # Get the most frequent base

        # Create a site ID column containing pangenome, position and allele information

        site_id = "Pangenome|pos-" + str(i + 1) + "|ID1|" + str(most_frequent_base)

        # Create pangenome column

        pangenome = "pangenome"

        # Create a locus type column

        locus_type = "CDS"
        # Append the corresponding columns
        data.append((site_id, pangenome, i + 1, most_frequent_base, locus_type))


    df = pd.DataFrame(data, columns=['Site ID', 'Genome', 'Position', 'Ref_allele', 'Locus_type'])

    return df

# Function that introduces the base in alignment for each sample


def mg_allele(df, f_file):

    print("Step 1: Loading")

    sequences = list(SeqIO.parse(f_file, 'fasta'))
    records = SeqIO.parse(f_file, 'fasta')
    fasta_ids = [record.id for record in records]
    alignment_length = len(sequences[0].seq)

    for k in range(len(fasta_ids)):

        df["Ref_allele_" + str(fasta_ids[k])] = None

        for i in range(alignment_length):
            column = [sequence.seq[i] for sequence in sequences]
            df["Ref_allele_" + str(fasta_ids[k])][i] = column[k]

    print("Step 1: Completed")
    return df

# Apply the functions to the input data


initial_df = create_dataframe(fasta_file)
second_df = mg_allele(initial_df, fasta_file)


# Function to get the gene_ids from the annotated genome #


def gene_loci(gff_file):

    print("Step 2: Loading")
    with open(gff_file, 'r') as file:

        info_list = []
        for line in file:
            if not line.startswith("##"):
                start, end = line.strip().split()[3], line.strip().split()[4]   # extract start and end pos
                gene_id = line.strip().split()[8][3:9]  # extract gene identifiers
                gene_len = int(end) - int(start) + 1

                # Append the corresponding columns
                info_list.append((start, end, gene_id, gene_len))

        info_data = pd.DataFrame(info_list, columns=['start', 'end', 'gene_id', 'gene_len'])

        print("Step 2: Completed")
        return info_data


# Apply the functions to the input data
gene_id_df = gene_loci(gff_file)

# Function to connect information from the two dataframes


def connecting_info(df1, df2, f_file):

    print("Step 3: Loading")

    sequences = list(SeqIO.parse(f_file, 'fasta'))
    records = SeqIO.parse(f_file, 'fasta')
    fasta_ids = [record.id for record in records]
    alignment_length = len(sequences[0].seq)

    for q in range(alignment_length):
        column = [sequence.seq[q] for sequence in sequences]

    # Initiate columns
    df1['Gene_id'] = None
    df1['gene_len'] = None
    df1['pos_in_gene'] = None
    df1['pos_in_gene_r'] = None

    for l in range(len(df2)):
        n = 0
        for i in range(int(df2['start'][l])-1, int(df2['end'][l])):
                df1['Gene_id'][i] = df2['gene_id'][l]
                df1['gene_len'][i] = df2['gene_len'][l]
                n += 1
                df1['pos_in_gene'][i] = n   # get the position of nucleotide in gene

    for k in range(len(fasta_ids)):     # iterate in number of genomes
        df1["N_" + str(fasta_ids[k])] = None

        for l in range(len(df2)):   # iterate for number of genes
            num = 0

            for i in range(int(df2['start'][l]) - 1, int(df2['end'][l])):   # iterate in gene
                column = [sequence.seq[i] for sequence in sequences]
                if column[k] == "N":
                    num += 1

                df1["N_" + str(fasta_ids[k])][i] = num
        print("Sample " + str(k) + " completed")


    print("Step 3: Completed")
    return df1


# Apply the connecting info function
df = connecting_info(second_df, gene_id_df, fasta_file)

# Initialize an empty list to store the individual dataframes
dfs = []

# Write dataframes to CSV files
df.to_csv('data.csv', index=False)

print(df)


