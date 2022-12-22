# Random mutations script #
# This is a python script that makes random single nucleotide mutations in the selected genomes

# Import packages
import csv
import pandas as pd
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Create empty lists where the data will be stored
sequences_list = []
contigs_list = []
mutated = []
record_list = []
mutated_pos = []

# Read the fasta file and store the corresponding info (seq, id)
for record in SeqIO.parse("contigs.bin3.fasta", "fasta"):
    sequences_list.append(record.seq)
    contigs_list.append(record.id)


# Create a function that introduced mutations based on a pre-selected cut-off
def mutate(sequence, mutation_rate):
    # Convert sequence to list
    dna_list = list(sequence)

    # Create a list to store the mutation positions
    mutation_list = []
    for n in range(len(sequence)):
        r = random.random()
        if r < mutation_rate:
            mutation_site = random.randint(0, len(dna_list) - 1)
            mutation_list.append(mutation_site)

            # Make the mutation depending on the nucleotide
            if dna_list[mutation_site] == 'A':
                dna_list[mutation_site] = random.choice(list('TCG'))
            elif dna_list[mutation_site] == 'T':
                dna_list[mutation_site] = random.choice(list('ACG'))
            elif dna_list[mutation_site] == 'C':
                dna_list[mutation_site] = random.choice(list('ATG'))
            else:
                dna_list[mutation_site] = random.choice(list('ATC'))

    return ''.join(dna_list), mutation_list


# Call the function for the list of contigs (fasta sequences) for 1% mutation rate
mutated_c = [mutate(s, 0.002) for s in sequences_list]
# The result contains a string of the new fasta sequence with the mutations, and a list of the mutation positions

# Create a list of the mutated contigs
for i in range(len(mutated_c)):
    mutated.append(mutated_c[i][0])

# Create a list of the mutated positions
for i in range(len(mutated_c)):
    mutated_pos.append(mutated_c[i][1])
# Create a SeqRecord using the id and mutated sequences
for i in range(len(mutated)):
    record_list.append((SeqRecord(Seq(mutated[i]), id=contigs_list[i], description='')))

# Output
# Write the fasta file (fasta)
SeqIO.write(record_list, "mutatedg_j.fasta", "fasta")

# Write the position file (fasta)
# Dictionary of lists
dict_list = {'CONTIG': contigs_list, 'Mutated.sites': mutated_pos}
df = pd.DataFrame(dict_list)
df.to_csv('positions_j.csv')
# <-->
