from Bio import AlignIO
import csv

# Load the aligned sequences from a file (replace 'alignment.fasta' with your file)
alignment = AlignIO.read("core_alignment.fasta", "fasta")

# Initialize an empty identity matrix
identity_matrix = {}

# Iterate through pairs of sequences
for i, seq1 in enumerate(alignment):
    for j, seq2 in enumerate(alignment):
        # Calculate the identity as the percentage of matching positions
        identity = sum(1 for a, b in zip(seq1.seq, seq2.seq) if a == b) / len(seq1.seq)

        # Store the identity in the matrix
        identity_matrix[(i, j)] = identity

# Define the output CSV file name
output_file = "IM_ID2_E.coli.csv"

# Open the CSV file for writing
with open(output_file, mode='w', newline='') as csv_file:
    writer = csv.writer(csv_file)

    # Write the header row with sequence IDs
    sequence_ids = [seq.id for seq in alignment]
    writer.writerow([''] + sequence_ids)

    # Write the identity values row by row
    for i, seq1 in enumerate(alignment):
        row_data = [sequence_ids[i]]
        for j in range(len(alignment)):
            row_data.append(identity_matrix[(i, j)])
        writer.writerow(row_data)

print(f"Identity matrix has been saved to {output_file}")






