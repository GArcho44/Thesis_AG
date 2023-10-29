from Bio import SeqIO

def subset_gff3_with_sequence(gff3_file, fasta_file, output_file):
    # Read contig IDs from the FASTA file
    contig_ids = set()
    for record in SeqIO.parse(fasta_file, "fasta"):
        contig_ids.add(record.id)

    # Subset GFF3 annotations and sequence based on contig IDs
    gff3_annotations = []
    sequence_records = []
    is_sequence_section = False
    with open(gff3_file, "r") as gff3:
        for line in gff3:
            if line.startswith("#"):
                if line.startswith("##FASTA"):
                    is_sequence_section = True
                if not is_sequence_section:
                    parts_header = line.split()
                    contig_id_h = parts_header[1]
                    if contig_id_h in contig_ids:
                        gff3_annotations.append(line)  # Store header lines
            else:
                if not is_sequence_section:
                    parts = line.split("\t")
                    contig_id = parts[0]
                    if contig_id in contig_ids:
                        gff3_annotations.append(line)  # Store annotation lines

    # Print the number of contigs
    print("You have", len(contig_ids), "contigs")

    # Write the filtered annotations and sequence to a new GFF3 file
    with open(output_file, "w") as output:
        output.writelines(gff3_annotations)
        output.write("##FASTA\n")
        output.writelines(sequence_records)

# Usage example
gff3_file = "../../data/E.coli/gffs/ID2_d334_s88.gff"
fasta_file = "../../data/E.coli/Genomes/ID2_d334_s88.fasta"
output_file = "../../data/P.dorei/input_gffs/ID2_d334_s88.gff"
subset_gff3_with_sequence(gff3_file, fasta_file, output_file)

# Print message
print("Procedure completed!")

def append_fasta_to_file(fasta_file, target_file):
    with open(fasta_file, 'r') as fasta_handle:
        fasta_content = fasta_handle.read()

    with open(target_file, 'a') as target_handle:
        target_handle.write(fasta_content)

# Usage example
fasta_file = "../../data/E.coli/Genomes/ID2_d334_s88.fasta"
target_file = "../../data/P.dorei/input_gffs/ID2_d334_s88.gff"
append_fasta_to_file(fasta_file, target_file)