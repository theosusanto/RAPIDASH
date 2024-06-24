import argparse
import gzip
import csv
from Bio import SeqIO
import re

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Extract selected records from a FASTA file based on TSV file.')
parser.add_argument('tsv_file', help='Path to the TSV file')
parser.add_argument('search_column', help='Column header to search in the TSV file')
parser.add_argument('fasta_gz_file', help='Path to the compressed FASTA file')
parser.add_argument('output_file', help='Path to the output file')
args = parser.parse_args()

# Extract command-line arguments
tsv_file = args.tsv_file
search_column = args.search_column
fasta_gz_file = args.fasta_gz_file
output_file = args.output_file

# Find the index for search_column
with open(tsv_file, 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    header_row = next(reader)
    try:
        column_index = header_row.index(search_column)
    except ValueError:
        print(f"Column '{search_column}' not found in the TSV file.")
        exit(1)

    # Extract the desired genes from the TSV file
    desired_genes = [row[column_index] for row in reader]

# Open the compressed file using gzip and parse it with Biopython
with gzip.open(fasta_gz_file, 'rt') as file:
    # Parse the FASTA file using Biopython SeqIO
    records = SeqIO.parse(file, 'fasta')

    # Create a list comprehension to filter and select the desired records
    selected_records = [record for record in records if re.search(r'(\w+)\.\d+', record.description).group(1) in desired_genes]

# Write the selected records to the output file using SeqIO.write()
SeqIO.write(selected_records, output_file, 'fasta')
