import gzip
import argparse
from Bio import SeqIO

# Function to calculate the length of a field identified in the description line
def get_field_length(field, field_prefix):
    if field.startswith(field_prefix):
        position_range = field.split(":")[1]
        start, end = map(int, position_range.split("-"))
        return end - start + 1
    else:
        return 0

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('fasta_gz_file', help='Path to the compressed FASTA file')
parser.add_argument('output_tsv_file', help='Path to the output TSV file')
args = parser.parse_args()

# Open the compressed file using gzip and parse it with Biopython
with gzip.open(args.fasta_gz_file, 'rt') as file, open(args.output_tsv_file, 'w') as output_file:
    # Parse the FASTA file using Biopython SeqIO
    records = SeqIO.parse(file, 'fasta')
    
    # Write the header to the output TSV file
    output_file.write("transcript\tl_tr\tl_utr5\tl_cds\tl_utr3\n")
    
    # Iterate over the records and extract information from the description
    for record in records:
        description = record.description
        fields = description.split("|")
        utr5_length = 0
        cds_length = 0 
        utr3_length = 0
        transcript_length = 0
        # Extract the length of the fields
        for field in fields:
            utr5_length += get_field_length(field, "UTR5:")
            cds_length += get_field_length(field, "CDS:")
            utr3_length += get_field_length(field, "UTR3:")
            transcript_length = utr5_length + cds_length + utr3_length
        transcriptID = description

        # Write the row to the output TSV file
        output_file.write(f"{transcriptID}\t{transcript_length}\t{utr5_length}\t{cds_length}\t{utr3_length}\n")
