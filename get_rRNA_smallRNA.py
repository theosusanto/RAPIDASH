import argparse
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def extract_gene_coordinates(annotation_file, gene_types):
    gene_coordinates = {}

    with gzip.open(annotation_file, 'rt') as file:
        for line in file:
            if line.startswith('#'):
                continue
            columns = line.split('\t')
            feature_type = columns[2]
            if feature_type in gene_types:
                chromosome = columns[0]
                feature = columns[2]
                start = int(columns[3])
                end = int(columns[4])
                strand = columns[6]
                attributes = columns[8]

                gene = {
                    'chromosome': chromosome,
                    'feature': feature,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'attributes': attributes
                }

                gene_name = extract_gene_name(attributes)
                gene_coordinates.setdefault(gene_name, []).append(gene)

    return gene_coordinates

def extract_gene_name(attributes):
    # Extract gene name from attributes string based on GFF format
    # You may need to adjust this function based on your GFF file format
    # Example implementation:
    name_start = attributes.find('gene=') + 5
    name_end = attributes.find(';', name_start)
    gene_name = attributes[name_start:name_end]
    return gene_name

def write_sequences_to_fasta(output_file, gene_coordinates, fna_file):
    with gzip.open(fna_file, 'rt') as fna, open(output_file, 'w') as fasta:
        sequences = SeqIO.to_dict(SeqIO.parse(fna, "fasta"))

        records = []
        for gene_name, gene_list in gene_coordinates.items():
            for gene in gene_list:
                chromosome = gene['chromosome']
                start = gene['start']
                end = gene['end']
                strand = gene['strand']
                attributes = gene['attributes']

                if chromosome in sequences:
                    sequence = sequences[chromosome].seq
                    gene_sequence = sequence[start - 1:end]

                    if strand == '-':
                        gene_sequence = gene_sequence.reverse_complement()

                    record_id = gene_name + '_' + gene['feature']  # Include feature type in the record ID
                    record = SeqRecord(gene_sequence, id=record_id, description='')
                    records.append(record)

        SeqIO.write(records, fasta, "fasta")

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Extract gene sequences from GFF annotation file and FNA file.')
parser.add_argument('gff_file', help='Path to the GFF annotation file (.gz)')
parser.add_argument('fna_file', help='Path to the FNA file (.gz)')
parser.add_argument('output_fasta', help='Path to the output FASTA file')
args = parser.parse_args()

gene_types = ['tRNA', 'rRNA', 'snoRNA', 'snRNA']

gene_coordinates = extract_gene_coordinates(args.gff_file, gene_types)
write_sequences_to_fasta(args.output_fasta, gene_coordinates, args.fna_file)
