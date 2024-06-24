import argparse
import glob
import os
import subprocess
import csv

parser = argparse.ArgumentParser(description = "Parse and collate counts from salmon quant.sf files")
parser.add_argument("--root_dir")
parser.add_argument("--tsv_file")
parser.add_argument("--out_file")
args = parser.parse_args()

log_dir = os.path.join(args.root_dir, "logsrna")
os.makedirs(log_dir, exist_ok=True)

tsv_file = args.tsv_file
dict_counts = {}

# check where salmon output is
stage_dirname = "stage3rna"
stage_dir = os.path.join(args.root_dir, stage_dirname)

sample_names = next(os.walk(stage_dir))[1]
quant_samples = [os.path.join(stage_dir, m, "quant.sf") for m in sample_names]


# create a dictionary of transcript names
with open(tsv_file, 'r') as file:
    reader = csv.reader(file, delimiter = '\t')
    header_row = next(reader)
    for row in reader:
        if len(row) > 0:
            trid = row[0].split(".")[0]
            dict_counts[trid] = [0]*len(quant_samples)

sample_colnames = []
for i, qt in enumerate(quant_samples):
    sample_colnames.append(qt.split("/")[-2])
    with open(qt, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        header_row = next(reader)
        for row in reader:
            if len(row) > 0:
                trid = row[0].split("|")[0].split(".")[0]
                counts = row[4].split(".")[0]+".0"
                if trid in dict_counts:
                    dict_counts[trid][i] = counts

# write an output with the collated counts
output_file_path = args.out_file
with open(output_file_path, 'w', newline = '') as file:
    writer = csv.writer(file, delimiter = '\t')
    writer.writerow(['tx_id'] + sample_colnames)
    for key, values in dict_counts.items():
        writer.writerow([key] + values)

