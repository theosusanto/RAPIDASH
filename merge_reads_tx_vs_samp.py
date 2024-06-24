import csv

tsv_file1 = 'reads_tx_vs_samp.txt'  # Replace with the path to the first TSV file
tsv_file2 = 'reads_tx_vs_samp_rna.txt'  # Replace with the path to the second TSV file
output_file = 'combined_file.txt'  # Replace with the desired output file path

data1 = {}
data2 = {}

header1 = []
header2 = []

# Read data from the first TSV file and store it in a dictionary
with open(tsv_file1, 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    header1 = next(reader)
    for row in reader:
        key = row[0]
        values = row[1:]
        data1[key] = values

# Read data from the second TSV file and store it in a dictionary
with open(tsv_file2, 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    header2 = next(reader)
    for row in reader:
        key = row[0]
        values = row[1:]
        data2[key] = values

merged_header = header1 + header2[1:]

# Combine the keys from both dictionaries, removing duplicates
combined_keys = set(data1.keys()).union(data2.keys())

# Combine the data from both files and write to the output file
with open(output_file, 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerow(merged_header)

    # Iterate over the combined keys
    for key in combined_keys:
        values1 = data1.get(key, ['0.0'] * (len(header1) - 1))
        values2 = data2.get(key, ['0.0'] * (len(header2) - 1))
        combined_values = values1 + values2
        writer.writerow([key] + combined_values)
