#!/bin/bash

module load anaconda
source activate riboprof_py310

input_file=$1
output_file=$2

python /home/tsusanto/git/ribosome_profiling/UTR_CDS_reads_count.py $input_file $output_file

