#!/bin/bash

# Modules
module load legacy/.base
module load python/2.7

filename=$1
utrcdslengthpath=$2

# Parse start of read relative to UTR, CDS
python /home/tsusanto/git/ribosome_profiling/bowtieSAM_parse_start_position.py ${filename} ${utrcdslengthpath} ${filename/.sam/_start.txt} >> ${filename/.sam/_start_distance.log}

# Calculate distance of read start to start and stop codons
python /home/tsusanto/git/ribosome_profiling/distanceFromStartStop_size.py ${filename/.sam/_start.txt} ${filename/.sam/_distance.txt}
