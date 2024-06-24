# modified 1/19/20

import glob, os, subprocess, argparse

# Example: python /labs/mbarna/users/gtiu/scripts/parse_start_distance_wrapper.py /labs/mbarna/users/gtiu/scripts/mm10_knownCanonical_UTR_CDS_length.txt

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sam_dir')
parser.add_argument('--out_dir')
parser.add_argument('--log_dir')
parser.add_argument('--script_path')
parser.add_argument("--utrcdslengthpath")
args = parser.parse_args()

stage9_suffix = '_RFP_stage9resort.sam'

for file in glob.glob(os.path.join(args.sam_dir, '*', f'*{stage9_suffix}')):    
    print(file)
    samp_name_short = os.path.basename(file).replace(stage9_suffix, '')
    os.makedirs(os.path.join(args.out_dir, samp_name_short))
    jobname = f'{samp_name_short}_startdist'
    sb_cmd = f"sbatch -J {jobname} -t 6:00:00 --mem=8000 -D {args.log_dir} -e %x-%j.e -o %x-%j.o -A mbarna -p batch {args.script_path} {file} {args.utrcdslengthpath}"
    subprocess.call(sb_cmd, shell = True)
