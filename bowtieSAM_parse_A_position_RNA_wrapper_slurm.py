import glob, os, subprocess, argparse, json

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sam_dir')
parser.add_argument('--sample_type_tag', default='rna')
parser.add_argument('--out_dir')
parser.add_argument('--log_dir')
parser.add_argument('--script_path')
parser.add_argument("--utrcdslengthpath", default='/labs/mbarna/users/gtiu/scripts/mm10_knownCanonical_GENCODE_VM20_highconf_UTR_CDS_length.txt')
args = parser.parse_args()

stage9_suffix = '_stage9resort.sam'

for file in glob.glob(os.path.join(args.sam_dir, '*', f'*{args.sample_type_tag}*{stage9_suffix}')):

    print(file)
    
    samp_name_short = os.path.basename(file).replace(stage9_suffix, '')
    
    out_file = file.replace(stage9_suffix, '_parseAposition_RNA.txt')
    
    jobname = f'{samp_name_short}_parseAposition_RNA'
    sb_cmd = f"sbatch -J {jobname} -t 1:00:00 --mem=8000 -D {args.log_dir} -e %x-%j.e -o %x-%j.o -A mbarna -p batch {args.script_path} {file} {args.utrcdslengthpath} {out_file}"       
    
    subprocess.call(sb_cmd, shell = True)

with open(os.path.join(args.log_dir, f"bowtieSAM_parse_A_position_RNA.config"), 'w') as config_file:
        configs = {"arguments": vars(args)}
        json.dump(configs, config_file, indent=4)