import glob, os, subprocess, argparse, json, re

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--parsedA_dir')
parser.add_argument('--log_dir')
parser.add_argument('--script_path', default='/home/tsusanto/git/ribosome_profiling/UTR_CDS_reads_count_submit.sh')
args = parser.parse_args()

parsedA_suffix_re = '_parseAposition_.+\.txt'

for file in glob.glob(os.path.join(args.parsedA_dir, '*', '*_parseAposition_*.txt')):

    print(file)
    
    samp_name_short = re.sub(parsedA_suffix_re, '', os.path.basename(file))
    
    print(samp_name_short)
    
    out_file = re.sub(parsedA_suffix_re, '_readcounts.txt', file)
    
    print(out_file)
    
    jobname = f'{samp_name_short}_count_reads'
    
    sb_cmd = f"sbatch -J {jobname} -t 1:00:00 --mem=8000 -D {args.log_dir} -e %x-%j.e -o %x-%j.o -A mbarna -p batch {args.script_path} {file} {out_file}"       
    
    subprocess.call(sb_cmd, shell = True)

#with open(os.path.join(args.log_dir, f"UTR_CDS_reads_count.config"), 'w') as config_file:
#        configs = {"arguments": vars(args)}
#        json.dump(configs, config_file, indent=4)
