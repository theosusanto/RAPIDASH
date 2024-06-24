import argparse
import os
import subprocess
import json
import glob

parser = argparse.ArgumentParser(description = "Run part 1 of ribosome profiling pipeline (through barcode splitting step)")
parser.add_argument("--root_dir")
parser.add_argument("--sb_name")
parser.add_argument("--sb_time", default = "12:0:0")
parser.add_argument("--sb_mem", default = "32G")
parser.add_argument("--sb_cpus", default = "4")
args = parser.parse_args()

# make log directory inside root directory
log_dir = os.path.join(args.root_dir, "logsrna")
os.makedirs(log_dir, exist_ok=True)

## make directories for each stage from 1 to 10 (stage0 has already been created so skip making it)
stage_dir = [os.path.join(args.root_dir, "stage"+str(i)+"rna") for i in range(0,10)]
for sd in stage_dir[1:10]:
    os.makedirs(sd, exist_ok=True)

## get lane names and prefixes from stage0 directory.
stage0_fqgz = glob.glob(os.path.join(stage_dir[0], "*_stage0raw.fastq.gz"))
stage0_read1_fqgz = [path for path in stage0_fqgz if path.endswith("1_stage0raw.fastq.gz")]
dir_fqgz = [os.path.basename(os.path.dirname(i)) for i in stage0_fqgz]
basename_prefix=[os.path.basename(i).replace("stage0raw.fastq.gz","") for i in stage0_fqgz]

## make lane directories in stage 1 directory
#for ln in lane_dir_short:
#    os.makedirs(os.path.join(stage_dir[1], ln), exist_ok=True)

stage1_fqgz = [os.path.join(stage_dir[1], basename_prefix[i]+"stage1trimmed.fastq.gz") for i in range(len(basename_prefix))]

# This is i7 adapter sequence
adapter3p = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
adapter5p = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

## default sbatch parameters
sb_defaults = "-e %x-%j.e -o %x-%j.o -A mbarna -p batch"
# %x: Job name, %j: jobid of the running job
# -e|--error <filename_pattern>. Connect batch script's standard error directly to the file name specified.
# -o|--output <filename_pattern>. Connect batch script's standard output directly to the file name specified. By default, both standard output and standard error are directed to the same file.
# -A|--acount <account>. mbarna
# -p|--partition <partition_names>. Request a specific partition for resource allocation. For scg4, there are 3: batch, nih_s10, and interactive. Refer to: https://login.scg.stanford.edu/tutorials/job_scripts/

## iterate over each stage0 fastq.gz input file

# cutadapt:
#   -j|--cores <cores>. Number of cores. 0 refers to autodetection.
#   -u|--cut <number>. Remove a fixed number of bases. Unconditionally remove bases from the beginning (if the given length is positive) or the end (if the given length is negative) of each read. It is applied before any adapter trimming. I removed the NN upstream of RPF.
#   --discard-untrimmed. Discard reads in which no adapter was found.
#   -m|--minimum-length <length>. Discard processed reads that are shorter than <length>.
#   -o|--output <path>. Specify output path.

for i in range(len(stage0_read1_fqgz)):
    stage0_read1_path = stage0_read1_fqgz[i]
    stage0_read2_path = stage0_read1_path.replace("_1_", "_2_")
    stage1_basename = os.path.basename(stage0_read1_path).replace("_stage0raw.fastq.gz","_stage1trimmed.fastq.gz")
    stage1_read1_path = os.path.join(stage_dir[1], stage1_basename)
    stage1_read2_path = stage1_read1_path.replace("_1_", "_2_")
    #construct shell code
    sh_code = '\n'.join([f"#!/bin/sh -l",
                         f"date +'%D %T'",
                         f"echo",
                         f"echo '### Starting step 0: trim reads'",
                         f"echo",
                         f"module load cutadapt/3.4",
                         f"Use Read 1 adapter sequence",
                         f"cutadapt -a {adapter3p} -A {adapter5p} -o {stage1_read1_path} -p {stage1_read2_path} {stage0_read1_path} {stage0_read2_path}",
                         f"echo '### DONE part 1'"])
    # construct sbatch command line
    sb_cmd = f"sbatch {sb_defaults} -D {log_dir} -J {basename_prefix[i]}{args.sb_name} -t {args.sb_time} --mem={args.sb_mem} -c {args.sb_cpus} <<EOF \n{sh_code}\nEOF"
    #print(sh_code)
    # call sbatch command line as subprocess, extract job ID
    sb_sub_msg = subprocess.check_output(sb_cmd, shell=True).decode('ascii')
    print(sb_sub_msg)
    job_id = sb_sub_msg.strip().replace("Submitted batch job ", "")

    # dump config file

    with open(os.path.join(log_dir, f"{basename_prefix[i]}{args.sb_name}_rnaseq-pt1-tts_{job_id}.config"), 'w') as config_file:
        git_version = str(subprocess.check_output(['git', 'rev-parse', 'HEAD']).strip())
        configs = {"git version": git_version, "arguments": vars(args), "sbatch command line": sb_cmd}

        json.dump(configs, config_file, indent=4)
