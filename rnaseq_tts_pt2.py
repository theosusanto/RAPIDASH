### pipeline part 2

import argparse
import glob
import os
import shutil
import subprocess
import json

parser = argparse.ArgumentParser(description = "Run part 2 of ribosome profiling pipeline (all steps performed on individual samples)")
parser.add_argument("--root_dir")
parser.add_argument("--sb_name")
parser.add_argument("--sb_time", default = "12:0:0")
parser.add_argument("--sb_mem", default = "32G")
parser.add_argument("--sb_cpus", default = "4")
args = parser.parse_args()

# expect logs directory to have been made in pipeline part 1
log_dir = os.path.join(args.root_dir, "logsrna")

# if not already made, make directories for stages 3 and beyond
stage_dir = [os.path.join(args.root_dir, "stage"+str(i))+"rna" for i in range(0,10)]
for sd in stage_dir[2:10]:
    os.makedirs(sd, exist_ok=True)
#i.e. '/labs/mbarna/users/adelexu/s27l-flag-riboprof/stage0'

# make sample directories in stage 1, and move stage 1 files there -- uncomment when running for first time
st1fq_suffix = "_stage1trimmed.fastq.gz"
for st1fq in glob.glob(f"{stage_dir[1]}/*_1{st1fq_suffix}"):
    sample_dir_short = os.path.basename(st1fq).replace("_1"+st1fq_suffix, "")
    sample_dir = os.path.join(stage_dir[1], sample_dir_short)
    os.makedirs(sample_dir, exist_ok=True)
    st1fq2 = st1fq.replace("_1_stage1trimmed", "_2_stage1trimmed")
    shutil.move(st1fq, sample_dir)
    shutil.move(st1fq2, sample_dir)

# find sample directories in stage 1 directory and make similar ones in stage 2-10 directories
sample_name = next(os.walk(stage_dir[1]))[1]
matched_name = [i for i in sample_name if "unmatched" not in i]
#i.e. '191217_lane3_2c4totrfp'
# uncomment if running for the first time
for m in matched_name:
    for stdr in stage_dir[2:10]:
        os.makedirs(os.path.join(stdr, m), exist_ok=True)

# names of input and output files at each step/stage
# note that step X takes in stage X and outputs stage X+1
stage1_fq_1 = [os.path.join(stage_dir[1], m, m+"_1_stage1trimmed.fastq.gz") for m in matched_name]
stage1_fq_2 = [os.path.join(stage_dir[1], m, m+"_2_stage1trimmed.fastq.gz") for m in matched_name]

rtsnRNA_index = "/labs/mbarna/index/rRNAtRNAsnRNAsnoRNA_human"
stage2_al_sam = [os.path.join(stage_dir[2], m, m+"_stage2alignedrtsnRNA.sam") for m in matched_name]
stage2_notal_fq = [os.path.join(stage_dir[2], m, m+"_stage2notalignedrtsnRNA.fastq.gz") for m in matched_name]

canonical_index = "/labs/mbarna/index/gencode.v43.pc_transcripts.EnsemblCanonical.protein_coding"
stage3_al_sam = [os.path.join(stage_dir[3], m, m+"_stage3alignedcanonical.sam") for m in matched_name]
stage3_al_fq = [os.path.join(stage_dir[3], m, m+"_stage3alignedcanonical.fastq.gz") for m in matched_name]
stage3_notal_fq = [os.path.join(stage_dir[3], m, m+"_stage3notalignedcanonical.fastq.gz") for m in matched_name]

stage4_al_bam = [os.path.join(stage_dir[4], m, m+"_stage4alignedcanonical.bam") for m in matched_name]
stage4_sorted = [os.path.join(stage_dir[4], m, m+"_stage4sorted.bam") for m in matched_name]
stage4_fixmate = [os.path.join(stage_dir[4], m, m+"_stage4fixmate.bam") for m in matched_name]
stage4_dedup = [os.path.join(stage_dir[4], m, m+"_stage4dedup.bam") for m in matched_name]
stage4_dedupsorted = [os.path.join(stage_dir[4], m, m+"_stage4dedupsorted.bam") for m in matched_name]

stage5_output = [os.path.join(stage_dir[5], m) for m in matched_name]

sb_defaults = "-e %x-%j.e -o %x-%j.o -A mbarna -p batch"

# iterate over each stage2 fastq.gz input file
for i in range(len(stage1_fq_1)):
    #construct shell code
    sh_code0 = '\n'.join([f"#!/bin/sh -l",
                          f"date +'%D %T'",
                          f"echo '### Starting step 1: align to rRNA, tRNA, snRNA'",
                          f"echo",
                          f"module add bowtie2/2.5.1",
                          f"bowtie2 -L 18 -x {rtsnRNA_index} -1 {stage1_fq_1[i]} -2 {stage1_fq_2[i]} -S {stage2_al_sam[i]} --un-conc-gz {stage2_notal_fq[i]}",
                          f"echo",
                          f"date +'%D %T'",
                          f"echo",
                          f"echo '### Starting step 2: align to canonical transcripts'",
                          f"echo",
                          f"bowtie2 -L 18 -x {canonical_index} -1 {stage2_notal_fq[i].replace('.fastq.gz', '.fastq.1.gz')} -2 {stage2_notal_fq[i].replace('.fastq.gz', '.fastq.2.gz')} -S {stage3_al_sam[i]} --al-conc-gz {stage3_al_fq[i]} --un-conc-gz {stage3_notal_fq[i]}",
                          f"date +'%D %T'",
                          f"echo",
                          f"echo '### Starting step 3: convert to BAM, markdup, sort, and index'",
                          f"echo",
                          f"module load samtools/1.17",
                          f"samtools view -bS {stage3_al_sam[i]} > {stage4_al_bam[i]}",
                          f"samtools sort {stage4_al_bam[i]} -o {stage4_sorted[i]}",
                          f"samtools index {stage4_sorted[i]}",
                          f"samtools fixmate -m {stage4_sorted[i]} {stage4_fixmate[i]}"
                          f"samtools markdup -r -s {stage4_fixmate[i]} {stage4_dedup[i]}",
                          f"samtools sort {stage4_dedup[i]} -o {stage4_dedupsorted[i]}",
                          f"samtools index {stage4_dedupsorted[i]}",
                          f"date +'%D %T'",
                          f"echo",
                          f"echo '### DONE'"])
    sh_code1 = '\n'.join([f"#!/bin/sh -l",
                         f"module load samtools/1.17",
                         f"samtools collate {stage4_sorted[i]} -O | samtools fixmate -m - - | samtools sort - -o {stage4_fixmate[i]}",
                         f"samtools markdup -r -s {stage4_fixmate[i]} {stage4_dedup[i]}",
                         f"samtools sort {stage4_dedup[i]} -o {stage4_dedupsorted[i]}",
                         f"samtools index {stage4_dedupsorted[i]}",
                         f"echo '### DONE'"])
    sh_code2 = '\n'.join([f"#!/bin/sh -l",
                        f"module load samtools/1.17",
                        f"samtools view -h {stage4_dedupsorted[i]} > {stage4_dedupsorted[i].replace('dedupsorted.bam', 'resort.sam')}",
                        f"echo '### DONE'"])
    sh_code = '\n'.join([f"#!/bin/sh -l",
                         f"module load salmon/1.4.0",
                         f"salmon quant -t /labs/mbarna/index/gencode.v43.pc_transcripts.EnsemblCanonical.protein_coding.fa.gz -l A -a {stage4_dedupsorted[i]} -o {stage5_output[i]}"])

    # construct sbatch command line
    sb_cmd = f"sbatch {sb_defaults} -D {log_dir} -J {matched_name[i]}_{args.sb_name} -t {args.sb_time} --mem={args.sb_mem} -c {args.sb_cpus} <<EOF \n{sh_code}\nEOF"

    # call sbatch command line as subprocess, extract job ID
    sb_sub_msg = subprocess.check_output(sb_cmd, shell=True).decode('ascii')
    print(matched_name[i], sb_sub_msg)
    job_id = sb_sub_msg.strip().replace("Submitted batch job ", "")

    # dump config file
    with open(os.path.join(log_dir, f"{matched_name[i]}{args.sb_name}_rnaseq-tts-pt2_{job_id}.config"), 'w') as config_file:
        git_version = str(subprocess.check_output(['git', 'rev-parse', 'HEAD']).strip())
        configs = {"git version": git_version, "arguments": vars(args), "sbatch command line": sb_cmd}
        json.dump(configs, config_file, indent=4)
