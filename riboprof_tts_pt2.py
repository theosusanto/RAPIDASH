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
log_dir = os.path.join(args.root_dir, "logs")

# if not already made, make directories for stages 3 and beyond
stage_dir = [os.path.join(args.root_dir, "stage"+str(i)) for i in range(0,10)]
for sd in stage_dir[3:10]:
    os.makedirs(sd, exist_ok=True)
#i.e. '/labs/mbarna/users/adelexu/s27l-flag-riboprof/stage0'

# rename files in stage2, such that the files with mismatched barcode is considered as "unmatched"
#for st2fq in glob.glob(f""):
#    fq_basename = os.path.basename(st2fq)
#    parts = fq_basename.split("_")
#    if "unmatched" in fq_basename:
#        continue
#    elif "barcode" in parts:
#        bc_i = parts.index(barcode)
#        common_phrase = "rep"
#        for string in parts:
#            if common_phrase in string:
#    else:
#        print(f"Unknown file {st2fq}")
        

# make sample directories in stage 2, and move stage 2 files there -- uncomment when running for first time
st2fq_suffix = "_stage2bcmatched.fastq"
for st2fq in glob.glob(f"{stage_dir[2]}/*{st2fq_suffix}"):
    sample_dir_short = os.path.basename(st2fq).replace(st2fq_suffix, "")
    sample_dir = os.path.join(stage_dir[2], sample_dir_short)
    os.makedirs(sample_dir, exist_ok=True)
    shutil.move(st2fq, sample_dir)

# find sample directories in stage 2 directory and make similar ones in stage 3-10 directories, excluding "unmatched"
sample_name = next(os.walk(stage_dir[2]))[1]
matched_name = [i for i in sample_name if "unmatched" not in i]
#i.e. '191217_lane3_2c4totrfp'

# uncomment if running for the first time
for m in matched_name:
    for stdr in stage_dir[3:10]:
        os.makedirs(os.path.join(stdr, m), exist_ok=True)

# names of input and output files at each step/stage
# note that step X takes in stage X and outputs stage X+1
stage2_fq = [os.path.join(stage_dir[2], m, m+"_stage2bcmatched.fastq") for m in matched_name]

stage3_log = [os.path.join(stage_dir[3], m, m+"_extract_log.txt") for m in matched_name]
stage3_fqgz = [os.path.join(stage_dir[3], m, m+"_stage3umibcextracted.fastq.gz") for m in matched_name]

stage4_fqgz = [os.path.join(stage_dir[4], m, m+"_stage4qualfiltered.fastq.gz") for m in matched_name]

rtsnRNA_index = "/labs/mbarna/index/rRNAtRNAsnRNAsnoRNA_human"
stage5_al_sam = [os.path.join(stage_dir[5], m, m+"_stage5alignedrtsnRNA.sam") for m in matched_name]
stage5_notal_fq = [os.path.join(stage_dir[5], m, m+"_stage5notalignedrtsnRNA.fastq") for m in matched_name]

canonical_index = "/labs/mbarna/index/gencode.v43.pc_transcripts.EnsemblCanonical.protein_coding"
stage6_al_sam = [os.path.join(stage_dir[6], m, m+"_stage6alignedcanonical.sam") for m in matched_name]
stage6_al_fq = [os.path.join(stage_dir[6], m, m+"_stage6alignedcanonical.fastq") for m in matched_name]
stage6_notal_fq = [os.path.join(stage_dir[6], m, m+"_stage6notalignedcanonical.fastq") for m in matched_name]

stage7_al_bam = [os.path.join(stage_dir[7], m, m+"_stage7alignedcanonical.bam") for m in matched_name]
stage7_sorted = [os.path.join(stage_dir[7], m, m+"_stage7sorted.bam") for m in matched_name]

stage8_stats = [os.path.join(stage_dir[8], m, m+"_stage8stats") for m in matched_name]
stage8_bam = [os.path.join(stage_dir[8], m, m+"_stage8dedup.bam") for m in matched_name]

stage9_bam = [os.path.join(stage_dir[9], m, m+"_stage9resort.bam") for m in matched_name]

sb_defaults = "-e %x-%j.e -o %x-%j.o -A mbarna -p batch"


# iterate over each stage2 fastq.gz input file
for i in range(len(stage2_fq)):

    #construct shell code
    sh_code = '\n'.join([f"#!/bin/sh -l",
                         f"date +'%D %T'",
                         f"echo",
                         f"echo '### Starting step 2: Extract UMIs and in-line barcodes'",
                         f"echo",
                         f"module load anaconda",
                         f"source activate riboprof_py310",
                         f"gzip {stage2_fq[i]}",
                         f"umi_tools extract --extract-method=string --bc-pattern=NNNNNCCCCC --log={stage3_log[i]} --3prime --stdin={stage2_fq[i]+'.gz'} --stdout={stage3_fqgz[i]}",
                         f"date +'%D %T'",
                         f"echo",
                         f"echo '### Starting step 3: Filter reads by quality'",
                         f"echo",
                         f"module load fastx_toolkit",
                         f"zcat {stage3_fqgz[i]} | fastq_quality_filter -Q33 -q 20 -p 70 -z -v -o {stage4_fqgz[i]}",
                         f"date +'%D %T'",
                         f"echo",
                         f"echo '### Starting step 4: align to rRNA, tRNA, snRNA'",
                         f"echo",
                         f"gunzip {stage4_fqgz[i]}",
                         f"module add bowtie2/2.5.1",
                         f"bowtie2 -L 18 -x {rtsnRNA_index} -q {stage4_fqgz[i].replace('.gz','')} -S {stage5_al_sam[i]} --un {stage5_notal_fq[i]}",
                         f"date +'%D %T'",
                         f"echo",
                         f"echo '### Starting step 5: align to canonical transcripts'",
                         f"echo",
                         f"bowtie2 --norc -L 18 -x {canonical_index} -q {stage5_notal_fq[i]} -S {stage6_al_sam[i]} --al {stage6_al_fq[i]} --un {stage6_notal_fq[i]}",
                         f"date +'%D %T'",
                         f"echo",
                         f"echo '### Starting step 6: convert to BAM, sort, and index'",
                         f"echo",
                         f"module load samtools",
                         f"samtools view -bS {stage6_al_sam[i]} > {stage7_al_bam[i]}",
                         f"samtools sort {stage7_al_bam[i]} -o {stage7_sorted[i]}",
                         f"samtools index {stage7_sorted[i]}",
                         f"date +'%D %T'",
                         f"echo",
                         f"echo '### Starting step 7: Dedup UMIs'",
                         f"echo",
                         f"umi_tools dedup -I {stage7_sorted[i]} --output-stats={stage8_stats[i]} -S {stage8_bam[i]}",
                         f"date +'%D %T'",
                         f"echo",
                         f"echo '### Starting step 8: re-sort and index BAM, convert back to SAM'",
                         f"echo",
                         f"samtools sort {stage8_bam[i]} -o {stage9_bam[i]}",
                         f"samtools index {stage9_bam[i]}",
                         f"samtools view -h {stage9_bam[i]} > {stage9_bam[i].replace('.bam', '.sam')}",
                         f"date +'%D %T'",
                         f"echo",
                         f"echo '### DONE part 2'"])

    # construct sbatch command line
    sb_cmd = f"sbatch {sb_defaults} -D {log_dir} -J {matched_name[i]}_{args.sb_name} -t {args.sb_time} --mem={args.sb_mem} -c {args.sb_cpus} <<EOF \n{sh_code}\nEOF"

    # call sbatch command line as subprocess, extract job ID
    sb_sub_msg = subprocess.check_output(sb_cmd, shell=True).decode('ascii')
    print(matched_name[i], sb_sub_msg)
    job_id = sb_sub_msg.strip().replace("Submitted batch job ", "")

    # dump config file
    with open(os.path.join(log_dir, f"{matched_name[i]}{args.sb_name}_riboprof-pt2_{job_id}.config"), 'w') as config_file:
        git_version = str(subprocess.check_output(['git', 'rev-parse', 'HEAD']).strip())
        configs = {"git version": git_version, "arguments": vars(args), "sbatch command line": sb_cmd}
        json.dump(configs, config_file, indent=4)
