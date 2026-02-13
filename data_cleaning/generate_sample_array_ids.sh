#!/usr/bin/env bash

# This creates a tab-delimited file with sample array IDs and 
# their corresponding R1 and R2 fastq file paths, Similar to kahmbati's scripts
# This will let us reference samples by array ID in downstream scripts to do things
# in Parallel.
# We'll need to change the paths here to wherever our fastq files are located and the
# output file we want to create.
# We also only need to run this once

fastq_dir="/grphome/grp_tb/fastqB"
output_file="sample_array_ids.txt"

i=1
for r1filename in "${fastq_dir}"/*.trim.paired_r1.fastq.gz; do
    r2filename="${r1filename/.trim.paired_r1.fastq.gz/.trim.paired_r2.fastq.gz}"
    sample_name=$(basename "${r1filename}" .trim.paired_r1.fastq.gz)
    echo "${sample_name}"$'\t'"${sample_name}.trim.paired_r1.fastq.gz"$'\t'"${sample_name}.trim.paired_r2.fastq.gz" >> $output_file
    ((i++))
done