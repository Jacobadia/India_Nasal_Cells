#!/usr/bin/env bash

# This creates a tab-delimited file with sample array IDs and 
# their corresponding R1 and R2 fastq file paths, Similar to kahmbati's scripts
# This will let us reference samples by array ID in downstream scripts to do things
# in Parallel.
# We'll need to change the paths here to wherever our fastq files are located and the
# output file we want to create.
# We also only need to run this once

$fastq_dir="$(pwd)/fastq"
$output_file="$(pwd)/sample_array_ids.txt"

i=1
for r1filename in $(ls $fastq_dir/*_R1.fastq.gz); do
    r2filename="${r1filename/_R1.fastq.gz/_R2.fastq.gz}"
    echo "${i}"$'\t'"${r1filename}"$'\t'"${r2filename}" >> $output_file
    ((i++))
done