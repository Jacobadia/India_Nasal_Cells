#/usr/bin/env bash

input_dir="/grphome/grp_tb/fastqT"

for file in "$input_dir"/*.fastq.gz; do
  if [[ -f "$file" ]]; then
    echo "Processing $file"
    maxlen=$(zcat "$file" | awk 'NR%4==2 {print length($0)}' | sort -nr | head -1)
    echo "Max read length in $file: $maxlen"
  else
    echo "No .fastq.gz files found in $input_dir"
  fi
done