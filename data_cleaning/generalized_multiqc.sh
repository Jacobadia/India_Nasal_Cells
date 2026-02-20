#!/usr/bin/env bash

# This is a generalized version of our run_multiqc.sh script
# To show how we used multiqc to generate Quality Control Reports
# The trimmed and untrimmed directories are passed as arguments to the script.
# They initially contain Raw FastQC out. This script is just for generating
# The figure, which appears in the output HTML file.

trimmed_dir="${1}"
trimmed_output_dir="${2}"
untrimmed_dir="${3}"
untrimmed_output_dir="${4}"

cd ${trimmed_dir}
multiqc . -o ${trimmed_output_dir}

cd ${untrimmed_dir}
multiqc . -o ${untrimmed_output_dir}