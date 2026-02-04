#!/usr/bin/env bash

# Run Trimmomatic with our parameters
srun run_trimmomatic.sh

# run the FastQC quality control script on the trimmed sequences
srun run_fastqc.sh

# This will aggregate all FastQC reports into a single MultiQC report
# Will will look at this report to determine if further cleaning is necessary
srun run_multiqc.sh