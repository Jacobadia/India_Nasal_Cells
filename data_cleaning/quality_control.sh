#!/usr/bin/env bash

#SBATCH --job-name=quality_control
#SBATCH --output=/grphome/grp_tb/processing_scripts/quality_control.out
#SBATCH --error=/grphome/grp_tb/processing_scripts/quality_control.err
#SBATCH --mail-user=bryant.steed@gmail.com
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --mem=150G

# Run Trimmomatic with our parameters
trimmomatic.sh

# run the FastQC quality control script on the trimmed sequences
# this should aggregate the summary.txt files too
fastqc_trimmed.sh

# run multiQC. If we decide to add this