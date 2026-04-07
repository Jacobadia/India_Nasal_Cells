#!/usr/bin/env bash

#SBATCH --job-name=quality_control2
#SBATCH --output=/India_Nasal_Cells/Part_1/DOCUMENTS/FASTQ2/ERROR/quality_control_%A_%a.out
#SBATCH --error=/India_Nasal_Cells/Part_1/DOCUMENTS/FASTQ2/OUTPUT_LOGS/quality_control_%A_%a.err
#SBATCH --mail-user=INSERT EMAIL HERE
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --array=1-79

CONDA_ENV="~/miniconda3/etc/profile.d/conda.sh"
source $CONDA_ENV

HOME_DIR="/India_Nasal_cells/Part_1"

INPUT_DIR="/India_Nasal_cells/Part_1/DATA/fastq" #fastq files
OUTPUT_DIR="/India_Nasal_cells/Part_1/DATA/alignment_input" #fastq files trimmed
OUTPUT_DIR2="/India_Nasal_cells/Part_1/DATA/fastqc_posttrim" #trimmed fastqc directory
SAMPLE_LIST="/India_Nasal_cells/Part_1/DATA/sample_array_ids.txt"

ERROR_DIR="/India_Nasal_Cells/Part_1/DOCUMENTS/FASTQ2/ERROR"
OUT_DIR="/India_Nasal_Cells/Part_1/DOCUMENTS/FASTQ2/OUTPUT_LOGS"

mkdir -p "$ERROR_DIR" "$OUT_DIR"

 #the file with the array of sample IDs in the data folder

# Trimming
CURRENT_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST" | tr -d '\r')
SAMPLE_ID=$(echo "$CURRENT_LINE" | awk '{print $1}')
R1_FILE=$(echo "$CURRENT_LINE" | awk '{print $2}')
R2_FILE=$(echo "$CURRENT_LINE" | awk '{print $3}')

R1_TRIM="${OUTPUT_DIR}/${SAMPLE_ID}.clean.R1.fastq.gz"
R2_TRIM="${OUTPUT_DIR}/${SAMPLE_ID}.clean.R2.fastq.gz"

echo "### Starting Pipeline for Sample: $SAMPLE_ID (Task $SLURM_ARRAY_TASK_ID) ###"

conda activate fastp

fastp \
  -i "${INPUT_DIR}/${R1_FILE}" \
  -I "${INPUT_DIR}/${R2_FILE}" \
  -o "${OUTPUT_DIR}/${SAMPLE_ID}.clean.R1.fastq.gz" \
  -O "${OUTPUT_DIR}/${SAMPLE_ID}.clean.R2.fastq.gz" \
  --detect_adapter_for_pe \
  --trim_poly_g \
  --trim_poly_x \
  --cut_front \
  --cut_tail \
  --cut_window_size 4 \
  --cut_mean_quality 20 \
  --length_required 36 \
  --thread $SLURM_CPUS_PER_TASK \
  --html "${OUTPUT_DIR}/${SAMPLE_ID}_fastp.html" \
  --json "${OUTPUT_DIR}/${SAMPLE_ID}_fastp.json"

conda deactivate

conda activate TB_pipeline

cd $OUTPUT_DIR #fastQ files but trimmed

fastqc "$R1_TRIM" "$R2_TRIM" \
    -t "$SLURM_CPUS_PER_TASK" \
    -o "$OUTPUT_DIR2"

conda deactivate

echo "Finished Sample $SAMPLE_ID"