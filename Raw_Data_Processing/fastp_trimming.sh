#!/usr/bin/env bash
#SBATCH --job-name=Clean_up              # Job name
#SBATCH --mail-type=BEGIN,END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=INSERT EMAIL HERE        # Where to send mail
#SBATCH --cpus-per-task=16                      # Number of CPU cores per task
#SBATCH --mem=64gb                             # Total memory
#SBATCH --time=6:00:00                         # Time limit hrs:min:sec
#SBATCH --output=/India_Nasal_Cells/Part1/DOCUMENTS/trimming_errors/trim_%A_%a.out   # Standard output and error log
#SBATCH --array=1-79

CONDA_ENV="~/miniconda3/etc/profile.d/conda.sh"

source $CONDA_ENV
conda activate fastp

INPUT_DIR="$1"
OUTPUT_DIR="$2"
SAMPLE_LIST="$3"
ERROR_LOG_DIR="/India_Nasal_Cells/Part1/DOCUMENTS/trimming_errors"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$ERROR_LOG_DIR"

CURRENT_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

SAMPLE_ID=$(echo "$CURRENT_LINE" | awk '{print $1}')
R1_FILE=$(echo "$CURRENT_LINE" | awk '{print $2}')
R2_FILE=$(echo "$CURRENT_LINE" | awk '{print $3}')

if [ ! -f "${INPUT_DIR}/${R1_FILE}" ]; then
    echo "Error: Input file ${INPUT_DIR}/${R1_FILE} not found."
    exit 1
fi

echo "Processing Array Job $SLURM_ARRAY_TASK_ID: Sample $SAMPLE_ID"

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

echo "Finished Sample $SAMPLE_ID"