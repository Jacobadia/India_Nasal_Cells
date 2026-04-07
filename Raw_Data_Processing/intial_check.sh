#!/usr/bin/env bash
#SBATCH --job-name=quality_control
#SBATCH --output=/India_Nasal_Cells/Part_1/DOCUMENTS/FASTQC1/ERROR/polyG_%A_%a.out
#SBATCH --error=/India_Nasal_Cells/Part_1/DOCUMENTS/FASTQC1/ERROR/polyG_%A_%a.out
#SBATCH --mail-user=INSERT EMAIL HERE
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --array=1-79

#----------
CONDA_ENV="~/miniconda3/etc/profile.d/conda.sh"

source $CONDA_ENV

HOME_DIR="/India_Nasal_Cells/Part_1"
INPUT_DIR="/India_Nasal_Cells/Part_1/DATA/fastq" #input fastq file directory
OUTPUT_DIR="/India_Nasal_Cells/Part_1/DATA/fastqc_pretrim" #fastqc folder
SAMPLE_LIST="/India_Nasal_Cells/Part_1/DATA/sample_idsPOLYG.txt"

mkdir -p "$OUTPUT_DIR"
mkdir -p "/India_Nasal_Cells/Part_1/DOCUMENTS/FASTQC1/ERROR"
mkdir -p "/India_Nasal_Cells/Part_1/DOCUMENTS/FASTQC1/OUTPUT_LOGS"

SAMPLE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
R1="${SAMPLE_ID}.R1.fastq.gz"
R2="${SAMPLE_ID}.R2.fastq.gz"

echo "Processing Task $SLURM_ARRAY_TASK_ID: $SAMPLE_ID"

cd $HOME_DIR

conda activate TB_pipeline

fastqc "${INPUT_DIR}/$R1" "${INPUT_DIR}/$R2" -t $SLURM_CPUS_PER_TASK -o "$OUTPUT_DIR"

cd $OUTPUT_DIR

unzip -o "${SAMPLE_ID}.R1_fastqc.zip"
unzip -o "${SAMPLE_ID}.R2_fastqc.zip"

echo "Finished Sample $SAMPLE_ID"

conda deactivate TB_pipeline