#!/usr/bin/env bash
#SBATCH --job-name=featurecounts_matrix               # Job name
#SBATCH --mail-type=BEGIN,END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=INSERT EMAIL HERE       # Where to send mail
#SBATCH --cpus-per-task=16                      # Number of CPU cores per task
#SBATCH --mem=64gb                             # Total memory
#SBATCH --time=12:00:00                         # Time limit hrs:min:sec
#SBATCH --output=/India_Nasal_Cells/Part_1/DATA/featureCounts_table/logs/featurecounts.out   # Standard output and error log
#SBATCH --array=1-79


CONDA_ENV="~/miniconda3/etc/profile.d/conda.sh"

source $CONDA_ENV
conda activate TB_pipeline #insert name of CONDA environment you installed the package on

INPUT_DIR="/India_Nasal_Cells/Part_1/DATA/alignments_output/alignments"
OUTPUT_DIR="/India_Nasal_Cells/Part_1/DATA/featureCounts_table"
GTF="/India_Nasal_Cells/Part_1/DATA/gencode.v49.primary_assembly.annotation.gtf"
LOG_DIR="/India_Nasal_Cells/Part_1/DATA/featureCounts_table/logs"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

cd $INPUT_DIR

echo "Starting FeatureCounts Matrix generation"

featureCounts \
	-T $SLURM_CPUS_PER_TASK \
	-s 2 \
	-p \
	-a "$GTF" \
       	-o "${OUTPUT_DIR}/gene_counts_matrix.txt" \
       	${INPUT_DIR}/*.sorted.bam

cd "/India_Nasal_Cells/Part_1"

Rscript cut_features.R $OUTPUT_DIR "gene_counts.tsv"

echo "Finished featureCounts successfully"