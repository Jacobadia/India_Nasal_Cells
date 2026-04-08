#!/usr/bin/env bash
#SBATCH --job-name=alignment_tb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=INSERT EMAIL HERE
#SBATCH --cpus-per-task=28
#SBATCH --mem=160gb
#SBATCH --time=72:00:00
#SBATCH --output=/India_Nasal_Cells/Part_1/DATA/alignments_output/logs/alignment_%A_%a.out
#SBATCH --array=1-79

# adjust all of the SBATCH parameters above to your specifications
# there are 79 samples in the DATA

echo "### Setting paths ###"

# make sure to have sample_ids, grch38 genome, and annotation files in from DATA folder

DATA_DIR="/India_Nasal_Cells/Part_1/DATA/alignment_input" #trimmed fastqs
SAMPLE_LIST="/India_Nasal_Cells/Part_1/DATA/sample_idsPOLYG.txt" #found in the DATA folder
OUTPUT_DIR="/India_Nasal_Cells/Part_1/DATA/alignments_output" #finalized bam files

CONDA_ENV="~/miniconda3/etc/profile.d/conda.sh"

genome_index_dir="/India_Nasal_Cells/Part_1/DATA/star_genome/index" #in DATA folder
annotations_gtf="/India_Nasal_Cells/Part_1/DATA/star_genome/gencode.v49.primary_assembly.annotation.gtf" #in DATA folder

echo "Data Path: $DATA_DIR"
echo "Sample Sheet: $SAMPLE_LIST"
echo "Results Dir: $OUTPUT_DIR"

mkdir -p "$OUTPUT_DIR/alignments"
mkdir -p "$OUTPUT_DIR/logs"
mkdir -p "$OUTPUT_DIR/counts"

# Get sample ID from array
SAMPLE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

SAMPLE_ID=$(echo "$SAMPLE_ID" | tr -d '\r')

read_1="${SAMPLE_ID}.clean.R1.fastq.gz"
read_2="${SAMPLE_ID}.clean.R2.fastq.gz"

echo "Processing Task $SLURM_ARRAY_TASK_ID: $SAMPLE_ID"
echo "Read 1: $read_1"
echo "Read 2: $read_2"

# Verify input files exist
if [[ ! -f "$DATA_DIR/$read_1" ]]; then
    echo "ERROR: Missing $read_1"
    exit 1
fi

if [[ ! -f "$DATA_DIR/$read_2" ]]; then
    echo "ERROR: Missing $read_2"
    exit 1
fi

# Temporary working directory
temp_dir="/India_Nasal_Cells/Part_1/aligner_temp/${SAMPLE_ID}"
mkdir -p "$temp_dir"

echo "Copying FASTQ files to temp directory"
cp "$DATA_DIR/$read_1" "$temp_dir/"
cp "$DATA_DIR/$read_2" "$temp_dir/"

cd "$temp_dir"

echo "### Loading STAR environment ###"
source $CONDA_ENV
conda activate star_aligner

echo "### Running STAR alignment ###"

# The index dir is the directory containing the STAR genome index files, it
# has things necessary for the alignment. like the suffix array etc ...
# runThreadN is the number of threads to use, we set this to be the number
# of cores requested in the SBATCH parameters.
# quantMode GeneCounts tells STAR to count the number of reads mapping to each gene
# outFilterMultimapNmax 1 tells STAR to only keep reads that map to a
# single location in the genome (i.e. unique mappers)
# readFilesIn are the two paired read files for this sample.
# Remember that SLURM is executing this asynchronously across all our
# samples.
# We use zcat to read the gzipped fastq files directly without needing to unzip them first.
# we output the alignment in a BAM format, and we use
# sam tools to convert it to the correct file format this is critical for
# getting our featureCounts table later on.

STAR \
  --genomeDir "$genome_index_dir" \
  --sjdbGTFfile "$annotations_gtf" \
  --runThreadN "${SLURM_CPUS_PER_TASK}" \
  --quantMode GeneCounts \
  --outFilterMultimapNmax 1 \
  --readFilesIn "$read_1" "$read_2" \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix "./"

echo "### Indexing BAM ###"
samtools index "$temp_dir/Aligned.sortedByCoord.out.bam"

echo "### Moving results ###"

mv "$temp_dir/Aligned.sortedByCoord.out.bam" \
   "$OUTPUT_DIR/alignments/${SAMPLE_ID}.sorted.bam"

mv "$temp_dir/Aligned.sortedByCoord.out.bam.bai" \
   "$OUTPUT_DIR/alignments/${SAMPLE_ID}.sorted.bam.bai"

mv "$temp_dir/Log.final.out" \
   "$OUTPUT_DIR/logs/${SAMPLE_ID}.txt"

mv "$temp_dir/ReadsPerGene.out.tab" \
   "$OUTPUT_DIR/counts/${SAMPLE_ID}.txt"

echo "### Cleaning temporary directory ###"
rm -rf "$temp_dir"

echo "### $SAMPLE_ID completed successfully ###"