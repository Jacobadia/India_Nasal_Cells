#!/usr/bin/env bash
#SBATCH --job-name=alignment_tb               # Job name
#SBATCH --mail-type=BEGIN,END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=bryant.steed@gmail.com
#SBATCH --mail-user=18briankim@gmail.com        # Where to send mail
#SBATCH --cpus-per-task=28                      # Number of CPU cores per task
#SBATCH --mem=160gb                             # Total memory
#SBATCH --time=48:00:00                         # Time limit hrs:min:sec
#SBATCH --output=/grphome/grp_tb/processing_scripts/alignment_tb.out   # Standard output and error log
#SBATCH --array=1-79                            # Create an array of jobs (one for each sample)



# uses mnt/scratch as tempdir, with projects as back up tempdir
# does all 79 samples
# no unmapped output
# need requeue- some samples may be preempted and need to be run again.

echo "### Setting paths ###"
# Set paths
data_path="/grphome/grp_tb/processing_scripts/results/trimmed"
sample_sheet="/grphome/grp_tb/processing_scripts/results/documents/sample_array_ids.txt"

echo "Data Path: $data_path"
echo "Sample Sheet: $sample_sheet"

# Read the specific line of the sample sheet corresponding to the task ID
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $sample_sheet)

# Parse the line
IFS=$'\t' read -r -a sample_array <<< "$sample"
sample_name="${sample_array[0]}"
read_1="${sample_array[1]}"
read_2="${sample_array[2]}"

echo "### Processing sample: $sample_name ###"
echo "Read 1: $read_1"
echo "Read 2: $read_2"

echo "### Checking available scratch space ###"

echo "$HOSTNAME"

temp_dir="/grphome/grp_tb/trimmed_temp/${sample_name}"

mkdir -p "$temp_dir"

cp /grphome/grp_tb/processing_scripts/results/trimmed/${read_1} "$temp_dir"
cp /grphome/grp_tb/processing_scripts/results/trimmed/${read_2} "$temp_dir"
echo "Temporary Directory: $temp_dir"

cd $temp_dir

echo "### unzipping fastq ###"
# Unzip fastq files
gunzip -f "$temp_dir/${read_1}"
gunzip -f "$temp_dir/${read_2}"

unzipped_read_1="${read_1%.gz}"
unzipped_read_2="${read_2%.gz}"

echo "### load STAR and related packages ###"
# Load STAR and samtools
conda activate star_aligner

echo "# running STAR #"
# Run STAR alignment

annotations_gtf="/grphome/grp_tb/star_genome/gencode.v49.primary_assembly.annotation.gtf"
genome_index_dir="/grphome/grp_tb/star_genome/index"
max_locations_mapped_to=1
max_mismatches_per_read=10

STAR \
--genomeDir "$genome_index_dir" \
--sjdbGTFfile "$annotations_gtf" \
--runThreadN "${SLURM_CPUS_PER_TASK}" \
--quantMode GeneCounts \
--outFilterMultimapNmax "$max_locations_mapped_to" \
--outFilterMismatchNmax "$max_mismatches_per_read" \
--readFilesIn "$temp_dir/$unzipped_read_1" "$temp_dir/$unzipped_read_2" \
--outFileNamePrefix "$temp_dir/"

echo "### cleaning up: deleting unzipped FASTQ files ###"
# Delete the unzipped FASTQ files to free up space
rm "$temp_dir/$unzipped_read_1"
rm "$temp_dir/$unzipped_read_2"
echo "### cleaned out fastq successfully ###"

echo "# Convert SAM to BAM #"
samtools view -bS "$temp_dir/Aligned.out.sam" > "$temp_dir/Aligned.out.bam"
echo "### SAM converted successfully ###"

echo "### cleaning up: deleting SAM file ###"
# rm "$temp_dir/Aligned.out.sam"
echo "### deleted SAM successfully ###"

echo "# Sort BAM and index"
samtools sort "$temp_dir/Aligned.out.bam" -o "$temp_dir/Aligned.out.sorted.bam"
samtools index "$temp_dir/Aligned.out.sorted.bam"
echo "### BAM sorted and indexed successfully ###"

echo "# tidying up #"
# Move outputs to the results directories
results_dir="/grphome/grp_tb/processing_scripts/results/star_output"
mv "$temp_dir/Aligned.out.sorted.bam" "$results_dir/alignments/${sample_name}.sorted.bam"
mv "$temp_dir/Aligned.out.sorted.bam.bai" "$results_dir/alignments/${sample_name}.sorted.bam.bai"
mv "$temp_dir/Log.final.out" "$results_dir/logs/${sample_name}.txt"
mv "$temp_dir/ReadsPerGene.out.tab" "$results_dir/counts/${sample_name}.txt"


echo "# $sample_name completed alignment #"

echo "## clean temporary directory ##"
rm -rf $temp_dir

echo "# $HOSTNAME or projects temp directory deleted successfully #"
