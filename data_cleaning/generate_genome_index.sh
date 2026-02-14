#!/usr/bin/env bash

#SBATCH --job-name=genome_index
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/grphome/grp_tb/processing_scripts/genome_index.out
#SBATCH --mail-user=bryant.steed@gmail.com
#SBATCH --mail-user=18briankim@gmail.com
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --mem=160G

genome_index_dir="/grphome/grp_tb/star_genome/index"
genome_fasta="/grphome/grp_tb/star_genome/GRCh38.primary_assembly.genome.fa"
annotations_gtf="/grphome/grp_tb/star_genome/gencode.v49.primary_assembly.annotation.gtf"
# The documentation recommends setting this to max(read_length) 
# they say that 99 is a good default for illumina reads that are 100 bp. 
# However, I ran the max read length 
# script on our samples and got 151 every time
overhang=150

conda_environment="tb_processing"

if [ ! -r "$genome_fasta" ]; then
    echo "Error: Genome FASTA file not found or not readable: $genome_fasta"
    exit 1
fi

if [ ! -r "$annotations_gtf" ]; then
    echo "Error: Annotations GTF file not found or not readable: $annotations_gtf"
    exit 1
fi

if [ ! -w "$genome_index_dir" ]; then
    echo "Error: Genome index directory not writable: $genome_index_dir"
    exit 1
fi

conda activate "${conda_environment}"

STAR \
--runThreadN "${SLURM_CPUS_PER_TASK}" \
--runMode genomeGenerate \
--genomeDir "${genome_index_dir}" \
--genomeFastaFiles "${genome_fasta}" \
--sjdbGTFfile "${annotations_gtf}" \
--sjdbOverhang "${overhang}"