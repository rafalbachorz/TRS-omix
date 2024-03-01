#!/bin/bash
#SBATCH --job-name=blastn_job
#SBATCH --output=blastn_job_%j.out
#SBATCH --error=blastn_job_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=12:00:00

# Load conda
source ~/miniconda3/etc/profile.d/conda.sh

# Activate your conda environment
conda activate blas
#Provide path to blast_db
blast_db="/home/hsalamaga/ncbi-blast-2.13.0+/bin/NCBI_nt_DB/nt.fa" 
# Check if query directory is provided as command-line argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <query_directory>"
    exit 1
fi

# Set query directory from command-line argument
query_dir="$1"

# Check if query directory exists
if [ ! -d "$query_dir" ]; then
    echo "Query directory does not exist: $query_dir"
    exit 1
fi
# Iterate over each query file in the directory
for query_file in "$query_dir"/*.fasta; do
    # Extract the filename without extension
    filename=$(basename -- "$query_file")
    filename_no_ext="${filename%.*}"

    # Define output file name
    output_file="$query_dir/${filename_no_ext}_blastn_out"

    # Run BLAST command for each query file
    blastn -query "$query_file" -db "$blast_db" -out "$output_file" -evalue 1e-5 -num_threads $SLURM_CPUS_PER_TASK -outfmt 6 -perc_identity 100
done
#Will output to filtered_sequences move the result files to blast_output, you don't have to rename them
