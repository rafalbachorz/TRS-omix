#!/bin/bash

conda init bash
# Activate your conda environment
conda activate blast
#blast_db="/media/data/blast_db/nt"

# Set query directory from command-line argument
query_dir="$QUERY_DIR"
blast_db="$BLAST_DB"
output_dir="$OUTPUT_DIR"

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
    output_file="$output_dir/${filename_no_ext}_blastn_out"

    # Run BLAST command for each query file
    blastn -query "$query_file" -db "$blast_db" -out "$output_file" -evalue 1e-5 -num_threads 16 -outfmt 6 -perc_identity 100 -max_target_seqs 50000
done
