#!/bin/bash

# Reference genome path
REF_GENOME="/home/ubuntu/Work/genomes/human_genome/release-113/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# Chromosome mapping file
CHROM_MAPPING="/home/ubuntu/Work/HRA006113/chrom_mapping.json"

# Base directory containing BAM files
BASE_DIR="/home/ubuntu/Work/HRA006113"

# Find all BAM files
find "$BASE_DIR" -name "*.bam" | while read -r bam_file; do
    # Get the directory name (sample ID)
    sample_dir=$(dirname "$bam_file")
    sample_id=$(basename "$sample_dir")
    
    echo "Processing $sample_id..."
    
    # Create output directory
    output_dir="$sample_dir"
    
    # Run MethylExtractor with parameters from launch.json
    /home/ubuntu/MethylExtractor/build/dynamic/MethylExtractor \
        --chrom-mapping "$CHROM_MAPPING" \
        --@ 16 \
        --o "$output_dir" \
        --hdf5-compression 9 \
        --q 30 \
        --p 20 \
        --CHG \
        --CHH \
        --c 0 \
        --chunk-size 1000000 \
        --no-cap-coverage \
        --split-context-files \
        "$REF_GENOME" \
        "$bam_file"
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "Successfully processed $sample_id"
    else
        echo "Error processing $sample_id"
    fi
    
    echo "----------------------------------------"
done 