#!/bin/bash

# List of BAM files that failed indexing
BAM_FILES=(
    "/home/ubuntu/Work/HRA006113/HRR1458952/HRR1458952.bam"
    "/home/ubuntu/Work/HRA006113/HRR1458954/HRR1458954.bam"
    "/home/ubuntu/Work/HRA006113/HRR1458959/HRR1458959.bam"
    "/home/ubuntu/Work/HRA006113/HRR1459043/HRR1459043.bam"
    "/home/ubuntu/Work/HRA006113/HRR1459044/HRR1459044.bam"
    "/home/ubuntu/Work/HRA006113/HRR1459045/HRR1459045.bam"
    "/home/ubuntu/Work/HRA006113/HRR1459046/HRR1459046.bam"
    "/home/ubuntu/Work/HRA006113/HRR1459047/HRR1459047.bam"
    "/home/ubuntu/Work/HRA006113/HRR1459048/HRR1459048.bam"
    "/home/ubuntu/Work/HRA006113/HRR1459049/HRR1459049.bam"
    "/home/ubuntu/Work/HRA006113/HRR1459051/HRR1459051.bam"
    "/home/ubuntu/Work/HRA006113/HRR1459052/HRR1459052.bam"
    "/home/ubuntu/Work/HRA006113/HRR1459053/HRR1459053.bam"
    "/home/ubuntu/Work/HRA006113/HRR1459054/HRR1459054.bam"
    "/home/ubuntu/Work/HRA006113/HRR1459055/HRR1459055.bam"
    "/home/ubuntu/Work/HRA006113/HRR1459056/HRR1459056.bam"
)

# Process each BAM file
for bam_file in "${BAM_FILES[@]}"; do
    echo "Processing $bam_file..."
    
    # Get the directory and filename
    dir=$(dirname "$bam_file")
    filename=$(basename "$bam_file")
    name="${filename%.bam}"
    
    # Sort the BAM file
    echo "Sorting $bam_file..."
    samtools sort -@ 4 -o "${dir}/${name}.sorted.bam" "$bam_file"
    
    # Index the sorted BAM file
    echo "Indexing ${dir}/${name}.sorted.bam..."
    samtools index "${dir}/${name}.sorted.bam"
    
    echo "Completed processing $bam_file"
    echo "----------------------------------------"
done

echo "All files have been processed!" 