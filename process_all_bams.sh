#!/bin/bash
#
# Batch-process every BAM under BASE_DIR with MethylExtractor.
# Output for each sample is written next to its BAM file.

set -euo pipefail

# Reference genome path
REF_GENOME="/work/genomes/human_genome/release-114/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# Chromosome mapping file
CHROM_MAPPING="/work/chrom_mapping.json"

# Base directory containing BAM files
BASE_DIR="/work/samples"

# MethylExtractor binary (built with `make`)
METHYL_EXTRACTOR="/home/ubuntu/MethylExtractor/build/dynamic/arm64/MethylExtractor"

# Find all BAM files
find "$BASE_DIR" -name "*.bam" | while read -r bam_file; do
    sample_dir=$(dirname "$bam_file")
    sample_id=$(basename "$sample_dir")

    echo "Processing $sample_id..."

    # Output directory (per sample)
    output_dir="$sample_dir"

    # CLI: MethylExtractor [options] <input.bam> [output_dir] [ref.fa]
    "$METHYL_EXTRACTOR" \
        --chrom-mapping "$CHROM_MAPPING" \
        --threads 16 \
        --min-mapq 30 \
        --min-phred 20 \
        --min-cov 0 \
        --compression 9 \
        --chunk-size 1000000 \
        --CHG \
        --CHH \
        --split \
        "$bam_file" \
        "$output_dir" \
        "$REF_GENOME"

    if [ $? -eq 0 ]; then
        echo "Successfully processed $sample_id"
        echo "Extraction manifest: ${output_dir}/${sample_id}.extraction_manifest.json"
        echo "Post-extraction Pass/Fail QC is run by MethylPipeline using the manifest and JSON sidecars."
    else
        echo "Error processing $sample_id"
    fi

    echo "----------------------------------------"
done
