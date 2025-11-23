# MethylExtractor

A high-performance tool for extracting DNA methylation data from bisulfite sequencing BAM files, featuring advanced overlapping read pair handling and comprehensive quality control.

## ✨ Key Features

- **🔬 Advanced Read Processing**: Proper handling of overlapping paired-end reads to prevent double-counting
- **⚡ High Performance**: Multi-threaded processing with deterministic results
- **📊 Comprehensive Analysis**: CpG, CHG, and CHH methylation contexts with detailed statistics
- **💾 Flexible Output**: HDF5, text, or combined formats with compression
- **🚀 Efficient Compression**: Zstd compression with gzip fallback for HDF5 output
- **🎯 Quality Control**: Rigorous filtering with MAPQ ≥ 30 and Phred ≥ 20 defaults

## Installation

### Prerequisites
- GCC compiler
- HTSlib (automatically installed by Makefile)
- HDF5 library (automatically installed by Makefile)
- **Zstd** (required for optimized compression)

### Build from Source
```bash
git clone <repository-url>
cd MethylExtractor
make
```

The Makefile automatically installs required dependencies (on Debian/Ubuntu systems) and compiles the optimized binary.

## Usage

### Basic Syntax
```bash
bin/MethylExtractor [options] <input.bam> <output_directory> [reference.fa]
```
*(Note: Binary is located in `build/dynamic/x64/MethylExtractor` or similar, depending on architecture)*

### Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `-h, --help` | Show help message and exit | N/A |
| `-t, --threads INT` | Number of processing threads | 16 |
| `-q, --min-mapq INT` | Minimum mapping quality (MAPQ) | 30 |
| `-p, --min-phred INT` | Minimum base quality (Phred score) | 20 |
| `-c, --min-cov INT` | Minimum coverage threshold | 4 |
| `-C, --cap-cov INT` | Cap coverage to this value (0 = disabled) | 0 |
| `-G, --CHG` | Include CHG methylation contexts | Disabled |
| `-H, --CHH` | Include CHH methylation contexts | Disabled |
| `-m, --chrom-mapping FILE` | Chromosome mapping configuration file | `chrom_mapping.json` |
| `-z, --compression INT` | Compression level (0-8). Zstd with gzip fallback | 6 |
| `-k, --chunk-size INT` | HDF5 chunk size for I/O optimization | 1,000,000 |
| `-f, --output-format STR` | Output format: `hdf5`, `txt`, or `both` | `hdf5` |
| `-s, --split` | Split output by methylation context | Disabled |
| `-o, --output-dir DIR` | Output directory | N/A |

### 🚀 Performance Optimization (New!)

For the best balance of speed and compression ratio, use **compression level 9**:

```bash
./MethylExtractor input.bam output_dir -z 9 -t 32
```

This triggers a special optimization:
1. Writes uncompressed HDF5 data extremely fast.
2. Automatically launches a multi-threaded `zstd` process to compress the file in the background.
3. Result: **~5x faster processing** with maximum compression.

### Chromosome Mapping File

When processing specific chromosomes, create a `chrom_mapping.json` file:

```json
{
  "reference": "/path/to/genome.fa",
  "chromosomes": [
    {
      "name": "chr1",
      "fasta": "chr1",
      "bam": "chr1",
      "extract": true
    }
  ]
}
```

## Output Formats

### HDF5 Format (Default)
Efficient binary format with compression:
```
Dataset: methylation_data
├── pos (uint32) - Genomic position (1-based)
├── mC (uint16) - Methylated cytosine count
├── uC (uint16) - Unmethylated cytosine count
└── tnc (uint8) - Trinucleotide context and strand
```

### Text Format
Human-readable tab-separated values:
```bash
# position strand methylated_count unmethylated_count context trinucleotide
161	+	11	9	CG	CGT
310	-	8	12	CG	CGA
```

### Output Files

| Format Option | Files Generated | Description |
|---------------|-----------------|-------------|
| `-f hdf5` | `chr1.h5`, `chr2.h5`, ... | Per-chromosome HDF5 files |
| `-f txt` | `chr1.txt`, `chr2.txt`, ... | Per-chromosome text files |
| `-f both` | Both HDF5 and text files | Combined output |
| `-s --split` | `chr1.CG.h5`, `chr1.CHG.h5`, ... | Separate files per context |

### Statistics Output

Each run generates a JSON statistics file with genome-wide metrics:

```json
{
  "num_positions": 1345382,
  "total_methylated": 5084976,
  "total_unmethylated": 18665653,
  "avg_methylation_level": 0.214,
  "avg_coverage": 17.65
}
```

## Examples

### Basic Usage
```bash
# Process entire BAM file with defaults
./MethylExtractor input.bam output_dir reference.fa
```

### High-Quality Analysis
```bash
# Strict quality control, all contexts, both formats
./MethylExtractor \
    -t 16 \
    -q 30 \
    -p 20 \
    -c 4 \
    -G \
    -H \
    -f both \
    -s \
    input.bam \
    output_dir \
    reference.fa
```

### Production Pipeline
```bash
# Optimized for large-scale analysis
./MethylExtractor \
    --threads 32 \
    --min-mapq 30 \
    --min-phred 20 \
    --min-cov 5 \
    --cap-cov 100 \
    --CHG \
    --CHH \
    --output-format hdf5 \
    --compression 9 \
    sample.bam \
    results/
```

### Custom Chromosome Selection
```bash
# Process specific chromosomes only
./MethylExtractor \
    -m custom_chromosomes.json \
    -t 8 \
    input.bam \
    output_dir
```

## Quality Control Features

MethylExtractor implements rigorous quality control:

- **Read-level filtering**: MAPQ ≥ 30, removes duplicates, secondary alignments, and multimappers
- **Base-level filtering**: Phred ≥ 20 quality scores
- **Bisulfite validation**: Ensures reads match expected conversion patterns
- **Overlapping read handling**: Prevents double-counting in paired-end data
- **Coverage normalization**: Optional capping for PCR bias correction

## Performance Guidelines

### Threading
- **Small datasets (< 10GB)**: 4-8 threads
- **Medium datasets (10-50GB)**: 8-16 threads
- **Large datasets (> 50GB)**: 16-32 threads

### Memory Usage
- Streaming I/O prevents loading entire BAM files
- Memory scales with chromosome size, not genome size
- HDF5 chunked storage optimizes I/O performance

## Dependencies

- **HTSlib**: BAM/SAM file processing
- **HDF5**: Binary data storage and compression
- **Zlib/Bzip2/LZMA**: Compression libraries
- **GCC**: Compiler with OpenMP support
- **Zstd**: High-performance compression with fallback to gzip

## Troubleshooting

### Common Issues
- **Empty output**: Check chromosome mapping and BAM headers
- **High memory usage**: Reduce thread count or process fewer chromosomes
- **Slow processing**: Use SSD storage and adjust chunk sizes. **Try `-z 9` for faster compression.**
- **Low methylation**: Verify bisulfite conversion and quality filters

### Performance Optimization
- Use SSD storage for input/output
- Adjust `--chunk-size` based on chromosome sizes
- Balance `--compression` level with I/O performance needs

## Documentation

For complete documentation, see `MethylExtractor_Documentation.html`.

## License

This project is licensed under the MIT License.