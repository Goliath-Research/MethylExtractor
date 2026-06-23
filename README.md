# MethylExtractor

A high-performance tool for extracting DNA methylation data from bisulfite sequencing BAM files, featuring advanced overlapping read pair handling and comprehensive quality control.

## ✨ Key Features

- **🔬 Advanced Read Processing**: Coordinate-based clipping of overlapping paired-end mates to prevent double-counting
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
make deps   # one-time: install system dependencies (Debian/Ubuntu, uses sudo)
make        # compile the optimized binary
```

`make deps` installs the required system packages (on Debian/Ubuntu); it is architecture-agnostic (apt installs for the native arch). It is a separate, explicit step so that a plain `make` never runs `sudo` on its own. `make` compiles the binary into `build/dynamic/<arch>/MethylExtractor`, where `<arch>` is auto-detected (`x64` on x86_64, `arm64` on aarch64).

The **HDF5 Zstd plugin** enables Zstd-compressed HDF5 output. Build it into the arch-correct project tree (`build/dynamic/<arch>/hdf5_zstd_plugin/`) with:

```bash
make plugin   # builds libH5Zzstd.so for this machine's architecture
```

`make plugin` prints the `HDF5_PLUGIN_PATH` to export if you want Zstd without a system install. The plugin is **not** copied to the system until you run `make install`. Without the plugin on `HDF5_PLUGIN_PATH`, HDF5 output falls back to gzip automatically.

### System install (optional)
```bash
make install
```

This copies the MethylExtractor binary to `/usr/local/bin/`, the HDF5 Zstd plugin to `/usr/local/hdf5/lib/plugin`, and appends `HDF5_PLUGIN_PATH` to your `~/.bashrc`. Before copying anything, **architecture is checked**: the binary and plugin must match the current machine (e.g. you cannot install an arm64 build on an x86_64 host). If there is a mismatch, `make install` aborts with an error—run `make` on the target machine first.

## Usage

### Basic Syntax
```bash
build/dynamic/x64/MethylExtractor [options] <input.bam> [output_directory] [reference.fa]
```
*(Note: Binary is located in `build/dynamic/x64/MethylExtractor` or similar, depending on architecture.)*

The output directory may be given positionally (2nd argument) or via `-o/--output-dir`; one of the two is required. The reference is optional and, if provided, overrides the one in `chrom_mapping.json`.

### Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `-h, --help` | Show help message and exit | N/A |
| `-t, --threads INT` | Number of worker threads | auto (CPU count) |
| `-q, --min-mapq INT` | Minimum mapping quality (MAPQ) | 30 |
| `-p, --min-phred INT` | Minimum base quality (Phred score) | 20 |
| `-c, --min-cov INT` | Minimum coverage threshold | 4 |
| `-C, --cap-cov INT` | Cap coverage to this value (0 = disabled) | 0 |
| `-G, --CHG` | Include CHG methylation contexts | Disabled |
| `-H, --CHH` | Include CHH methylation contexts | Disabled |
| `-m, --chrom-mapping FILE` | Chromosome mapping configuration file | `chrom_mapping.json` |
| `-z, --compression INT` | Compression level (0=none, 1-19=Zstd, gzip fallback) | 6 |
| `-k, --chunk-size INT` | HDF5 chunk size for I/O optimization | 1,000,000 |
| `-f, --output-format STR` | Output format: `hdf5`, `txt`, or `both` | `hdf5` |
| `-s, --split` | Split output by methylation context | Disabled |
| `-o, --output-dir DIR` | Output directory | N/A |

### Compression

HDF5 output is compressed in-place with the Zstd HDF5 filter when it is available
(see `make install`, which sets up `HDF5_PLUGIN_PATH`). If the Zstd filter is not
available at runtime, MethylExtractor falls back to gzip/deflate automatically.

```bash
build/dynamic/x64/MethylExtractor input.bam output_dir -z 9 -t 32
```

- `-z 0` disables compression (fastest write, largest files).
- `-z 1-19` selects the Zstd level (higher = smaller, slower). With the gzip
  fallback the level is capped at 9.
- Output files remain valid, self-describing HDF5 in all cases.

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

Each output file (`{chrom}-{context}.h5` or `{chrom}.h5`) gets a companion JSON sidecar with legacy top-level fields plus an extended QC block for [MethylPipeline](https://github.com/) post-extraction Pass/Fail evaluation.

**Legacy fields** (unchanged for backward compatibility):

```json
{
  "num_positions": 1345382,
  "total_methylated": 5084976,
  "total_unmethylated": 18665653,
  "avg_methylation_level": 0.214,
  "avg_coverage": 17.65
}
```

**Extended context QC** (`schema_name`: `methylextractor.context_qc`) adds:

- `metadata` — chromosome, context, filters (`min_mapq`, `min_phred`, `min_cov`, `cap_cov`)
- `sites` — `sites_in_reference`, `sites_passing_min_cov`, `fraction_sites_covered`
- `coverage` — mean/median/p10/p90 on all reference sites
- `methylation` — global level and `by_strand` (+/−)
- `read_filtering` — per-chromosome read drop counts and retention rate

### Extraction manifest (sample-level)

After all chromosomes complete, MethylExtractor writes `{sample_id}.extraction_manifest.json` in the output directory (`schema_name`: `methylextractor.extraction_manifest`). MethylPipeline should consume this file (and/or per-context JSON sidecars) to run post-extraction guardrails — Pass/Fail is **not** computed in MethylExtractor.

The manifest contains:

- `metadata` — sample ID, BAM path, reference, contexts extracted (CG always; CHG/CHH when `-G`/`-H` used), filter parameters
- `summary` — genome-wide weighted CpG coverage/methylation, CHG/CHH methylation when applicable, `cpg_fraction_sites_covered`
- `read_filtering` — aggregated read filter stats across chromosomes
- `per_chromosome` — per-context stats and read filtering for each chromosome

Example path: `/work/samples/1401-042825-50082/1401-042825-50082.extraction_manifest.json`

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
- **Overlapping read handling**: Coordinate-based mate clipping prevents double-counting in paired-end data (the left mate yields the overlap to the right mate so each reference position is counted once)
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
- **Slow processing**: Use SSD storage, increase `--threads`, and adjust chunk sizes. Lower the `-z` level (or use `-z 0`) to trade compression ratio for write speed.
- **Low methylation**: Verify bisulfite conversion and quality filters

### Performance Optimization
- Use SSD storage for input/output
- Adjust `--chunk-size` based on chromosome sizes
- Balance `--compression` level with I/O performance needs

## Documentation

For complete documentation, see `MethylExtractor_Documentation.html`.

## License

Copyright © Epimethyl Analytics. All rights reserved. This software is proprietary;
see [`LICENSE`](LICENSE) for terms.

MethylExtractor is a fork of [MethylDackel](https://github.com/dpryan79/MethylDackel),
which is distributed under the MIT License; portions derived from it remain under
those terms (see the third-party notice in [`LICENSE`](LICENSE)).