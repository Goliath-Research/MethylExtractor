# MethylExtractor

MethylExtractor is a command-line tool for extracting methylation data from BAM files and writing it to HDF5 or text format.

## Project Structure

```
MethylExtractor/
├── src/
│   └── MethylExtractor.c
├── include/
│   └── (header files for tshlib and hdf5 if needed)
├── lib/
│   └── (static/dynamic libs for tshlib and hdf5 if custom)
├── build/
│   ├── static/
│   ├── debug/
│   ├── dynamic/
│   └── docker/
├── Dockerfile
├── Makefile
└── README.md
```

## Building the Project

### Prerequisites

- GCC compiler
- Make
- HDF5 library and development headers
- HTSlib library and development headers

### Build Instructions

1. **Static Build**:
   ```bash
   make static
   ```
   This will create a statically linked binary at `build/static/MethylExtractor`.

2. **Dynamic Build**:
   ```bash
   make dynamic
   ```
   This will create a dynamically linked binary at `build/dynamic/MethylExtractor`.

3. **Containerized Build**:
   ```bash
   docker build -t MethylExtractor:latest .
   ```
   This will build a Docker image named `MethylExtractor:latest`.

## Usage

Run the tool with the following command:

```bash
./build/dynamic/MethylExtractor [options] <input.bam> [ref.fa]
```

### Required Arguments:
- `<input.bam>`: Input BAM file containing methylation data

### Optional Arguments:
- `[ref.fa]`: Reference FASTA file (optional, will override reference in chrom_mapping.json if provided)

### Options:
```
  -h, --help                Show help message
  -q, --min-mapq INT        Minimum mapping quality [30]
  -p, --min-phred INT       Minimum base quality [20]
  -c, --cap-cov INT         Cap coverage (optional, default: 0)
  -G, --CHG                 Process CHG context
  -H, --CHH                 Process CHH context
  -m, --chrom-mapping FILE  Chromosome mapping file [chrom_mapping.json]
  -z, --compression INT     HDF5 compression level [6]
  -k, --chunk-size INT      HDF5 chunk size [1000000]
  -f, --output-format STR   Output format (hdf5, txt, both) [hdf5]
  -s, --split-context-files Split output by context
  -o, --output-dir DIR      Output directory (required)
```

### Example Usage:

1. Basic usage with default options:
```bash
./build/dynamic/MethylExtractor -o output_dir input.bam
```

2. Process all contexts with custom quality thresholds:
```bash
./build/dynamic/MethylExtractor -G -H -q 20 -p 15 -o output_dir input.bam
```

3. Output in text format with split context files:
```bash
./build/dynamic/MethylExtractor -f txt -s -o output_dir input.bam
```

4. Use custom reference file:
```bash
./build/dynamic/MethylExtractor -o output_dir input.bam ref.fa
```

## License

[Add your license information here] 