# MethylExtractor

A tool for extracting methylation data from BAM files.

## Installation

```bash
git clone https://github.com/yourusername/MethylExtractor.git
cd MethylExtractor
make
```

## Usage

```bash
./MethylExtractor [options] <input.bam> [reference.fa]
```

### Options

- `-h, --help`: Show help message
- `-t, --threads INT`: Number of threads [1]
- `-q, --min-mapq INT`: Minimum mapping quality [0]
- `-p, --min-phred INT`: Minimum base quality [0]
- `-c, --cap-cov INT`: Cap coverage to this value [0]
- `-G, --CHG`: Include CHG context
- `-H, --CHH`: Include CHH context
- `-m, --chrom-mapping FILE`: Chromosome mapping file
- `-z, --compression INT`: HDF5 compression level [6]
- `-k, --chunk-size INT`: HDF5 chunk size [1000]
- `-f, --output-format STR`: Output format (hdf5, txt, both) [hdf5]
- `-s, --split`: Split output by context
- `-o, --output-dir DIR`: Output directory

### Output Formats

- `hdf5`: HDF5 format (default)
- `txt`: Text format
- `both`: Both HDF5 and text formats

### Examples

Basic usage:
```bash
./MethylExtractor input.bam reference.fa
```

With options:
```bash
./MethylExtractor -t 4 -q 20 -p 20 -G -H -f both -o output_dir input.bam reference.fa
```

## Output Files

The tool generates the following output files:

- `output.h5`: HDF5 file containing methylation data
- `output.txt`: Text file containing methylation data (if text format is selected)
- `output.CG.txt`: Text file containing CG context data (if split by context)
- `output.CHG.txt`: Text file containing CHG context data (if split by context)
- `output.CHH.txt`: Text file containing CHH context data (if split by context)

## Dependencies

- HTSlib
- HDF5
- khash
- pthread

## License

This project is licensed under the MIT License - see the LICENSE file for details. 