# EpiExtractor

EpiExtractor is a command-line tool for extracting methylation data from BAM files and writing it to HDF5 format.

## Project Structure

```
epiextractor/
├── src/
│   └── epiextractor.c
├── include/
│   └── (header files for tshlib and hdf5 if needed)
├── lib/
│   └── (static/dynamic libs for tshlib and hdf5 if custom)
├── build/
│   ├── static/
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
   This will create a statically linked binary at `build/static/epiextractor`.

2. **Dynamic Build**:
   ```bash
   make dynamic
   ```
   This will create a dynamically linked binary at `build/dynamic/epiextractor`.

3. **Containerized Build**:
   ```bash
   docker build -t epiextractor:latest .
   ```
   This will build a Docker image named `epiextractor:latest`.

## Usage

Run the tool with the following command:

```bash
./build/dynamic/epiextractor --hdf5-output-dir output --hdf5-compression 6 --hdf5-chunk-size 100000 --threads 8 reference.fa input.bam
```

For help and additional options:

```bash
./build/dynamic/epiextractor --help
```

dizada@EM-DIZADA:~$ samtools faidx /home/dizada/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
dizada@EM-DIZADA:~$ cut -f 1 /home/dizada/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai
1
2
3
4
5
Mt
Pt
dizada@EM-DIZADA:~$ samtools view -H /home/dizada/A14.clara_parabrics.duplicates_marked.bam | grep '^@SQ' | cut -f 2 | cut -d':' -f 2
1
2
3
4
5
Mt
Pt

## License

[Add your license information here] 