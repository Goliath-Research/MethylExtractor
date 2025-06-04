# MethylExtractor

MethylExtractor is a command-line tool for extracting methylation data from BAM files and writing it to HDF5 format.

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
./build/dynamic/MethylExtractor --hdf5-output-dir output --hdf5-compression 6 --hdf5-chunk-size 100000 --threads 8 reference.fa input.bam
```

For help and additional options:

```bash
./build/dynamic/MethylExtractor --help
```
sudo apt install samtools
sudo apt install parallel
find /home/ubuntu/Work/HRA006113 -type f -name "*.bam" | parallel samtools index {}

## License

[Add your license information here] 