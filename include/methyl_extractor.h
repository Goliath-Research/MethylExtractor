#ifndef METHYL_EXTRACTOR_H
#define METHYL_EXTRACTOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/stat.h>
#include <ctype.h>
#include <errno.h>
#include <hdf5/serial/hdf5.h>
#include <pthread.h>
#include <sys/sysinfo.h>
#include <time.h>
#include "cjson/cJSON.h"
#include <htslib/hts_defs.h>

#define DEFAULT_MAX_CHR 24
#define DEFAULT_HDF5_COMPRESSION 6
#define DEFAULT_HDF5_CHUNK_SIZE 1000000
#define DEFAULT_THREADS 16
#define DEFAULT_CHUNK_SIZE 1000000
#define DEFAULT_MIN_MAPQ 30
#define DEFAULT_MIN_PHRED 20
#define DEFAULT_MIN_COV 4
#define DEFAULT_CAP_COVERAGE 0
#define DEFAULT_FLAGS (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY)

#define TNC_A 0
#define TNC_C 1
#define TNC_G 2
#define TNC_T 3
#define TNC_N 4

#define CONTEXT_CPG 1
#define CONTEXT_CHG 2
#define CONTEXT_CHH 3

#define ZSTD_FILTER 32015

typedef enum
{
    OUTPUT_NONE = 0,
    OUTPUT_HDF5 = 1,
    OUTPUT_TXT = 2,
    OUTPUT_BOTH = 3
} OutputFormat;

typedef struct
{
    unsigned tnc : 5;
    unsigned context : 2;
    unsigned strand : 1;
} tnc_bitfield_t;

typedef struct
{
    uint32_t pos;
    uint16_t mC;
    uint16_t uC;
    tnc_bitfield_t tnc;
    uint8_t _pad[1];
} MethylRecord;

// Private per-thread counters (dense genomic arrays)
typedef struct
{
    uint32_t *mC;
    uint32_t *uC;
    size_t size;
} PrivateCounts;

typedef struct
{
    const char *bam_file;
    const char *out_dir;
    int tid;
    char *chr;
    uint32_t chr_len;
    char *chr_seq;
    int min_mapq;
    int min_phred;
    int min_cov;
    int cap_cov;
    int keep_chg;
    int keep_chh;
    int compression;
    int hdf5_chunk_size;
    uint32_t chunk_size;
    uint32_t start_pos;
    uint32_t end_pos;
    OutputFormat output_format;
    int split_context_files;
} ThreadArg;

// Extended for private counting
typedef struct
{
    ThreadArg base;
    PrivateCounts counts;
    size_t start_site_idx;
    size_t end_site_idx;
} RegionArg;

typedef struct
{
    char fasta[64];
    char bam[64];
    char name[64];
    int extract;
} ChromMapEntry;

typedef struct
{
    size_t num_positions;
    uint64_t total_methylated, total_unmethylated;
    double avg_methylation_level, avg_coverage;
} MethylStats;

// Function Prototypes

// utils.c
void log_time(const char *format, ...);
int make_directory(const char *path);
uint8_t encode_nucleotide(char n);
char decode_nucleotide(uint8_t n);
void decode_trinucleotide(uint8_t tnc, char *trinucl);
const char *get_context_string(int context);
uint8_t encode_trinucleotide_context(const char *chr_seq, int pos, int chr_len, char strand);
int isCpG(const char *seq, int pos, int seqlen);
int isCHG(const char *seq, int pos, int seqlen);
int isCHH(const char *seq, int pos, int seqlen);
int get_context(const char *chr_seq, int chr_len, int pos, int8_t *strand_ctx, uint8_t *tnc, int keep_chg, int keep_chh);
size_t count_methylation_sites(const char *chr_seq, uint32_t chr_len, int keep_chg, int keep_chh);
void initialize_buffer(MethylRecord *buffer, size_t site_count, const char *chr_seq, uint32_t chr_len, int keep_chg, int keep_chh);

// output_formats.c
MethylStats calculate_statistics(MethylRecord *filtered_buffer, size_t n_records);
int write_statistics_json(const char *base_filename, MethylStats stats);
size_t flush_buffer(const char *filename, MethylRecord *buffer, size_t n_records,
                    int compression, int chunk_size, int append_mode,
                    int min_cov, int cap_cov,
                    OutputFormat output_format);
void cleanup_hdf5(void);

// bam_processing.c
int getRealStrand(bam1_t *b);
void *process_region_direct(void *arg);
void process_chromosome(ThreadArg *targ);
int load_chrom_mapping(const char *filename, ChromMapEntry **entries, int *n_entries, char **reference_file);

#endif // METHYL_EXTRACTOR_H
