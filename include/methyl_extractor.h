#ifndef METHYL_EXTRACTOR_H
#define METHYL_EXTRACTOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>
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
#define DEFAULT_CHROM_PARALLEL 2
#define DEFAULT_MAX_RSS_GB 32
#define DEFAULT_BGZF_THREADS 2
#define DEFAULT_CHUNK_SIZE 1000000
#define DEFAULT_MIN_MAPQ 30
#define DEFAULT_MIN_PHRED 20
#define DEFAULT_MIN_COV 4
#define DEFAULT_CAP_COVERAGE 0
#define DEFAULT_TILE_SIZE 4
#define MIN_TILE_SIZE 2
#define MAX_TILE_SIZE 8
#define DEFAULT_FLAGS (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY)

#define READ_LEVEL_PATTERN_ENCODING "bitmask_msb_first"
#define READ_LEVEL_SCHEMA_VERSION "1.0.0"
#define READ_LEVEL_GROUP "read_level_patterns"

#define TNC_A 0
#define TNC_C 1
#define TNC_G 2
#define TNC_T 3
#define TNC_N 4

#define CONTEXT_CPG 1
#define CONTEXT_CHG 2
#define CONTEXT_CHH 3

#define ZSTD_FILTER 32015

#define MAX_CONTEXTS_PER_CHR 3
#define COV_HIST_BINS 256

#define EXTRACTION_CONTEXT_QC_SCHEMA "methylextractor.context_qc"
#define EXTRACTION_MANIFEST_SCHEMA "methylextractor.extraction_manifest"
#define EXTRACTION_TIMING_SCHEMA "methylextractor.timing"
#define EXTRACTION_CONTEXT_QC_VERSION "1.0.0"
#define EXTRACTION_MANIFEST_VERSION "1.0.0"
#define EXTRACTION_TIMING_VERSION "1.0.0"

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

typedef struct
{
    uint32_t *mC;
    uint32_t *uC;
    size_t size;
} PrivateCounts;

typedef struct
{
    uint64_t reads_seen;
    uint64_t reads_used;
    uint64_t reads_dropped_unmapped;
    uint64_t reads_dropped_secondary;
    uint64_t reads_dropped_qc_fail;
    uint64_t reads_dropped_duplicate;
    uint64_t reads_dropped_supplementary;
    uint64_t reads_dropped_low_mapq;
    uint64_t reads_dropped_multimap;
    uint64_t reads_dropped_no_strand;
    uint64_t bases_skipped_overlap_clip;
    uint64_t bases_skipped_low_phred;
    uint64_t bases_skipped_non_cytosine;
    uint64_t bases_counted_methylated;
    uint64_t bases_counted_unmethylated;
} FilterStats;

typedef struct
{
    size_t num_positions;
    uint64_t total_methylated, total_unmethylated;
    double avg_methylation_level, avg_coverage;
    size_t sites_in_reference;
    size_t sites_with_any_coverage;
    size_t sites_below_min_cov;
    size_t sites_capped;
    double avg_coverage_all_sites;
    double coverage_median;
    double coverage_p10;
    double coverage_p90;
    double methylation_plus;
    double methylation_minus;
    uint64_t methylated_plus;
    uint64_t unmethylated_plus;
    uint64_t methylated_minus;
    uint64_t unmethylated_minus;
} MethylStats;

typedef struct
{
    char chromosome[64];
    uint64_t site_enum_ms;
    uint64_t bam_scan_ms;
    uint64_t merge_ms;
    uint64_t hdf5_write_ms;
    uint64_t qc_json_ms;
    uint64_t total_ms;
    uint64_t rss_est_bytes;
} ChromosomeTiming;

typedef struct
{
    uint64_t total_ms;
    uint64_t bam_ms;
    uint64_t write_ms;
    int chrom_parallel;
    int region_threads;
    int max_rss_gb;
    int n_chromosomes;
} SampleTiming;

typedef struct
{
    const char *chromosome;
    const char *context;
    const char *output_path;
    int min_mapq;
    int min_phred;
    int min_cov;
    int cap_cov;
    ChromosomeTiming *timing;
} ExtractionMeta;

typedef struct
{
    char context[8];
    char output_path[1024];
    MethylStats stats;
} ContextExport;

typedef struct
{
    char chromosome[64];
    int n_contexts;
    ContextExport contexts[MAX_CONTEXTS_PER_CHR];
    FilterStats filter_stats;
} ChromosomeExport;

typedef struct
{
    const char *bam_file;
    const char *out_dir;
    const char *reference;
    int min_mapq;
    int min_phred;
    int min_cov;
    int cap_cov;
    int keep_chg;
    int keep_chh;
    int split_context_files;
    int read_level;
    int tile_size;
} SampleRunInfo;

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
    int num_threads;
    int read_level;
    int tile_size;
    int bgzf_threads;
    bam_hdr_t *hdr;
    ChromosomeExport *chr_export;
    ChromosomeTiming *timing;
} ThreadArg;

typedef struct
{
    int context;
    int tile_size;
    size_t n_cpg;
    uint32_t *cpg_pos;
    int *pos_to_cpg;
    size_t n_tiles;
    uint32_t *tile_start_pos;
    uint32_t *tile_cpg_positions;
} ReadLevelTiles;

typedef struct
{
    ReadLevelTiles tiles;
    uint32_t *hist;
    int active;
} ReadLevelContext;

typedef struct
{
    ReadLevelContext ctx[MAX_CONTEXTS_PER_CHR + 1];
    int n_active;
} ReadLevelChromData;

typedef struct
{
    ThreadArg base;
    PrivateCounts counts;
    hts_idx_t *idx;
    size_t start_site_idx;
    size_t end_site_idx;
    FilterStats filter_stats;
    ReadLevelChromData *rl;
    size_t *rl_tile_start;
    size_t *rl_tile_end;
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
    uint64_t bins[COV_HIST_BINS];
    uint64_t overflow;
    uint64_t total_sites;
} CoverageHistogram;

// utils.c
void log_time(const char *format, ...);
uint64_t now_ms(void);
uint64_t elapsed_ms(uint64_t start_ms);
int make_directory(const char *path);
void export_lock(void);
void export_unlock(void);
uint64_t estimate_chromosome_rss_bytes(uint32_t chr_len, int keep_chg, int keep_chh,
                                       int read_level, int split_context_files);
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
void merge_filter_stats(FilterStats *dst, const FilterStats *src);
void coverage_histogram_init(CoverageHistogram *hist);
void coverage_histogram_add(CoverageHistogram *hist, uint32_t coverage);
double coverage_histogram_percentile(const CoverageHistogram *hist, double percentile);
MethylStats calculate_statistics(MethylRecord *filtered_buffer, size_t n_records);
MethylStats analyze_buffer(MethylRecord *buffer, size_t n_records, int min_cov);
size_t flush_buffer(const char *filename, MethylRecord *buffer, size_t n_records,
                    int compression, int chunk_size, int append_mode,
                    int min_cov, int cap_cov, OutputFormat output_format,
                    const ExtractionMeta *meta, const FilterStats *filter_stats,
                    MethylStats *out_stats);
void cleanup_hdf5(void);

// extraction_export.c
int write_context_qc_json(const char *base_filename, const ExtractionMeta *meta,
                          const FilterStats *filter_stats, MethylStats stats);
int write_extraction_manifest(const SampleRunInfo *run,
                              const ChromosomeExport *chromosomes, int n_chromosomes);
int write_extraction_timing(const char *out_dir, const SampleTiming *sample,
                            const ChromosomeTiming *chromosomes, int n_chromosomes);
void path_basename(const char *path, char *out, size_t out_len);

// bam_processing.c
int getRealStrand(bam1_t *b);
void *process_region_direct(void *arg);
void process_chromosome(ThreadArg *targ);
void process_chromosomes_parallel(ThreadArg *args, int n_chroms, int chrom_parallel,
                                  uint64_t max_rss_bytes);
int load_chrom_mapping(const char *filename, ChromMapEntry **entries, int *n_entries, char **reference_file);

// read_level.c
void free_read_level_tiles(ReadLevelTiles *tiles);
int build_read_level_tiles(const char *chr_seq, uint32_t chr_len, int context,
                           int tile_size, ReadLevelTiles *out);
size_t read_level_tile_range_for_region(const ReadLevelTiles *tiles,
                                        uint32_t region_start,
                                        uint32_t region_end,
                                        size_t *out_start, size_t *out_end);
void read_level_accumulate(const ReadLevelTiles *tiles, uint32_t *hist,
                           size_t tile_start, size_t tile_end,
                           const char *chr_seq, uint8_t *seq, uint8_t *qual,
                           uint32_t *cigar, int n_cigar, int32_t read_start,
                           int32_t clip_from, int strand, int min_phred,
                           uint32_t chr_len);
int write_read_level_patterns_h5(const char *filename, const char *context,
                                 int tile_size, int min_tile_reads,
                                 int compression,
                                 const ReadLevelTiles *tiles,
                                 const uint32_t *hist);
void free_read_level_chrom_data(ReadLevelChromData *data);

#endif // METHYL_EXTRACTOR_H
