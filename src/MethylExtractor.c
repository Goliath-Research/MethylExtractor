#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/stat.h>
#include <ctype.h>
#include <hdf5/serial/hdf5.h>
#include <pthread.h>
#include <sys/sysinfo.h>
#include <time.h>
#include "cjson/cJSON.h"

// Memory size constants
#define KB (1024ULL)
#define MB (KB * 1024ULL)
#define GB (MB * 1024ULL)
#define TB (GB * 1024ULL)

// Memory allocation constants
#define MIN_BUFFER_SIZE (1ULL * GB)
#define MAX_BUFFER_PERCENT 0.5
#define OPTIMAL_BUFFER_PERCENT 0.25
#define MIN_BUFFER_PERCENT 0.001
#define SAMPLE_READ_COUNT 100000

// Coverage thresholds
#define HIGH_COVERAGE_THRESHOLD 100
#define MEDIUM_COVERAGE_THRESHOLD 50

// Chunk size constants
#define HIGH_COVERAGE_CHUNK (1ULL * GB)      // 1GB for high coverage
#define MEDIUM_COVERAGE_CHUNK (2ULL * GB)    // 2GB for medium coverage
#define NORMAL_COVERAGE_CHUNK (4ULL * GB)    // 4GB for normal coverage

// Add processing chunk size constant
#define PROCESSING_CHUNK_SIZE (1ULL * GB)  // Process in 1GB chunks

// Optional CUDA support
#ifdef USE_CUDA
#include <cuda_runtime.h>
#include <cuda.h>
#endif

// System resource detection
static long get_total_memory() 
{
    struct sysinfo si;
    if (sysinfo(&si) == 0)
        return si.totalram * si.mem_unit;
    return 432L * 1024 * 1024 * 1024; // Default to 432GB if sysinfo fails
}

static int get_optimal_thread_count() 
{
    int cpu_count = sysconf(_SC_NPROCESSORS_ONLN);
    return cpu_count > 0 ? cpu_count : 64; // Default to 64 if detection fails
}

static size_t get_optimal_chunk_size() 
{
    long total_mem = get_total_memory();
    // Use 10% of total memory for chunk size, but cap at 32GB
    size_t chunk = (size_t)(total_mem * 0.1);
    return chunk > 32L*1024*1024*1024 ? 32L*1024*1024*1024 : chunk;
}

// Dynamic defaults based on system resources
#define DEFAULT_MAX_CHR 24
#define DEFAULT_HDF5_COMPRESSION 6
#define DEFAULT_HDF5_CHUNK_SIZE get_optimal_chunk_size()
#define DEFAULT_THREADS get_optimal_thread_count()
#define DEFAULT_CHUNK_SIZE get_optimal_chunk_size()
#define MAX_REGIONS_PER_CHR (get_optimal_thread_count() * 4)  // 4 regions per CPU

// Quality thresholds (these can stay fixed)
#define DEFAULT_MIN_MAPQ 30
#define DEFAULT_MIN_PHRED 20
#define DEFAULT_MIN_COV 4
#define DEFAULT_CAP_COVERAGE 1
#define DEFAULT_MIN_METH 0
#define DEFAULT_MAX_METH 100
#define DEFAULT_FLAGS (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY)
#define TNC_A 0
#define TNC_C 1
#define TNC_G 2
#define TNC_T 3
#define TNC_N 4
#define CONTEXT_CPG 1
#define CONTEXT_CHG 2
#define CONTEXT_CHH 3
#define STRAND_MASK 0x80
#define CONTEXT_MASK 0x03
#define MAX_CHR_NAME 2
#define TNC_MASK 0x1F
#define STRAND_BITS_MASK 0x60
#define SIGN_MASK 0x80

// Hash table for position-to-buffer-index mapping
KHASH_MAP_INIT_INT64(pos, size_t)
KHASH_SET_INIT_STR(str)

#define ZSTD_FILTER 32015  // Zstandard filter ID

// Output format types
typedef enum 
{
    OUTPUT_NONE = 0,
    OUTPUT_HDF5 = 1,
    OUTPUT_TXT = 2,
    OUTPUT_BOTH = 3
} OutputFormat;

typedef struct 
{
    unsigned tnc     : 5;
    unsigned context : 2;
    unsigned strand  : 1;
} tnc_bitfield_t;

typedef struct
{
    uint32_t pos;
    uint16_t mC;
    uint16_t uC;
    tnc_bitfield_t tnc; // Now stores strand + context + TNC as bitfield
    uint8_t _pad[1]; // explicit padding for alignment
} MethylRecord;

typedef struct
{
    char *chr;
    uint32_t start_pos;
    uint32_t end_pos;
    int tid;
    uint32_t chr_len;
    const char *bam_file;
    const char *out_dir;
    int min_phred;
    int min_mapq;
    int min_cov;
    int cap_cov;
    int min_meth;
    int max_meth;
    int compression;
    int chunk_size;
    int append_mode;
    OutputFormat output_format;
    int split_context_files;
    MethylRecord **context_buffers;  // Array of context buffers
    khash_t(pos) **context_pos_maps; // Array of position maps for each context
    samFile *fp;
    bam_hdr_t *header;
    int use_gpu;
    char *chr_seq;
    int keep_chg;
    int keep_chh;
    int hdf5_compression;
    int hdf5_chunk_size;
    pthread_mutex_t *buffer_mutex;
} ThreadArg;

typedef struct 
{
    samFile *in;
    hts_itr_t *iter;
    bam_hdr_t *hdr;
    ThreadArg *targ;
} mplp_data_t;

typedef struct 
{
    char fasta[64];
    char bam[64];
    char name[64];
    int extract;
} ChromMapEntry;

// Add new structure for memory requirements
typedef struct 
{
    size_t min_buffer_size;    // Minimum buffer size needed
    size_t optimal_buffer_size; // Optimal buffer size based on data
    size_t max_buffer_size;    // Maximum buffer size allowed
    size_t chunk_size;         // Processing chunk size
    int region_count;          // Number of regions to process in parallel
} MemoryRequirements;

// Add new structure for BAM statistics
typedef struct {
    double avg_coverage;
    size_t total_reads;
    size_t total_bases;
    size_t max_read_length;
    int has_coverage_info;
} BamStats;

// Add global memory requirements
static MemoryRequirements global_mem_req = {0};

static void log_time(const char *format, ...) 
{
    time_t rawtime;
    struct tm *timeinfo;
    char time_str[64];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(time_str, sizeof(time_str), "%Y-%m-%d %H:%M:%S", timeinfo);
    
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[%s] ", time_str);
    vfprintf(stderr, format, args);
    va_end(args);
}

int make_directory(const char *path)
{
    struct stat st = {0};
    if (stat(path, &st) == -1)
        return mkdir(path, 0700);
    return 0;
}

static inline uint8_t encode_nucleotide(char n)
{
    switch (toupper(n))
    {
    case 'A':
        return TNC_A;
    case 'C':
        return TNC_C;
    case 'G':
        return TNC_G;
    case 'T':
        return TNC_T;
    default:
        return TNC_N;
    }
}

static inline char decode_nucleotide(uint8_t n)
{
    switch (n)
    {
    case TNC_A:
        return 'A';
    case TNC_C:
        return 'C';
    case TNC_G:
        return 'G';
    case TNC_T:
        return 'T';
    default:
        return 'N';
    }
}

static inline void decode_trinucleotide(uint8_t tnc, char *trinucl)
{
    // Decode MethylDackel-style TNC index (0-24) to first and second bases (middle base not stored)
    uint8_t n2 = (tnc / 5) % 4; // Middle base index
    uint8_t n3 = tnc % 5;       // Last base index
    trinucl[0] = 'C';           // Central base is always C in MethylDackel naming
    trinucl[1] = decode_nucleotide(n2);
    trinucl[2] = decode_nucleotide(n3);
    trinucl[3] = '\0';
}

static inline const char *get_context_string(int context)
{
    if (context == CONTEXT_CPG)
        return "CG";
    else if (context == CONTEXT_CHG)
        return "CHG";
    else if (context == CONTEXT_CHH)
        return "CHH";
    return "???";
}

static inline uint8_t encode_trinucleotide_context(const char *chr_seq, int pos, int chr_len, char strand)
{
    // MethylDackel-style: encode based on surrounding bases considering strand direction
    int direction = (strand == '+') ? 1 : -1;
    uint8_t rv = 0;
    char base;

    // Last base: column
    if ((direction > 0 && pos + 2 >= chr_len) || (direction < 0 && pos <= 1))
        rv = 4;
    else
    {
        base = chr_seq[pos + 2 * direction];
        if (direction < 0)
            base = (base == 'A') ? 'T'  : (base == 'T') ? 'A'
                                        : (base == 'C') ? 'G'
                                        : (base == 'G') ? 'C'
                                        : 'N';
        switch (toupper(base))
        {
        case 'A':
            rv = 0;
            break;
        case 'C':
            rv = 1;
            break;
        case 'G':
            rv = 2;
            break;
        case 'T':
            rv = 3;
            break;
        default:
            rv = 4;
            break;
        }
    }

    // Middle base
    if ((direction > 0 && pos + 1 >= chr_len) || (direction < 0 && pos == 0))
        rv += 20;
    else
    {
        base = chr_seq[pos + direction];
        if (direction < 0)
            base = (base == 'A') ? 'T'  : (base == 'T') ? 'A'
                                        : (base == 'C') ? 'G'
                                        : (base == 'G') ? 'C'
                                        : 'N';
        switch (toupper(base))
        {
        case 'A':
            rv += 0;
            break;
        case 'C':
            rv += 5;
            break;
        case 'G':
            rv += 10;
            break;
        case 'T':
            rv += 15;
            break;
        default:
            rv += 20;
            break;
        }
    }
    return rv; // 0-24
}

static inline int8_t encode_strand_context(char strand, int context)
{
    return strand == '-' ? -context : context;
}

static inline int isCpG(const char *seq, int pos, int seqlen)
{
    if (pos >= seqlen)
        return 0;
    if (toupper(*(seq + pos)) == 'C')
    {
        if (pos + 1 == seqlen)
            return 0;
        if (toupper(*(seq + pos + 1)) == 'G')
            return 1;
        return 0;
    }
    else if (toupper(*(seq + pos)) == 'G')
    {
        if (pos == 0)
            return 0;
        if (toupper(*(seq + pos - 1)) == 'C')
            return -1;
        return 0;
    }
    return 0;
}

static inline int isCHG(const char *seq, int pos, int seqlen)
{
    if (pos >= seqlen)
        return 0;
    if (toupper(*(seq + pos)) == 'C')
    {
        if (pos + 2 >= seqlen)
            return 0;
        if (toupper(*(seq + pos + 2)) == 'G')
            return 1;
        return 0;
    }
    else if (toupper(*(seq + pos)) == 'G')
    {
        if (pos <= 1)
            return 0;
        if (toupper(*(seq + pos - 2)) == 'C')
            return -1;
        return 0;
    }
    return 0;
}

static inline int isCHH(const char *seq, int pos, int seqlen)
{
    if (pos >= seqlen)
        return 0;
    if (toupper(*(seq + pos)) == 'C')
        return 1;
    else if (toupper(*(seq + pos)) == 'G')
        return -1;
    return 0;
}

int get_context(const char *chr_seq, int chr_len, int pos, int8_t *strand_ctx, uint8_t *tnc, int keep_chg, int keep_chh)
{
    int context = 0;
    if (isCpG(chr_seq, pos, chr_len))
        context = CONTEXT_CPG;
    else if (isCHG(chr_seq, pos, chr_len))
        context = CONTEXT_CHG;
    else if (isCHH(chr_seq, pos, chr_len))
        context = CONTEXT_CHH;

    char c2 = toupper(chr_seq[pos]);
    char strand = (c2 == 'C') ? '+' : '-';
    *tnc = encode_trinucleotide_context(chr_seq, pos, chr_len, strand);
    if (c2 == 'C')
        *strand_ctx = encode_strand_context('+', context);
    else if (c2 == 'G')
        *strand_ctx = encode_strand_context('-', context);
    else
        *strand_ctx = 0;

    return context;
}

size_t count_methylation_sites(const char *chr_seq, uint32_t chr_len, int keep_chg, int keep_chh)
{
    size_t count = 0;
    for (uint32_t pos = 0; pos < chr_len; pos++)
        if (isCpG((char *)chr_seq, pos, chr_len) ||
            (keep_chg && isCHG((char *)chr_seq, pos, chr_len)) ||
            (keep_chh && isCHH((char *)chr_seq, pos, chr_len)))
            count++;
    return count;
}

void initialize_buffer(
    MethylRecord *buffer, 
    size_t site_count, 
    const char *chr_seq, 
    uint32_t chr_len, 
    int keep_chg, 
    int keep_chh
)
{
    size_t idx = 0;
    for (uint32_t pos = 0; pos < chr_len && idx < site_count; pos++)
    {
        int8_t strand_ctx;
        uint8_t tnc_val;
        int ctx = get_context(chr_seq, chr_len, pos, &strand_ctx, &tnc_val, keep_chg, keep_chh);
        if (ctx)
        {
            char strand = (strand_ctx > 0) ? '+' : '-';
            buffer[idx].pos = pos + 1;
            buffer[idx].mC = 0;
            buffer[idx].uC = 0;
            buffer[idx].tnc.tnc = tnc_val;
            buffer[idx].tnc.context = ctx;
            buffer[idx].tnc.strand = (strand == '+') ? 0 : 1;
            idx++;
        }
    }
}

// CPU version of methylation processing
static void process_methylation_cpu(
    const uint8_t *seq,
    const uint8_t *qual,
    const int *strands,
    const uint32_t *positions,
    MethylRecord **context_buffers,
    khash_t(pos) **context_pos_maps,
    const int min_phred,
    const int min_mapq,
    const size_t n_records,
    const uint32_t start_pos
) 
{
    // Pre-compute position keys for faster lookup
    uint64_t *keys = malloc(n_records * sizeof(uint64_t));
    if (!keys) return;
    
    for (size_t i = 0; i < n_records; i++) 
    {
        keys[i] = ((uint64_t)positions[i] << 32) | (positions[i] - 1);
    }

    // Process records in batches for better cache utilization
    const size_t BATCH_SIZE = 1024;
    for (size_t batch_start = 0; batch_start < n_records; batch_start += BATCH_SIZE) 
    {
        size_t batch_end = (batch_start + BATCH_SIZE < n_records) ? batch_start + BATCH_SIZE : n_records;
        
        for (size_t i = batch_start; i < batch_end; i++) 
        {
            if (qual[i] < min_phred) continue;

            // Try each context in order of likelihood (CPG most common)
            for (int ctx = 1; ctx <= 3; ctx++) 
            {
                if (!context_buffers[ctx] || !context_pos_maps[ctx]) continue;
                
                khint_t iter = kh_get(pos, context_pos_maps[ctx], keys[i]);
                if (iter != kh_end(context_pos_maps[ctx])) 
                {
                    size_t idx = kh_val(context_pos_maps[ctx], iter);
                    if ((seq[i] == 'C' && (strands[i] == 1 || strands[i] == 3)) ||
                        (seq[i] == 'G' && (strands[i] == 2 || strands[i] == 4)))
                        context_buffers[ctx][idx].mC++;
                    else if ((seq[i] == 'T' && (strands[i] == 1 || strands[i] == 3)) ||
                             (seq[i] == 'A' && (strands[i] == 2 || strands[i] == 4)))
                        context_buffers[ctx][idx].uC++;
                    break;  // Found the correct context, no need to check others
                }
            }
        }
    }
    
    free(keys);
}

// CPU version of filtering and coverage calculation
static void filter_and_calculate_coverage_cpu(
    const MethylRecord *input_buffer,
    MethylRecord *output_buffer,
    const int min_cov,
    const int cap_cov,
    const double avg_cov,
    const int min_meth,
    const int max_meth,
    size_t *n_filtered
) 
{
    size_t j = 0;
    for (size_t i = 0; i < *n_filtered; i++) 
    {
        int total = input_buffer[i].mC + input_buffer[i].uC;
        if (total >= min_cov) 
        {
            double meth_level = total > 0 ? 100.0 * ((double)input_buffer[i].mC / total) : 0.0;
            if (meth_level >= min_meth && meth_level <= max_meth) 
            {
                if (cap_cov && total > avg_cov) 
                {
                    double prop = (double)input_buffer[i].mC / total;
                    output_buffer[j].mC = (uint16_t)round((avg_cov * prop));
                    output_buffer[j].uC = (uint16_t)(avg_cov - output_buffer[j].mC);
                }
                else
                    output_buffer[j] = input_buffer[i];
                j++;
            }
        }
    }
    *n_filtered = j;
}

#ifdef USE_CUDA
// Initialize CUDA device
static int init_cuda_device() 
{
    int device_count;
    cudaError_t error = cudaGetDeviceCount(&device_count);
    if (error != cudaSuccess || device_count == 0) 
    {
        fprintf(stderr, "No CUDA-capable device found\n");
        return -1;
    }

    // Select the first device
    error = cudaSetDevice(0);
    if (error != cudaSuccess) 
    {
        fprintf(stderr, "Failed to set CUDA device\n");
        return -1;
    }

    // Get device properties
    cudaDeviceProp prop;
    error = cudaGetDeviceProperties(&prop, 0);
    if (error != cudaSuccess) 
    {
        fprintf(stderr, "Failed to get device properties\n");
        return -1;
    }

    // Set device flags for better performance
    cudaSetDeviceFlags(cudaDeviceMapHost | cudaDeviceScheduleAuto);

    // Print device information
    fprintf(stderr, "Using GPU: %s\n", prop.name);
    fprintf(stderr, "Total GPU memory: %.2f GB\n", prop.totalGlobalMem / (1024.0 * 1024.0 * 1024.0));
    fprintf(stderr, "Number of multiprocessors: %d\n", prop.multiProcessorCount);
    fprintf(stderr, "Max threads per block: %d\n", prop.maxThreadsPerBlock);

    return 0;
}

// Calculate optimal block size for kernel
static int get_optimal_block_size() 
{
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    
    // Use 256 threads per block as default, but adjust based on device capabilities
    int block_size = 256;
    if (prop.maxThreadsPerBlock < block_size) 
        block_size = prop.maxThreadsPerBlock;
    
    // Ensure block size is a multiple of 32 (warp size)
    block_size = (block_size / 32) * 32;
    
    return block_size;
}

// CUDA kernel for processing methylation data
__global__ void process_methylation_kernel(
    const uint8_t *seq,
    const uint8_t *qual,
    const int *strands,
    const uint32_t *positions,
    MethylRecord *buffer,
    const int min_phred,
    const int min_mapq,
    const size_t n_records
) 
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_records) 
        return;

    // Process methylation data in parallel
    if (qual[idx] >= min_phred) 
    {
        int base = seq[idx];
        int strand = strands[idx];
        uint32_t pos = positions[idx];

        // Apply methylation logic
        if ((base == 'C' && (strand == 1 || strand == 3)) ||
            (base == 'G' && (strand == 2 || strand == 4)))
            atomicAdd(&buffer[pos].mC, 1);
        else if ((base == 'T' && (strand == 1 || strand == 3)) ||
                 (base == 'A' && (strand == 2 || strand == 4)))
            atomicAdd(&buffer[pos].uC, 1);
    }
}

// CUDA kernel for filtering and coverage calculation
__global__ void filter_and_calculate_coverage_kernel(
    const MethylRecord *input_buffer,
    MethylRecord *output_buffer,
    const int min_cov,
    const int cap_cov,
    const double avg_cov,
    const int min_meth,
    const int max_meth,
    size_t *n_filtered
) 
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= *n_filtered) 
        return;

    int total = input_buffer[idx].mC + input_buffer[idx].uC;
    if (total >= min_cov) 
    {
        double meth_level = total > 0 ? 100.0 * ((double)input_buffer[idx].mC / total) : 0.0;
        if (meth_level >= min_meth && meth_level <= max_meth) 
        {
            if (cap_cov && total > avg_cov) 
            {
                double prop = (double)input_buffer[idx].mC / total;
                output_buffer[idx].mC = (uint16_t)round((avg_cov * prop));
                output_buffer[idx].uC = (uint16_t)(avg_cov - output_buffer[idx].mC);
            }
            else
                output_buffer[idx] = input_buffer[idx];
        }
    }
}
#endif

// Extract coverage information from BAM header
static BamStats extract_bam_stats(bam_hdr_t *header, const char *chr_name) 
{
    BamStats stats = {0};
    stats.has_coverage_info = 0;
    
    // Check for coverage information in @CO lines
    char *coverage_str = NULL;
    char *coverage_key = "coverage=";
    char *coverage_end = NULL;
    
    for (int i = 0; i < header->n_targets; i++) 
    {
        if (strcmp(header->target_name[i], chr_name) == 0) 
        {
            // Look for coverage information in comments
            char *comment = header->text;
            while ((comment = strstr(comment, "@CO")) != NULL) 
            {
                if ((coverage_str = strstr(comment, coverage_key)) != NULL) 
                {
                    coverage_str += strlen(coverage_key);
                    stats.avg_coverage = strtod(coverage_str, &coverage_end);
                    if (coverage_end != coverage_str) 
                    {
                        stats.has_coverage_info = 1;
                        break;
                    }
                }
                comment += 3; // Move past "@CO"
            }
            break;
        }
    }
    
    return stats;
}

// Update calculate_memory_requirements to use header information
static MemoryRequirements calculate_memory_requirements(
    const char *bam_file,
    const char *chr_name,
    uint32_t chr_len,
    int min_mapq,
    int min_phred
) 
{
    MemoryRequirements req = {0};
    long total_mem = get_total_memory();
    int cpu_count = get_optimal_thread_count();
    
    fprintf(stderr, "Memory calculation for %s (length: %u)\n", chr_name, chr_len);
    fprintf(stderr, "Total system memory: %.2f GB\n", (double)total_mem / GB);
    
    // Set more aggressive defaults
    req.min_buffer_size = 1ULL * GB;  // 1GB minimum
    req.optimal_buffer_size = total_mem / 4;  // Use 25% of total memory
    req.max_buffer_size = total_mem / 2;  // Allow up to 50% of total memory
    req.chunk_size = 1ULL * GB;  // Default to 1GB chunks
    req.region_count = cpu_count * 4;  // More regions for better parallelization
    
    // Open BAM file to get statistics
    samFile *in = sam_open(bam_file, "r");
    if (!in) 
    {
        fprintf(stderr, "Failed to open BAM file for memory calculation\n");
        return req;
    }
    
    bam_hdr_t *header = sam_hdr_read(in);
    if (!header) 
    {
        fprintf(stderr, "Failed to read BAM header\n");
        sam_close(in);
        return req;
    }
    
    // Try to get coverage information from header first
    BamStats stats = extract_bam_stats(header, chr_name);
    double avg_coverage = 0.0;
    
    if (!stats.has_coverage_info) 
    {
        // If no coverage info in header, use a conservative estimate
        avg_coverage = 30.0;  // Assume 30x coverage if unknown
        fprintf(stderr, "No coverage info in header, using conservative estimate of 30x\n");
    }
    else
    {
        avg_coverage = stats.avg_coverage;
        fprintf(stderr, "Using coverage information from BAM header: %.2fx\n", avg_coverage);
    }
    
    // Calculate memory requirements based on coverage, but with hard limits
    size_t estimated_bases = (size_t)(chr_len * avg_coverage);
    fprintf(stderr, "Estimated bases to process: %zu\n", estimated_bases);
    
    // Calculate optimal buffer size with hard limits
    req.optimal_buffer_size = estimated_bases;
    if (req.optimal_buffer_size > req.max_buffer_size)
        req.optimal_buffer_size = req.max_buffer_size;
    if (req.optimal_buffer_size < req.min_buffer_size)
        req.optimal_buffer_size = req.min_buffer_size;
    
    // Adjust chunk size based on coverage and available memory
    if (avg_coverage > HIGH_COVERAGE_THRESHOLD)
        req.chunk_size = HIGH_COVERAGE_CHUNK;  // 1GB for high coverage
    else if (avg_coverage > MEDIUM_COVERAGE_THRESHOLD)
        req.chunk_size = MEDIUM_COVERAGE_CHUNK;  // 2GB for medium coverage
    else
        req.chunk_size = NORMAL_COVERAGE_CHUNK;  // 4GB for normal coverage
    
    // Ensure chunk size doesn't exceed optimal buffer size
    if (req.chunk_size > req.optimal_buffer_size)
        req.chunk_size = req.optimal_buffer_size;
    
    // Ensure chunk size is at least 1GB
    if (req.chunk_size < 1ULL * GB)
        req.chunk_size = 1ULL * GB;
    
    fprintf(stderr, "Final memory settings:\n");
    fprintf(stderr, "  Min buffer size: %.2f GB\n", (double)req.min_buffer_size / GB);
    fprintf(stderr, "  Optimal buffer size: %.2f GB\n", (double)req.optimal_buffer_size / GB);
    fprintf(stderr, "  Max buffer size: %.2f GB\n", (double)req.max_buffer_size / GB);
    fprintf(stderr, "  Chunk size: %.2f GB\n", (double)req.chunk_size / GB);
    fprintf(stderr, "  Region count: %d\n", req.region_count);
    
    // Cleanup
    sam_hdr_destroy(header);
    sam_close(in);
    
    return req;
}

// Update process_chromosome_region to use more efficient memory allocation
static void *process_chromosome_region(void *arg) 
{
    ThreadArg *targ = (ThreadArg *)arg;
    bam1_t *b = bam_init1();
    size_t n_records = 0;
    int ret;
    hts_itr_t *iter = NULL;  // Initialize to NULL
    samFile *fp = NULL;      // Local BAM file pointer
    
    // Allocate larger initial buffers to reduce reallocations
    size_t initial_size = targ->chunk_size;
    uint8_t *seq = malloc(initial_size);
    uint8_t *qual = malloc(initial_size);
    int *strands = malloc(initial_size * sizeof(int));
    uint32_t *positions = malloc(initial_size * sizeof(uint32_t));
    
    if (!seq || !qual || !strands || !positions) 
    {
        fprintf(stderr, "Failed to allocate sequence data buffers\n");
        goto cleanup;
    }

    // Open BAM file for this region
    fp = sam_open(targ->bam_file, "r");
    if (!fp) 
    {
        fprintf(stderr, "Failed to open BAM file: %s\n", targ->bam_file);
        goto cleanup;
    }

    // Create index for BAM file
    hts_idx_t *idx = sam_index_load(fp, targ->bam_file);
    if (!idx) 
    {
        fprintf(stderr, "Failed to load BAM index\n");
        goto cleanup;
    }

    // Create iterator for this region
    iter = sam_itr_queryi(idx, targ->tid, targ->start_pos, targ->end_pos);
    hts_idx_destroy(idx);
    if (!iter) 
    {
        fprintf(stderr, "Failed to create iterator for region %s:%u-%u\n", 
                targ->chr, targ->start_pos, targ->end_pos);
        goto cleanup;
    }

    // Process reads using iterator
    while ((ret = sam_itr_next(fp, iter, b)) >= 0) 
    {
        // Check if read is mapped and passes quality filters
        if (b->core.flag & BAM_FUNMAP || b->core.qual < targ->min_mapq)
            continue;

        // Get sequence and quality data
        uint8_t *bseq = bam_get_seq(b);
        uint8_t *bqual = bam_get_qual(b);
        int len = b->core.l_qseq;

        // Check if we need to resize buffers
        if (n_records + len > initial_size) 
        {
            // Double the size
            size_t new_size = initial_size * 2;
            uint8_t *new_seq = realloc(seq, new_size);
            uint8_t *new_qual = realloc(qual, new_size);
            int *new_strands = realloc(strands, new_size * sizeof(int));
            uint32_t *new_positions = realloc(positions, new_size * sizeof(uint32_t));
            
            if (!new_seq || !new_qual || !new_strands || !new_positions) 
            {
                // If realloc fails, process what we have and continue
                if (n_records > 0) 
                {
                    process_methylation_cpu(
                        seq, qual, strands, positions, targ->context_buffers,
                        targ->context_pos_maps, targ->min_phred, targ->min_mapq, 
                        n_records, targ->start_pos
                    );
                    n_records = 0;
                }
                continue;
            }
            
            seq = new_seq;
            qual = new_qual;
            strands = new_strands;
            positions = new_positions;
            initial_size = new_size;
        }

        // Copy sequence and quality data
        for (int i = 0; i < len; i++) 
        {
            seq[n_records + i] = seq_nt16_str[bam_seqi(bseq, i)];
            qual[n_records + i] = bqual[i];
            strands[n_records + i] = b->core.flag & BAM_FREVERSE ? 2 : 1;
            positions[n_records + i] = b->core.pos + i;
        }
        n_records += len;

        // Process in chunks to avoid memory issues
        if (n_records >= targ->chunk_size) 
        {
            process_methylation_cpu(
                seq, qual, strands, positions, targ->context_buffers,
                targ->context_pos_maps, targ->min_phred, targ->min_mapq, 
                n_records, targ->start_pos
            );
            n_records = 0;  // Reset counter after processing
        }
    }

    // Process any remaining records
    if (n_records > 0) 
    {
        process_methylation_cpu(
            seq, qual, strands, positions, targ->context_buffers,
            targ->context_pos_maps, targ->min_phred, targ->min_mapq, 
            n_records, targ->start_pos
        );
    }

cleanup:
    if (b) bam_destroy1(b);
    if (seq) free(seq);
    if (qual) free(qual);
    if (strands) free(strands);
    if (positions) free(positions);
    if (fp) sam_close(fp);
    if (iter) hts_itr_destroy(iter);
    return NULL;
}

void flush_buffer(
    const char *out_path,
    MethylRecord *buffer,
    size_t n_records,
    int compression,
    int chunk_size,
    int append_mode,
    int min_cov,
    int cap_cov,
    int min_meth,
    int max_meth,
    OutputFormat output_format
)
{
    if (n_records == 0)
        return;

    // Filter and calculate coverage
    MethylRecord *filtered_buffer = malloc(n_records * sizeof(MethylRecord));
    if (!filtered_buffer)
    {
        fprintf(stderr, "Failed to allocate filtered buffer\n");
        return;
    }

    size_t n_filtered = n_records;
    memcpy(filtered_buffer, buffer, n_records * sizeof(MethylRecord));

    // Calculate average coverage
    double avg_cov = 0.0;
    for (size_t i = 0; i < n_records; i++)
        avg_cov += filtered_buffer[i].mC + filtered_buffer[i].uC;
    avg_cov /= n_records;

    // Filter records
    filter_and_calculate_coverage_cpu(
        filtered_buffer,
        buffer,
        min_cov,
        cap_cov,
        avg_cov,
        min_meth,
        max_meth,
        &n_filtered
    );

    if (n_filtered == 0)
    {
        free(filtered_buffer);
        return;
    }

    // Output based on format
    if (output_format == OUTPUT_HDF5 || output_format == OUTPUT_BOTH)
    {
        // Create HDF5 file
        hid_t file_id = H5Fcreate(out_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (file_id < 0)
        {
            fprintf(stderr, "Failed to create HDF5 file: %s\n", out_path);
            free(filtered_buffer);
            return;
        }

        // Create dataset
        hsize_t dims[1] = {n_filtered};
        hid_t space_id = H5Screate_simple(1, dims, NULL);
        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        
        // Calculate appropriate chunk size (must be <= dataset size)
        hsize_t chunk_dims[1] = {(hsize_t)chunk_size};
        if (chunk_dims[0] > n_filtered)
            chunk_dims[0] = n_filtered;
        
        H5Pset_chunk(plist_id, 1, chunk_dims);
        H5Pset_deflate(plist_id, compression);

        // Create compound type for MethylRecord
        hid_t type_id = H5Tcreate(H5T_COMPOUND, sizeof(MethylRecord));
        H5Tinsert(type_id, "pos", HOFFSET(MethylRecord, pos), H5T_NATIVE_UINT32);
        H5Tinsert(type_id, "mC", HOFFSET(MethylRecord, mC), H5T_NATIVE_UINT16);
        H5Tinsert(type_id, "uC", HOFFSET(MethylRecord, uC), H5T_NATIVE_UINT16);
        H5Tinsert(type_id, "tnc", HOFFSET(MethylRecord, tnc), H5T_NATIVE_UINT8);

        // Create dataset and write data
        hid_t dset_id = H5Dcreate(
            file_id, 
            "methylation_data", 
            type_id, 
            space_id, 
            H5P_DEFAULT, 
            plist_id, 
            H5P_DEFAULT
        );
        if (dset_id < 0)
        {
            fprintf(stderr, "Failed to create dataset in HDF5 file\n");
            H5Tclose(type_id);
            H5Pclose(plist_id);
            H5Sclose(space_id);
            H5Fclose(file_id);
            free(filtered_buffer);
            return;
        }

        H5Dwrite(dset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

        // Clean up HDF5 resources
        H5Dclose(dset_id);
        H5Tclose(type_id);
        H5Pclose(plist_id);
        H5Sclose(space_id);
        H5Fclose(file_id);
    }

    if (output_format == OUTPUT_TXT || output_format == OUTPUT_BOTH)
    {
        // Create text file
        FILE *fp = fopen(out_path, append_mode ? "a" : "w");
        if (!fp)
        {
            fprintf(stderr, "Failed to create text file: %s\n", out_path);
            free(filtered_buffer);
            return;
        }

        // Write header
        fprintf(fp, "pos\tmC\tuC\tcontext\tstrand\ttnc\n");

        // Write data
        for (size_t i = 0; i < n_filtered; i++)
        {
            char trinucl[4];
            decode_trinucleotide(buffer[i].tnc.tnc, trinucl);
            fprintf(fp, "%u\t%u\t%u\t%s\t%c\t%s\n",
                buffer[i].pos,
                buffer[i].mC,
                buffer[i].uC,
                get_context_string(buffer[i].tnc.context),
                buffer[i].tnc.strand ? '-' : '+',
                trinucl
            );
        }

        fclose(fp);
    }

    free(filtered_buffer);
}

void process_chromosome(ThreadArg *targ)
{
    log_time("Starting processing of chromosome %s\n", targ->chr);
    
    // Use the global memory requirements
    targ->chunk_size = global_mem_req.chunk_size;
    
    // Pre-calculate sites for each context
    size_t sites_per_context[4] = {0}; // Index 0 unused, 1=CPG, 2=CHG, 3=CHH
    
    // Count sites per context
    for (uint32_t pos = 0; pos < targ->chr_len; pos++) 
    {
        int8_t strand_ctx;
        uint8_t tnc_val;
        int ctx = get_context(targ->chr_seq, targ->chr_len, pos, &strand_ctx, &tnc_val, targ->keep_chg, targ->keep_chh);
        if (ctx > 0 && ctx <= 3)
            sites_per_context[ctx]++;
    }
    
    // Pre-allocate buffers for each context
    MethylRecord *context_buffers[4] = {NULL}; // Index 0 unused, 1=CPG, 2=CHG, 3=CHH
    for (int ctx = 1; ctx <= 3; ctx++) 
    {
        if ((ctx == CONTEXT_CPG) ||
            (ctx == CONTEXT_CHG && targ->keep_chg) ||
            (ctx == CONTEXT_CHH && targ->keep_chh)) 
        {
            context_buffers[ctx] = calloc(sites_per_context[ctx], sizeof(MethylRecord));
            if (!context_buffers[ctx]) 
            {
                fprintf(stderr, "Failed to allocate buffer for %s context\n", get_context_string(ctx));
                goto cleanup;
            }
            
            // Initialize buffer with positions and context information
            size_t idx = 0;
            for (uint32_t pos = 0; pos < targ->chr_len && idx < sites_per_context[ctx]; pos++) 
            {
                int8_t strand_ctx;
                uint8_t tnc_val;
                int current_ctx = get_context(targ->chr_seq, targ->chr_len, pos, &strand_ctx, &tnc_val, targ->keep_chg, targ->keep_chh);
                if (current_ctx == ctx) 
                {
                    context_buffers[ctx][idx].pos = pos + 1;
                    context_buffers[ctx][idx].mC = 0;
                    context_buffers[ctx][idx].uC = 0;
                    context_buffers[ctx][idx].tnc.tnc = tnc_val;
                    context_buffers[ctx][idx].tnc.context = ctx;
                    context_buffers[ctx][idx].tnc.strand = (strand_ctx > 0) ? 0 : 1;
                    idx++;
                }
            }
        }
    }

    // Create position maps for each context
    khash_t(pos) *context_pos_maps[4] = {NULL}; // Index 0 unused, 1=CPG, 2=CHG, 3=CHH
    for (int ctx = 1; ctx <= 3; ctx++) 
    {
        if (context_buffers[ctx]) 
        {
            context_pos_maps[ctx] = kh_init(pos);
            for (size_t i = 0; i < sites_per_context[ctx]; i++) 
            {
                uint64_t key = ((uint64_t)targ->tid << 32) | (context_buffers[ctx][i].pos - 1);
                int ret;
                khint_t iter = kh_put(pos, context_pos_maps[ctx], key, &ret);
                kh_val(context_pos_maps[ctx], iter) = i;
            }
        }
    }

    uint32_t chunk_size = targ->chunk_size;
    if (chunk_size > targ->chr_len)
        chunk_size = targ->chr_len;

    int n_regions = 1;  // Initialize to 1 as minimum value
    n_regions = (int)ceil((double)targ->chr_len / chunk_size);
    if (n_regions > global_mem_req.region_count)
    {
        n_regions = global_mem_req.region_count;
        chunk_size = (targ->chr_len + n_regions - 1) / n_regions;
    }

    // Allocate thread arguments array
    ThreadArg *region_args = NULL;  // Initialize to NULL
    region_args = malloc(n_regions * sizeof(ThreadArg));
    if (!region_args) 
    {
        fprintf(stderr, "Failed to allocate region arguments\n");
        goto cleanup;
    }

    pthread_mutex_t buffer_mutex;
    pthread_mutex_init(&buffer_mutex, NULL);

    // First, validate the chromosome name
    if (!targ->chr || strlen(targ->chr) == 0) 
    {
        fprintf(stderr, "Invalid chromosome name\n");
        free(region_args);
        pthread_mutex_destroy(&buffer_mutex);
        goto cleanup;
    }

    for (int i = 0; i < n_regions; i++)
    {
        // Initialize each thread argument structure
        memset(&region_args[i], 0, sizeof(ThreadArg));
        
        // Copy all non-pointer fields
        region_args[i].bam_file = targ->bam_file;
        region_args[i].out_dir = targ->out_dir;
        region_args[i].tid = targ->tid;
        
        // Make a deep copy of the chromosome name and validate it
        region_args[i].chr = strdup(targ->chr);
        if (!region_args[i].chr || strlen(region_args[i].chr) == 0) 
        {
            fprintf(stderr, "Failed to allocate or validate chromosome name copy for region %d\n", i);
            // Clean up previously allocated regions
            for (int j = 0; j < i; j++) 
            {
                free(region_args[j].chr);
            }
            free(region_args);
            pthread_mutex_destroy(&buffer_mutex);
            goto cleanup;
        }

        region_args[i].chr_len = targ->chr_len;
        region_args[i].chr_seq = targ->chr_seq;
        region_args[i].min_mapq = targ->min_mapq;
        region_args[i].min_phred = targ->min_phred;
        region_args[i].min_cov = targ->min_cov;
        region_args[i].cap_cov = targ->cap_cov;
        region_args[i].min_meth = targ->min_meth;
        region_args[i].max_meth = targ->max_meth;
        region_args[i].keep_chg = targ->keep_chg;
        region_args[i].keep_chh = targ->keep_chh;
        region_args[i].hdf5_compression = targ->hdf5_compression;
        region_args[i].hdf5_chunk_size = targ->hdf5_chunk_size;
        region_args[i].chunk_size = targ->chunk_size;
        region_args[i].output_format = targ->output_format;
        region_args[i].split_context_files = targ->split_context_files;

        // Set shared resources
        region_args[i].buffer_mutex = &buffer_mutex;
        region_args[i].context_buffers = context_buffers;
        region_args[i].context_pos_maps = context_pos_maps;

        uint32_t start_pos = i * chunk_size;
        uint32_t end_pos = (i == n_regions - 1) ? targ->chr_len : ((i + 1) * chunk_size);

        if (end_pos > targ->chr_len)
            end_pos = targ->chr_len;

        region_args[i].start_pos = start_pos;
        region_args[i].end_pos = end_pos;
    }

    pthread_t *threads = NULL;  // Initialize to NULL
    threads = malloc(n_regions * sizeof(pthread_t));
    if (!threads)
    {
        fprintf(stderr, "Failed to allocate threads\n");
        for (int i = 0; i < n_regions; i++) 
            free(region_args[i].chr);
        free(region_args);
        pthread_mutex_destroy(&buffer_mutex);
        goto cleanup;
    }

    int active_threads = 0;
    for (int i = 0; i < n_regions; i++)
    {
        if (pthread_create(&threads[i], NULL, (void *(*)(void *))process_chromosome_region, &region_args[i]) != 0)
        {
            fprintf(stderr, "Failed to create thread for %s\n", region_args[i].chr);
            continue;
        }
        active_threads++;
        for (int j = 0; j <= i; j++)
            if (pthread_join(threads[j], NULL) == 0)
            {
                active_threads--;
                break;
            }
        
        int max_region_threads = 8;
        while (active_threads >= max_region_threads)
            for (int j = 0; j <= i; j++)
                if (pthread_join(threads[j], NULL) == 0)
                {
                    active_threads--;
                    break;
                }
    }

    pthread_mutex_destroy(&buffer_mutex);

    // Process each context
    for (int ctx = 1; ctx <= 3; ctx++) 
    {
        if (context_buffers[ctx] && 
            ((ctx == CONTEXT_CPG) ||
             (ctx == CONTEXT_CHG && targ->keep_chg) ||
             (ctx == CONTEXT_CHH && targ->keep_chh))) 
        {
            log_time("Processing %s context for chromosome %s\n", get_context_string(ctx), targ->chr);
            
            // Output file name - use appropriate extension based on output format
            char out_path[1024];
            const char *ext = (targ->output_format == OUTPUT_TXT) ? ".txt" : ".h5";  // For OUTPUT_BOTH, use .h5
            snprintf(
                out_path, 
                sizeof(out_path), 
                "%s/%s-%s%s", 
                targ->out_dir, 
                targ->chr, 
                get_context_string(ctx),
                ext
            );
            
            flush_buffer(
                out_path,
                context_buffers[ctx],
                sites_per_context[ctx],
                targ->hdf5_compression,
                targ->hdf5_chunk_size,
                0,
                targ->min_cov,
                targ->cap_cov,
                targ->min_meth,
                targ->max_meth,
                targ->output_format
            );
            
            log_time("Finished processing %s context for chromosome %s\n", get_context_string(ctx), targ->chr);
        }
    }

cleanup:
    // Clean up region arguments
    if (region_args) {
        for (int i = 0; i < n_regions; i++)
            if (region_args[i].chr) 
                free(region_args[i].chr);
        free(region_args);
    }
    if (threads) free(threads);
    
    // Clean up context buffers and position maps
    for (int ctx = 1; ctx <= 3; ctx++) 
    {
        if (context_buffers[ctx])
            free(context_buffers[ctx]);
        if (context_pos_maps[ctx])
            kh_destroy(pos, context_pos_maps[ctx]);
    }
    
    log_time("Finished processing chromosome %s\n", targ->chr);
}

void cleanup_hdf5(void)
{
    H5close();
}

int load_chrom_mapping(const char *filename, ChromMapEntry **entries, int *n_entries, char **reference_file) 
{
    FILE *fp = fopen(filename, "r");
    if (!fp) 
        return -1;
    fseek(fp, 0, SEEK_END);
    long len = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    char *data = malloc(len + 1);
    size_t nread = fread(data, 1, len, fp);
    if (nread != len) 
    {
        free(data);
        fclose(fp);
        return -2; // Error: could not read the expected number of bytes
    }
    data[len] = 0;
    fclose(fp);

    cJSON *json = cJSON_Parse(data);
    free(data);
    if (!json) 
        return -3;

    // Get reference file path
    cJSON *ref = cJSON_GetObjectItem(json, "reference");
    if (ref && cJSON_IsString(ref) && ref->valuestring) 
        *reference_file = strdup(ref->valuestring);

    cJSON *chroms = cJSON_GetObjectItem(json, "chromosomes");
    if (!chroms || !cJSON_IsArray(chroms)) 
    {
        cJSON_Delete(json);
        return -4;
    }
    int count = cJSON_GetArraySize(chroms);
    *entries = calloc(count, sizeof(ChromMapEntry));
    *n_entries = 0;
    for (int i = 0; i < count; ++i) 
    {
        cJSON *item = cJSON_GetArrayItem(chroms, i);
        if (!cJSON_IsObject(item)) 
            continue;
        
        ChromMapEntry *e = &(*entries)[*n_entries];
        
        // Initialize all strings to empty
        e->fasta[0] = '\0';
        e->bam[0] = '\0';
        e->name[0] = '\0';
        e->extract = 0;  // Default to false
        
        cJSON *fasta = cJSON_GetObjectItem(item, "fasta");
        cJSON *bam = cJSON_GetObjectItem(item, "bam");
        cJSON *name = cJSON_GetObjectItem(item, "name");
        cJSON *extract = cJSON_GetObjectItem(item, "extract");
        
        // Validate and copy each field
        if (fasta && cJSON_IsString(fasta) && fasta->valuestring) 
        {
            strncpy(e->fasta, fasta->valuestring, sizeof(e->fasta) - 1);
            e->fasta[sizeof(e->fasta) - 1] = '\0';
        }
        
        if (bam && cJSON_IsString(bam) && bam->valuestring) 
        {
            strncpy(e->bam, bam->valuestring, sizeof(e->bam) - 1);
            e->bam[sizeof(e->bam) - 1] = '\0';
        }
        
        if (name && cJSON_IsString(name) && name->valuestring) 
        {
            strncpy(e->name, name->valuestring, sizeof(e->name) - 1);
            e->name[sizeof(e->name) - 1] = '\0';
        }
        
        // Handle optional extract field
        if (extract && cJSON_IsBool(extract))
            e->extract = cJSON_IsTrue(extract);
        
        // Validate that we have all required fields
        if (e->fasta[0] == '\0' || e->bam[0] == '\0' || e->name[0] == '\0') 
        {
            fprintf(stderr, "Warning: Skipping chromosome entry with missing required fields\n");
            continue;
        }
        
        (*n_entries)++;
    }
    cJSON_Delete(json);
    return 0;
}

static void print_usage(const char *prog) 
{
    fprintf(stderr, "Usage: %s [options] <input.bam> <output_dir>\n", prog);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h, --help                Show this help message\n");
    fprintf(stderr, "  -q, --min-mapq INT        Minimum mapping quality [%d]\n", DEFAULT_MIN_MAPQ);
    fprintf(stderr, "  -p, --min-phred INT       Minimum base quality [%d]\n", DEFAULT_MIN_PHRED);
    fprintf(stderr, "  -c, --min-cov INT         Minimum coverage [%d]\n", DEFAULT_MIN_COV);
    fprintf(stderr, "  -C, --cap-cov INT         Cap coverage at this value [0 = no cap]\n");
    fprintf(stderr, "  -m, --min-meth FLOAT      Minimum methylation level [%d]\n", DEFAULT_MIN_METH);
    fprintf(stderr, "  -M, --max-meth FLOAT      Maximum methylation level [%d]\n", DEFAULT_MAX_METH);
    fprintf(stderr, "  -z, --compression INT     HDF5 compression level [%d]\n", DEFAULT_HDF5_COMPRESSION);
    fprintf(stderr, "  -k, --chunk-size INT      HDF5 chunk size [%zu]\n", DEFAULT_HDF5_CHUNK_SIZE);
    fprintf(stderr, "  -a, --append              Append to existing files\n");
    fprintf(stderr, "  -f, --output-format STR   Output format (hdf5, txt, both) [hdf5]\n");
    fprintf(stderr, "  -S, --split-context-files Split output by context\n");
    fprintf(stderr, "  -x, --chrom-mapping FILE  Chromosome mapping file\n");
    fprintf(stderr, "  -O, --output-dir DIR      Output directory\n");
    fprintf(stderr, "  -G, --CHG                 Process CHG context\n");
    fprintf(stderr, "  -H, --CHH                 Process CHH context\n");
#ifdef USE_CUDA
    fprintf(stderr, "  -g, --gpu                 Use GPU acceleration\n");
#endif
    fprintf(stderr, "\n");
}

int main(int argc, char *argv[])
{
#ifdef USE_CUDA
    // Initialize CUDA device
    if (init_cuda_device() != 0) 
        fprintf(stderr, "Failed to initialize CUDA device. Running in CPU-only mode.\n");
#endif

    fprintf(stderr, "Program: MethylExtractor\nParameters:\n");
    for (int i = 0; i < argc; i++)
        fprintf(stderr, "  Arg %d: %s\n", i, argv[i]);

    log_time("Starting processing...\n");

    int hdf5_compression = DEFAULT_HDF5_COMPRESSION;
    int hdf5_chunk_size = DEFAULT_HDF5_CHUNK_SIZE;
    uint32_t chunk_size = DEFAULT_CHUNK_SIZE;
    int num_threads = DEFAULT_THREADS;
    int keep_chg = 0;
    int keep_chh = 0;
    int min_mapq = DEFAULT_MIN_MAPQ;
    int min_phred = DEFAULT_MIN_PHRED;
    int min_cov = DEFAULT_MIN_COV;
    int cap_cov = DEFAULT_CAP_COVERAGE;
    int min_meth = DEFAULT_MIN_METH;
    int max_meth = DEFAULT_MAX_METH;
    OutputFormat output_format = OUTPUT_HDF5;  // Default to HDF5 output
    const char *out_dir = NULL;
    int split_context_files = 0;
    const char *chrom_mapping_file = NULL;
    int use_gpu = 0;
    char *reference_file = NULL;

    static struct option long_options[] = 
    {
        {"help", no_argument, 0, 'h'},
        {"min-mapq", required_argument, 0, 'q'},
        {"min-phred", required_argument, 0, 'p'},
        {"min-cov", required_argument, 0, 'c'},
        {"cap-cov", required_argument, 0, 'C'},
        {"min-meth", required_argument, 0, 'm'},
        {"max-meth", required_argument, 0, 'M'},
        {"compression", required_argument, 0, 'z'},
        {"chunk-size", required_argument, 0, 'k'},
        {"append", no_argument, 0, 'a'},
        {"output-format", required_argument, 0, 'f'},
        {"split-context-files", no_argument, 0, 'S'},
        {"chrom-mapping", required_argument, 0, 'x'},
        {"CHG", no_argument, 0, 'G'},
        {"CHH", no_argument, 0, 'H'},
        {"output-dir", required_argument, 0, 'O'},
#ifdef USE_CUDA
        {"gpu", no_argument, 0, 'g'},
#endif
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "hq:p:c:C:m:M:z:k:af:x:GHSO:g", long_options, NULL)) != -1) 
    {
        switch (opt)
        {
        case 'G':
            keep_chg = 1;
            break;
        case 'H':
            keep_chh = 1;
            break;
        case 'S':
            split_context_files = 1;
            break;
        case 'O':
            out_dir = optarg;
            break;
        case 'f':
            if (strcmp(optarg, "hdf5") == 0)
                output_format = OUTPUT_HDF5;
            else if (strcmp(optarg, "txt") == 0)
                output_format = OUTPUT_TXT;
            else if (strcmp(optarg, "both") == 0)
                output_format = OUTPUT_BOTH;
            else
            {
                fprintf(stderr, "Invalid output format '%s'. Must be one of: hdf5, txt, both\n", optarg);
                return 1;
            }
            break;
        case 'z':
            hdf5_compression = atoi(optarg);
            break;
        case 'k':
            hdf5_chunk_size = atoi(optarg);
            if (hdf5_chunk_size < 1)
            {
                fprintf(stderr, "HDF5 chunk size must be positive\n");
                return 1;
            }
            break;
        case 's':
            chunk_size = atoi(optarg);
            if (chunk_size < 1)
            {
                fprintf(stderr, "Chunk size must be positive\n");
                return 1;
            }
            break;
        case 't':
            num_threads = atoi(optarg);
            if (num_threads < 1)
            {
                fprintf(stderr, "Number of threads must be positive\n");
                return 1;
            }
            break;
        case 'q':
            min_mapq = atoi(optarg);
            if (min_mapq < 0)
            {
                fprintf(stderr, "Minimum mapping quality must be non-negative\n");
                return 1;
            }
            break;
        case 'p':
            min_phred = atoi(optarg);
            if (min_phred < 0)
            {
                fprintf(stderr, "Minimum Phred score must be non-negative\n");
                return 1;
            }
            break;
        case 'c':
            min_cov = atoi(optarg);
            if (min_cov < 0)
            {
                fprintf(stderr, "Minimum coverage must be positive\n");
                return 1;
            }
            break;
        case 'C':
            cap_cov = atoi(optarg);
            break;
        case 'm':
            min_meth = atoi(optarg);
            if (min_meth < 0 || min_meth > 100)
            {
                fprintf(stderr, "Minimum methylation level must be between 0 and 100\n");
                return 1;
            }
            break;
        case 'M':
            max_meth = atoi(optarg);
            if (max_meth < 0 || max_meth > 100)
            {
                fprintf(stderr, "Maximum methylation level must be between 0 and 100\n");
                return 1;
            }
            break;
        case 'x':
            chrom_mapping_file = optarg;
            break;
        case 'h':
            print_usage(argv[0]);
            return 0;
        case '?':
        default:
            print_usage(argv[0]);
            return 1;
        }
    }
    if (argc - optind < 2)
    {
        print_usage(argv[0]);
        return 1;
    }

    if (min_meth > max_meth)
    {
        fprintf(stderr, "Minimum methylation level (%d) must not exceed maximum methylation level (%d)\n", min_meth, max_meth);
        return 1;
    }

    // If there are two arguments, the first is the reference file and the second is the BAM file
    const char *bam_file = NULL;
    const char *ref_file = NULL;
    if (argc - optind == 2)
    {
        ref_file = argv[optind];
        bam_file = argv[optind + 1];
    }
    else
        bam_file = argv[optind];

    if (num_threads == DEFAULT_THREADS)
    {
        num_threads = sysconf(_SC_NPROCESSORS_ONLN);
        if (num_threads < 1)
            num_threads = DEFAULT_THREADS;
    }
    if (make_directory(out_dir) != 0)
    {
        fprintf(stderr, "Failed to create output directory: %s\n", out_dir);
        return 1;
    }

    // Load chromosome mapping
    ChromMapEntry *chroms = NULL;
    int n_chroms = 0;
    if (!chrom_mapping_file || load_chrom_mapping(chrom_mapping_file, &chroms, &n_chroms, &reference_file) != 0) 
    {
        fprintf(stderr, "Failed to load chromosome mapping from %s\n", chrom_mapping_file);
        return 1;
    }

    // Use reference file from JSON if provided
    if (!ref_file && reference_file)
    {
        ref_file = reference_file;
        fprintf(stderr, "Using reference file from chromosome mapping: %s\n", ref_file);
    }

    if (!ref_file)
    {
        fprintf(stderr, "Reference file not provided\n");
        return 1;
    }

    // Load reference file
    faidx_t *fai = fai_load(ref_file);
    if (!fai)
    {
        fprintf(stderr, "Failed to load reference: %s\n", ref_file);
        return 1;
    }
    samFile *in = sam_open(bam_file, "r");
    if (!in)
    {
        fprintf(stderr, "Failed to open BAM: %s\n", bam_file);
        fai_destroy(fai);
        return 1;
    }
    bam_hdr_t *header = sam_hdr_read(in);
    if (!header)
    {
        fprintf(stderr, "Failed to read BAM header\n");
        sam_close(in);
        fai_destroy(fai);
        return 1;
    }
    sam_close(in);

    // Calculate memory requirements once using the first chromosome
    if (n_chroms > 0) 
    {
        ChromMapEntry *first_chrom = &chroms[0];
        int tid = bam_name2id(header, first_chrom->bam);
        if (tid >= 0) 
        {
            global_mem_req = calculate_memory_requirements(
                bam_file,
                first_chrom->bam,
                header->target_len[tid],
                min_mapq,
                min_phred
            );
            fprintf(stderr, "Global memory settings (will be used for all chromosomes):\n");
            fprintf(stderr, "  Min buffer size: %.2f GB\n", (double)global_mem_req.min_buffer_size / GB);
            fprintf(stderr, "  Optimal buffer size: %.2f GB\n", (double)global_mem_req.optimal_buffer_size / GB);
            fprintf(stderr, "  Max buffer size: %.2f GB\n", (double)global_mem_req.max_buffer_size / GB);
            fprintf(stderr, "  Chunk size: %.2f GB\n", (double)global_mem_req.chunk_size / GB);
            fprintf(stderr, "  Region count: %d\n", global_mem_req.region_count);
        }
    }

    ThreadArg *thread_args = malloc(header->n_targets * sizeof(ThreadArg));
    if (!thread_args)
    {
        fprintf(stderr, "Failed to allocate thread arguments\n");
        sam_hdr_destroy(header);
        fai_destroy(fai);
        return 1;
    }
    int valid_chr_count = 0;
    for (int i = 0; i < n_chroms; ++i) 
    {
        ChromMapEntry *entry = &chroms[i];
        // Find BAM tid for entry->bam
        int tid = bam_name2id(header, entry->bam);
        if (tid < 0) 
        {
            fprintf(stderr, "BAM does not contain chromosome %s (looking for %s in BAM)\n", entry->bam, entry->bam);
            fprintf(stderr, "Available chromosomes in BAM:\n");
            for (int j = 0; j < header->n_targets; j++) {
                fprintf(stderr, "  %s\n", header->target_name[j]);
            }
            continue;
        }
        // Fetch sequence for entry->fasta
        int seq_len = 0;
        char *seq = faidx_fetch_seq(fai, entry->fasta, 0, header->target_len[tid], &seq_len);
        if (!seq || seq_len <= 0) 
        {
            fprintf(stderr, "Failed to fetch sequence for %s\n", entry->fasta);
            if (seq) free(seq);
            continue;
        }
        // Set up ThreadArg as before, but use entry->name for output
        thread_args[valid_chr_count].bam_file = bam_file;
        thread_args[valid_chr_count].out_dir = out_dir;
        thread_args[valid_chr_count].tid = tid;
        thread_args[valid_chr_count].chr = strdup(entry->name);  // Make a copy of the name
        thread_args[valid_chr_count].chr_len = header->target_len[tid];
        thread_args[valid_chr_count].min_mapq = min_mapq;
        thread_args[valid_chr_count].min_phred = min_phred;
        thread_args[valid_chr_count].min_cov = min_cov;
        thread_args[valid_chr_count].cap_cov = cap_cov;
        thread_args[valid_chr_count].min_meth = min_meth;
        thread_args[valid_chr_count].max_meth = max_meth;
        thread_args[valid_chr_count].keep_chg = keep_chg;
        thread_args[valid_chr_count].keep_chh = keep_chh;
        thread_args[valid_chr_count].hdf5_compression = hdf5_compression;
        thread_args[valid_chr_count].hdf5_chunk_size = hdf5_chunk_size;
        thread_args[valid_chr_count].chunk_size = chunk_size;
        thread_args[valid_chr_count].chr_seq = seq;
        thread_args[valid_chr_count].output_format = output_format;
        thread_args[valid_chr_count].split_context_files = split_context_files;
        thread_args[valid_chr_count].use_gpu = use_gpu;
        valid_chr_count++;
    }
    free(chroms);
    pthread_t *threads = malloc(valid_chr_count * sizeof(pthread_t));
    if (!threads)
    {
        fprintf(stderr, "Failed to allocate threads\n");
        for (int i = 0; i < valid_chr_count; i++)
            free(thread_args[i].chr_seq);
        free(thread_args);
        sam_hdr_destroy(header);
        fai_destroy(fai);
        return 1;
    }

    int active_threads = 0;
    for (int i = 0; i < valid_chr_count; i++)
    {
        if (pthread_create(&threads[i], NULL, (void *(*)(void *))process_chromosome, &thread_args[i]) != 0)
        {
            fprintf(stderr, "Failed to create thread for %s\n", thread_args[i].chr);
            continue;
        }
        active_threads++;
        for (int j = 0; j <= i; j++)
            if (pthread_join(threads[j], NULL) == 0)
            {
                active_threads--;
                break;
            }
        
        int max_region_threads = 8;
        while (active_threads >= max_region_threads)
            for (int j = 0; j <= i; j++)
                if (pthread_join(threads[j], NULL) == 0)
                {
                    active_threads--;
                    break;
                }
    }

    cleanup_hdf5();
    
    for (int i = 0; i < valid_chr_count; i++)
    {
        free(thread_args[i].chr_seq);
        free(thread_args[i].chr);  // Free the copied chromosome name
    }
    free(threads);
    free(thread_args);
    sam_hdr_destroy(header);
    fai_destroy(fai);

    // Clean up reference file string if it was allocated
    if (reference_file)
        free(reference_file);

    log_time("Processing complete. MethylExtractor has finished.\n");
        
    return 0;
}
