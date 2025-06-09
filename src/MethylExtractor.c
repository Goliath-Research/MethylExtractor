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
#include <stdint.h>  // For int64_t

// Memory size constants
#define KB (1024ULL)
#define MB (KB * 1024ULL)
#define GB (MB * 1024ULL)
#define TB (GB * 1024ULL)

// Memory allocation constants
#define MIN_BUFFER_SIZE (1ULL * MB)      // Initial buffer size
#define MAX_BUFFER_SIZE (1ULL * GB)      // Maximum buffer size
#define BUFFER_GROWTH_FACTOR 2           // Buffer growth multiplier

// Memory allocation constants
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
    return cpu_count > 0 ? cpu_count : 24; // Default to 64 if detection fails
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
KHASH_INIT(pos, uint64_t, size_t, 1, kh_int64_hash_func, kh_int64_hash_equal)  // Use uint64_t keys
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
    size_t chunk_size;
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
} MemoryRequirements;

// Add global memory requirements
static MemoryRequirements global_mem_req = {0};

// Structure to hold chromosome processing results
typedef struct {
    char *chr_name;
    int tid;
    MethylRecord **context_buffers;  // Array of 3 buffers (CPG, CHG, CHH)
    size_t *context_sizes;          // Size of each context buffer
    khash_t(pos) **context_pos_maps; // Array of position maps for each context
    int processed;                  // Flag to indicate processing is complete
    pthread_mutex_t mutex;          // Mutex for thread safety
} ChromosomeResult;

// Structure for the processing queue
typedef struct {
    ChromosomeResult *results;      // Array of chromosome results
    int total_chromosomes;          // Total number of chromosomes
    int processed_count;            // Count of processed chromosomes
    pthread_mutex_t queue_mutex;    // Mutex for queue operations
    pthread_cond_t queue_cond;      // Condition variable for queue synchronization
} ProcessingQueue;

// Structure for thread arguments
typedef struct 
{
    ProcessingQueue *queue;         // Pointer to the processing queue
    const char *bam_file;          // BAM file path
    const char *out_dir;           // Output directory
    const char *ref_file;          // Reference file path
    const char *chr_name;          // BAM chromosome name
    const char *fasta_chr_name;    // FASTA chromosome name
    const char *display_name;      // Display name for output files
    int min_mapq;                  // Minimum mapping quality
    int min_phred;                 // Minimum Phred score
    int keep_chg;                  // Whether to process CHG context
    int keep_chh;                  // Whether to process CHH context
    int thread_id;                 // Thread identifier
    pthread_mutex_t fai_mutex;     // Per-thread mutex for FASTA operations
    int use_gpu;                   // GPU flag
    int split_context_files;       // Split context files flag
    OutputFormat output_format;    // Output format
    int tid;                       // Chromosome ID
    uint32_t chunk_size;          // Processing chunk size
    samFile *fp;                   // BAM file pointer
    bam_hdr_t *header;            // BAM header
    hts_idx_t *idx;               // BAM index
    faidx_t *fai;                 // FASTA index
    char *chr_seq;                // Chromosome sequence
    uint32_t chr_len;             // Chromosome length
    ChromosomeResult *result;     // Result structure
    int min_cov;                  // Minimum coverage
    int cap_cov;                  // Coverage cap
    int min_meth;                 // Minimum methylation
    int max_meth;                 // Maximum methylation
    int hdf5_compression;         // HDF5 compression level
    int hdf5_chunk_size;          // HDF5 chunk size
} ThreadContext;

// Function declarations
static void *process_chromosome_thread(void *arg);

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
    int keep_chh,
    khash_t(pos) **context_pos_maps
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

            // Add position to the appropriate context map
            if (context_pos_maps && context_pos_maps[ctx]) 
            {
                // Construct key as uint64_t: (pos << 32) | (pos - 1)
                uint64_t key = ((uint64_t)pos << 32) | ((uint64_t)(pos - 1));
                int ret;
                khint_t k = kh_put(pos, context_pos_maps[ctx], key, &ret);
                if (ret == -1) 
                {
                    // Hash table is full, resize it
                    kh_resize(pos, context_pos_maps[ctx], kh_size(context_pos_maps[ctx]) * 2);
                    k = kh_put(pos, context_pos_maps[ctx], key, &ret);
                }
                if (ret >= 0) 
                {
                    kh_val(context_pos_maps[ctx], k) = idx;
                }
            }

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
                
                uint64_t key = ((uint64_t)positions[i] << 32) | (positions[i] - 1);
                khint_t iter = kh_get(pos, context_pos_maps[ctx], key);
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
                 (base == 'A' && (strand == 2 || strands[i] == 4)))
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

// Update process_chromosome_region signature
static void process_chromosome_region(
    ThreadContext *ctx,
    samFile *in,
    bam_hdr_t *header,
    hts_idx_t *idx,
    faidx_t *fai,
    char *chr_seq,
    uint32_t chr_len,
    ChromosomeResult *result
);

static void process_chromosome_region(
    ThreadContext *ctx,
    samFile *in,
    bam_hdr_t *header,
    hts_idx_t *idx,
    faidx_t *fai,
    char *chr_seq,
    uint32_t chr_len,
    ChromosomeResult *result
) 
{
    bam1_t *b = bam_init1();
    size_t n_records = 0;
    int ret;
    hts_itr_t *iter = NULL;

    // Allocate buffers
    size_t initial_size = MIN_BUFFER_SIZE;
    uint8_t *seq = malloc(initial_size);
    uint8_t *qual = malloc(initial_size);
    int *strands = malloc(initial_size * sizeof(int));
    uint32_t *positions = malloc(initial_size * sizeof(uint32_t));

    if (!seq || !qual || !strands || !positions) 
    {
        fprintf(stderr, "Failed to allocate sequence data buffers\n");
        goto cleanup;
    }

    // Create iterator for this region
    iter = sam_itr_queryi(idx, ctx->tid, 0, chr_len);
    if (!iter) 
    {
        fprintf(stderr, "Failed to create iterator for chromosome %s\n", ctx->chr_name);
        goto cleanup;
    }

    // Process reads
    while ((ret = sam_itr_next(in, iter, b)) >= 0) 
    {
        // Process read
        if (b->core.flag & BAM_FUNMAP || b->core.qual < ctx->min_mapq)
            continue;

        // Get sequence and quality data
        uint8_t *bseq = bam_get_seq(b);
        uint8_t *bqual = bam_get_qual(b);
        int len = b->core.l_qseq;

        // Check if we need to resize buffers
        if (n_records + len > initial_size) 
        {
            size_t new_size = initial_size * 2;
            uint8_t *new_seq = realloc(seq, new_size);
            uint8_t *new_qual = realloc(qual, new_size);
            int *new_strands = realloc(strands, new_size * sizeof(int));
            uint32_t *new_positions = realloc(positions, new_size * sizeof(uint32_t));

            if (!new_seq || !new_qual || !new_strands || !new_positions) 
            {
                fprintf(stderr, "Failed to resize buffers\n");
                goto cleanup;
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

        // Process in chunks
        if (n_records >= ctx->chunk_size) 
        {
            process_methylation_cpu(
                seq, qual, strands, positions, result->context_buffers,
                result->context_pos_maps, ctx->min_phred, ctx->min_mapq, 
                n_records, 0
            );
            n_records = 0;
        }
    }

    // Process remaining records
    if (n_records > 0) 
    {
        process_methylation_cpu(
            seq, qual, strands, positions, result->context_buffers,
            result->context_pos_maps, ctx->min_phred, ctx->min_mapq, 
            n_records, 0
        );
    }

cleanup:
    if (b) bam_destroy1(b);
    if (seq) free(seq);
    if (qual) free(qual);
    if (strands) free(strands);
    if (positions) free(positions);
    if (iter) hts_itr_destroy(iter);
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

static void process_chromosome(ThreadContext *ctx)
{
    log_time("Starting processing of chromosome %s\n", ctx->chr_name);
    
    // Use the global memory requirements
    ctx->chunk_size = global_mem_req.chunk_size;
    
    // Pre-calculate sites for each context
    size_t sites_per_context[4] = {0}; // Index 0 unused, 1=CPG, 2=CHG, 3=CHH
    
    // Count sites per context
    for (uint32_t pos = 0; pos < ctx->chr_len; pos++) 
    {
        int8_t strand_ctx;
        uint8_t tnc_val;
        int ctx_type = get_context(ctx->chr_seq, ctx->chr_len, pos, &strand_ctx, &tnc_val, ctx->keep_chg, ctx->keep_chh);
        if (ctx_type > 0 && ctx_type <= 3)
            sites_per_context[ctx_type]++;
    }
    
    // Pre-allocate buffers for each context
    for (int ctx_type = 1; ctx_type <= 3; ctx_type++) 
    {
        if ((ctx_type == CONTEXT_CPG) ||
            (ctx_type == CONTEXT_CHG && ctx->keep_chg) ||
            (ctx_type == CONTEXT_CHH && ctx->keep_chh)) 
        {
            // Initialize position map first
            ctx->result->context_pos_maps[ctx_type] = kh_init(pos);
            if (!ctx->result->context_pos_maps[ctx_type]) 
            {
                fprintf(stderr, "Failed to initialize position map for %s context\n", get_context_string(ctx_type));
                goto cleanup;
            }

            // Allocate buffer
            ctx->result->context_buffers[ctx_type] = calloc(sites_per_context[ctx_type], sizeof(MethylRecord));
            if (!ctx->result->context_buffers[ctx_type]) 
            {
                fprintf(stderr, "Failed to allocate buffer for %s context\n", get_context_string(ctx_type));
                kh_destroy(pos, ctx->result->context_pos_maps[ctx_type]);
                ctx->result->context_pos_maps[ctx_type] = NULL;
                goto cleanup;
            }
            
            // Initialize buffer with positions and context information
            initialize_buffer(
                ctx->result->context_buffers[ctx_type],
                sites_per_context[ctx_type],
                ctx->chr_seq,
                ctx->chr_len,
                ctx->keep_chg,
                ctx->keep_chh,
                ctx->result->context_pos_maps
            );
            
            ctx->result->context_sizes[ctx_type] = sites_per_context[ctx_type];
        }
    }

    // Process the chromosome data
    process_chromosome_region(ctx, ctx->fp, ctx->header, ctx->idx, ctx->fai, ctx->chr_seq, ctx->chr_len, ctx->result);

    // Process each context
    for (int ctx_type = 1; ctx_type <= 3; ctx_type++) 
    {
        if (ctx->result->context_buffers[ctx_type] && 
            ((ctx_type == CONTEXT_CPG) ||
             (ctx_type == CONTEXT_CHG && ctx->keep_chg) ||
             (ctx_type == CONTEXT_CHH && ctx->keep_chh)))
        {
            log_time("Processing %s context for chromosome %s\n", get_context_string(ctx_type), ctx->display_name);
            
            // Output file name - use appropriate extension based on output format
            char out_path[1024];
            const char *ext = (ctx->output_format == OUTPUT_TXT) ? ".txt" : ".h5";  // For OUTPUT_BOTH, use .h5
            snprintf(
                out_path, 
                sizeof(out_path), 
                "%s/%s-%s%s", 
                ctx->out_dir, 
                ctx->display_name,  // Use display name for output files
                get_context_string(ctx_type),
                ext
            );
            
            flush_buffer(
                out_path,
                ctx->result->context_buffers[ctx_type],
                ctx->result->context_sizes[ctx_type],
                ctx->hdf5_compression,
                ctx->hdf5_chunk_size,
                0,
                ctx->min_cov,
                ctx->cap_cov,
                ctx->min_meth,
                ctx->max_meth,
                ctx->output_format
            );
            
            log_time("Finished processing %s context for chromosome %s\n", get_context_string(ctx_type), ctx->display_name);
        }
    }

cleanup:
    // Clean up context buffers and position maps
    for (int ctx_type = 1; ctx_type <= 3; ctx_type++) 
    {
        if (ctx->result->context_buffers[ctx_type])
            free(ctx->result->context_buffers[ctx_type]);
        if (ctx->result->context_pos_maps[ctx_type])
            kh_destroy(pos, ctx->result->context_pos_maps[ctx_type]);
    }
    
    log_time("Finished processing chromosome %s\n", ctx->chr_name);
}

void cleanup_hdf5(void)
{
    H5close();
}

int load_chrom_mapping(const char *filename, ChromMapEntry **entries, int *n_entries, char **reference_file) 
{
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Failed to open chromosome mapping file: %s\n", filename);
        return -1;
    }

    // Read the entire file
    fseek(fp, 0, SEEK_END);
    long len = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    
    char *data = malloc(len + 1);
    if (!data) {
        fprintf(stderr, "Failed to allocate memory for file data\n");
        fclose(fp);
        return -2;
    }
    
    size_t nread = fread(data, 1, len, fp);
    if (nread != len) {
        fprintf(stderr, "Failed to read chromosome mapping file\n");
        free(data);
        fclose(fp);
        return -2;
    }
    data[len] = 0;
    fclose(fp);

    // Parse JSON
    cJSON *json = cJSON_Parse(data);
    free(data);
    if (!json) 
    {
        fprintf(stderr, "Failed to parse chromosome mapping JSON\n");
        return -3;
    }

    // Get reference file path
    cJSON *ref = cJSON_GetObjectItem(json, "reference");
    if (ref && cJSON_IsString(ref) && ref->valuestring) 
    {
        *reference_file = strdup(ref->valuestring);
        fprintf(stderr, "Found reference file: %s\n", *reference_file);
    }
    else
        fprintf(stderr, "Warning: No reference file specified in mapping\n");

    // Get chromosomes array
    cJSON *chroms = cJSON_GetObjectItem(json, "chromosomes");
    if (!chroms || !cJSON_IsArray(chroms)) 
    {
        fprintf(stderr, "No chromosomes array found in mapping\n");
        cJSON_Delete(json);
        return -4;
    }

    int count = cJSON_GetArraySize(chroms);
    fprintf(stderr, "Found %d chromosome mappings\n", count);
    
    // First pass: count valid entries
    int valid_count = 0;
    for (int i = 0; i < count; ++i) {
        cJSON *item = cJSON_GetArrayItem(chroms, i);
        if (!cJSON_IsObject(item)) continue;
        
        cJSON *fasta = cJSON_GetObjectItem(item, "fasta");
        cJSON *bam = cJSON_GetObjectItem(item, "bam");
        cJSON *name = cJSON_GetObjectItem(item, "name");
        
        if (!fasta || !cJSON_IsString(fasta) || !fasta->valuestring ||
            !bam || !cJSON_IsString(bam) || !bam->valuestring ||
            !name || !cJSON_IsString(name) || !name->valuestring)
            continue;
        
        if (strlen(fasta->valuestring) >= 64 ||
            strlen(bam->valuestring) >= 64 ||
            strlen(name->valuestring) >= 64)
            continue;
        
        valid_count++;
    }
    
    if (valid_count == 0) 
    {
        fprintf(stderr, "Error: No valid chromosome mappings found\n");
        cJSON_Delete(json);
        return -6;
    }
    
    // Allocate space for valid entries only
    *entries = calloc(valid_count, sizeof(ChromMapEntry));
    if (!*entries) 
    {
        fprintf(stderr, "Failed to allocate memory for chromosome entries\n");
        cJSON_Delete(json);
        return -5;
    }
    
    // Second pass: copy valid entries
    *n_entries = 0;
    for (int i = 0; i < count; ++i) {
        cJSON *item = cJSON_GetArrayItem(chroms, i);
        if (!cJSON_IsObject(item)) continue;
        
        cJSON *fasta = cJSON_GetObjectItem(item, "fasta");
        cJSON *bam = cJSON_GetObjectItem(item, "bam");
        cJSON *name = cJSON_GetObjectItem(item, "name");
        cJSON *extract = cJSON_GetObjectItem(item, "extract");
        
        if (!fasta || !cJSON_IsString(fasta) || !fasta->valuestring ||
            !bam || !cJSON_IsString(bam) || !bam->valuestring ||
            !name || !cJSON_IsString(name) || !name->valuestring)
            continue;
        
        if (strlen(fasta->valuestring) >= 64 ||
            strlen(bam->valuestring) >= 64 ||
            strlen(name->valuestring) >= 64)
            continue;
        
        ChromMapEntry *e = &(*entries)[*n_entries];
        
        strcpy(e->fasta, fasta->valuestring);
        strcpy(e->bam, bam->valuestring);
        strcpy(e->name, name->valuestring);
        
        if (extract && cJSON_IsBool(extract))
            e->extract = cJSON_IsTrue(extract);
        
        (*n_entries)++;
    }
    
    cJSON_Delete(json);
    
    fprintf(stderr, "Successfully loaded %d valid chromosome mappings\n", *n_entries);
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

    // Parse command line arguments
    const char *bam_file = NULL;
    const char *ref_file = NULL;
    const char *out_dir = NULL;
    const char *chrom_mapping_file = NULL;
    int min_mapq = DEFAULT_MIN_MAPQ;
    int min_phred = DEFAULT_MIN_PHRED;
    int min_cov = DEFAULT_MIN_COV;
    int min_meth = DEFAULT_MIN_METH;
    int max_meth = DEFAULT_MAX_METH;
    uint32_t chunk_size = DEFAULT_CHUNK_SIZE;
    int num_threads = DEFAULT_THREADS;
    OutputFormat output_format = OUTPUT_HDF5;  // Default to HDF5 output
    int split_context_files = 0;
    int use_gpu = 0;

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
            use_gpu = 1;
            break;
        case 'H':
            use_gpu = 1;
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
            use_gpu = 1;
            break;
        case 'k':
            chunk_size = atoi(optarg);
            if (chunk_size < 1)
            {
                fprintf(stderr, "Chunk size must be positive\n");
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
            use_gpu = 1;
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
    char *json_ref_file = NULL;
    if (!chrom_mapping_file || load_chrom_mapping(chrom_mapping_file, &chroms, &n_chroms, &json_ref_file) != 0) 
    {
        fprintf(stderr, "Failed to load chromosome mapping from %s\n", chrom_mapping_file);
        return 1;
    }

    // Use reference file from command line if provided, otherwise use the one from JSON
    if (!ref_file && json_ref_file)
    {
        ref_file = json_ref_file;
        fprintf(stderr, "Using reference file from chromosome mapping: %s\n", ref_file);
    }
    else if (ref_file)
    {
        fprintf(stderr, "Using reference file from command line: %s\n", ref_file);
    }

    if (!ref_file)
    {
        fprintf(stderr, "Reference file not provided\n");
        return 1;
    }

    // Open BAM file
    samFile *in = sam_open(bam_file, "r");
    if (!in) 
    {
        fprintf(stderr, "Failed to open BAM file: %s\n", bam_file);
        return 1;
    }

    // Load BAM header
    bam_hdr_t *header = sam_hdr_read(in);
    if (!header) 
    {
        fprintf(stderr, "Failed to read BAM header\n");
        sam_close(in);
        return 1;
    }

    // Create index for BAM file
    hts_idx_t *idx = sam_index_load(in, bam_file);
    if (!idx) 
    {
        fprintf(stderr, "Failed to load BAM index\n");
        sam_hdr_destroy(header);
        sam_close(in);
        return 1;
    }

    // Initialize processing queue
    ProcessingQueue queue = {0};
    queue.results = malloc(n_chroms * sizeof(ChromosomeResult));
    if (!queue.results) 
    {
        fprintf(stderr, "Failed to allocate memory for chromosome results\n");
        hts_idx_destroy(idx);
        sam_hdr_destroy(header);
        sam_close(in);
        return 1;
    }
    queue.total_chromosomes = n_chroms;
    pthread_mutex_init(&queue.queue_mutex, NULL);
    pthread_cond_init(&queue.queue_cond, NULL);

    // Initialize each ChromosomeResult structure
    for (int i = 0; i < n_chroms; i++) {
        queue.results[i].chr_name = NULL;
        queue.results[i].tid = -1;
        queue.results[i].context_buffers = NULL;
        queue.results[i].context_sizes = NULL;
        queue.results[i].context_pos_maps = NULL;
        queue.results[i].processed = 0;
        pthread_mutex_init(&queue.results[i].mutex, NULL);
    }

    // Initialize thread contexts and create threads
    ThreadContext *contexts = malloc(num_threads * sizeof(ThreadContext));
    if (!contexts) {
        fprintf(stderr, "Failed to allocate memory for thread contexts\n");
        return 1;
    }

    pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
    if (!threads) {
        fprintf(stderr, "Failed to allocate memory for threads\n");
        free(contexts);
        return 1;
    }

    fprintf(stderr, "Processing %d chromosomes using %d threads\n", n_chroms, num_threads);

    // Initialize thread contexts and create threads
    int thread_count = 0;
    for (int i = 0; i < n_chroms; i++) {
        contexts[thread_count].queue = &queue;
        contexts[thread_count].bam_file = bam_file;
        contexts[thread_count].out_dir = out_dir;
        contexts[thread_count].ref_file = ref_file;
        contexts[thread_count].min_mapq = min_mapq;
        contexts[thread_count].min_phred = min_phred;
        contexts[thread_count].use_gpu = use_gpu;
        contexts[thread_count].split_context_files = split_context_files;
        contexts[thread_count].output_format = output_format;
        contexts[thread_count].min_cov = min_cov;
        contexts[thread_count].cap_cov = 0;  // Default to no cap
        contexts[thread_count].min_meth = min_meth;
        contexts[thread_count].max_meth = max_meth;
        contexts[thread_count].hdf5_compression = 0;  // Default to no compression
        contexts[thread_count].hdf5_chunk_size = 0;  // Default to no chunk size
        contexts[thread_count].thread_id = thread_count;  // Set thread ID
        contexts[thread_count].keep_chg = 1;  // Always process CHG
        contexts[thread_count].keep_chh = 1;  // Always process CHH
        contexts[thread_count].chr_name = chroms[i].bam;
        contexts[thread_count].fasta_chr_name = chroms[i].fasta;
        contexts[thread_count].display_name = chroms[i].name;

        // Initialize mutex
        if (pthread_mutex_init(&contexts[thread_count].fai_mutex, NULL) != 0) {
            fprintf(stderr, "Failed to initialize mutex for thread %d\n", thread_count);
            // Clean up already created threads and mutexes
            for (int j = 0; j < thread_count; j++) {
                pthread_join(threads[j], NULL);
                pthread_mutex_destroy(&contexts[j].fai_mutex);
            }
            free(threads);
            free(contexts);
            return 1;
        }

        if (pthread_create(&threads[thread_count], NULL, process_chromosome_thread, &contexts[thread_count]) != 0) {
            fprintf(stderr, "Failed to create thread %d\n", thread_count);
            // Clean up already created threads and mutexes
            for (int j = 0; j < thread_count; j++) {
                pthread_join(threads[j], NULL);
                pthread_mutex_destroy(&contexts[j].fai_mutex);
            }
            pthread_mutex_destroy(&contexts[thread_count].fai_mutex);
            free(threads);
            free(contexts);
            return 1;
        }
        thread_count++;
    }

    // Wait for all threads to complete
    for (int i = 0; i < thread_count; i++) {
        pthread_join(threads[i], NULL);
        pthread_mutex_destroy(&contexts[i].fai_mutex);
    }

    // Cleanup
    free(threads);
    free(contexts);
    free(queue.results);
    hts_idx_destroy(idx);
    sam_hdr_destroy(header);
    sam_close(in);

    return 0;
}

static void *process_chromosome_thread(void *arg) 
{
    ThreadContext *ctx = (ThreadContext *)arg;

    // Skip processing if no chromosome assigned
    if (!ctx->chr_name) {
        fprintf(stderr, "Thread %d: No chromosome assigned, skipping\n", ctx->thread_id);
        return NULL;
    }

    int len = 0;  // For faidx_fetch_seq

    //fprintf(stderr, "Thread %d starting: chr_name=%s, fasta_chr_name=%s\n", 
    //        ctx->thread_id, ctx->chr_name, ctx->fasta_chr_name);

    // Open BAM file
    ctx->fp = sam_open(ctx->bam_file, "r");
    if (!ctx->fp) 
    {
        fprintf(stderr, "Thread %d: Failed to open BAM file: %s\n", ctx->thread_id, ctx->bam_file);
        return NULL;
    }

    // Load BAM header
    ctx->header = sam_hdr_read(ctx->fp);
    if (!ctx->header) 
    {
        fprintf(stderr, "Thread %d: Failed to read BAM header\n", ctx->thread_id);
        sam_close(ctx->fp);
        return NULL;
    }

    //fprintf(stderr, "Thread %d: BAM header loaded, looking for chromosome %s\n", 
    //        ctx->thread_id, ctx->chr_name);

    // Get chromosome ID and length
    ctx->tid = sam_hdr_name2tid(ctx->header, ctx->chr_name);
    if (ctx->tid < 0) 
    {
        fprintf(stderr, "Thread %d: Chromosome %s not found in BAM header\n", ctx->thread_id, ctx->chr_name);
        sam_hdr_destroy(ctx->header);
        sam_close(ctx->fp);
        return NULL;
    }

    ctx->chr_len = ctx->header->target_len[ctx->tid];
    //fprintf(stderr, "Thread %d: Processing chromosome %s (length: %u)\n", ctx->thread_id, ctx->chr_name, ctx->chr_len);

    // Load FASTA index
    pthread_mutex_lock(&ctx->fai_mutex);
    ctx->fai = fai_load(ctx->ref_file);
    if (!ctx->fai) 
    {
        fprintf(stderr, "Thread %d: Failed to load FASTA index\n", ctx->thread_id);
        pthread_mutex_unlock(&ctx->fai_mutex);
        sam_hdr_destroy(ctx->header);
        sam_close(ctx->fp);
        return NULL;
    }

    // Fetch chromosome sequence with correct type
    ctx->chr_seq = faidx_fetch_seq(ctx->fai, ctx->fasta_chr_name, 0, ctx->chr_len - 1, &len);
    if (len > 0)
        ctx->chr_len = (uint32_t)len;  // Convert to uint32_t
    else 
    {
        fprintf(stderr, "Thread %d: Failed to fetch sequence for chromosome %s\n", ctx->thread_id, ctx->fasta_chr_name);
        fai_destroy(ctx->fai);
        sam_hdr_destroy(ctx->header);
        sam_close(ctx->fp);
        return NULL;
    }
    pthread_mutex_unlock(&ctx->fai_mutex);

    // Create BAM index
    ctx->idx = sam_index_load(ctx->fp, ctx->bam_file);
    if (!ctx->idx) 
    {
        fprintf(stderr, "Thread %d: Failed to load BAM index\n", ctx->thread_id);
        free(ctx->chr_seq);
        fai_destroy(ctx->fai);
        sam_hdr_destroy(ctx->header);
        sam_close(ctx->fp);
        return NULL;
    }

    // Initialize chromosome result
    pthread_mutex_lock(&ctx->queue->queue_mutex);
    ctx->result = &ctx->queue->results[ctx->queue->processed_count];
    ctx->result->chr_name = strdup(ctx->chr_name);
    ctx->result->tid = ctx->tid;
    ctx->result->context_buffers = malloc(3 * sizeof(MethylRecord *));
    ctx->result->context_sizes = malloc(3 * sizeof(size_t));
    ctx->result->context_pos_maps = malloc(3 * sizeof(khash_t(pos) *));
    if (!ctx->result->chr_name || !ctx->result->context_buffers || !ctx->result->context_sizes || !ctx->result->context_pos_maps) 
    {
        fprintf(stderr, "Thread %d: Failed to allocate memory for chromosome result\n", ctx->thread_id);
        pthread_mutex_unlock(&ctx->queue->queue_mutex);
        hts_idx_destroy(ctx->idx);
        free(ctx->chr_seq);
        fai_destroy(ctx->fai);
        sam_hdr_destroy(ctx->header);
        sam_close(ctx->fp);
        return NULL;
    }

    // Initialize context buffers and position maps
    for (int i = 0; i < 3; i++) 
    {
        ctx->result->context_buffers[i] = NULL;
        ctx->result->context_sizes[i] = 0;
        ctx->result->context_pos_maps[i] = kh_init(pos);
        if (!ctx->result->context_pos_maps[i]) 
        {
            fprintf(stderr, "Thread %d: Failed to initialize position map for context %d\n", ctx->thread_id, i);
            // Clean up already allocated resources
            for (int j = 0; j < i; j++) 
            {
                if (ctx->result->context_pos_maps[j]) {
                    kh_destroy(pos, ctx->result->context_pos_maps[j]);
                }
            }
            free(ctx->result->context_pos_maps);
            free(ctx->result->context_sizes);
            free(ctx->result->context_buffers);
            free(ctx->result->chr_name);
            pthread_mutex_unlock(&ctx->queue->queue_mutex);
            hts_idx_destroy(ctx->idx);
            free(ctx->chr_seq);
            fai_destroy(ctx->fai);
            sam_hdr_destroy(ctx->header);
            sam_close(ctx->fp);
            return NULL;
        }
        // Initialize the hash table with a reasonable size
        kh_resize(pos, ctx->result->context_pos_maps[i], 1024);
    }
    pthread_mutex_init(&ctx->result->mutex, NULL);
    pthread_mutex_unlock(&ctx->queue->queue_mutex);

    // Process chromosome
    process_chromosome(ctx);

    // Cleanup
    for (int i = 0; i < 3; i++) {
        if (ctx->result->context_pos_maps[i]) {
            kh_destroy(pos, ctx->result->context_pos_maps[i]);
        }
    }
    hts_idx_destroy(ctx->idx);
    free(ctx->chr_seq);
    fai_destroy(ctx->fai);
    sam_hdr_destroy(ctx->header);
    sam_close(ctx->fp);

    // Mark chromosome as processed
    pthread_mutex_lock(&ctx->queue->queue_mutex);
    ctx->result->processed = 1;
    ctx->queue->processed_count++;
    pthread_cond_broadcast(&ctx->queue->queue_cond);
    pthread_mutex_unlock(&ctx->queue->queue_mutex);

    return NULL;
}
