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
#include <errno.h>
#include <hdf5/serial/hdf5.h>
#include <pthread.h>
#include <sys/sysinfo.h>
#include <time.h>
#include "cjson/cJSON.h"

// Global debug variables and file for monitoring buffer[41].mC changes
int debug_buffer_41_changed = 0;
int debug_buffer_41_value = 0;
FILE *debug_file = NULL;

#define DEFAULT_MAX_CHR 24
#define DEFAULT_HDF5_COMPRESSION 6
#define DEFAULT_HDF5_CHUNK_SIZE 1000000
#define DEFAULT_HDF5_CHUNK_SIZE_INT 1000000
#define DEFAULT_THREADS 16
#define DEFAULT_CHUNK_SIZE 1000000
#define DEFAULT_MIN_MAPQ 30
#define DEFAULT_MIN_PHRED 20
#define DEFAULT_MIN_COV 4
#define DEFAULT_CAP_COVERAGE 1
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
#define MAX_REGIONS_PER_CHR 8
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
    int hdf5_compression;
    int hdf5_chunk_size;
    uint32_t chunk_size;
    uint32_t start_pos;
    uint32_t end_pos;
    MethylRecord *buffer;
    size_t buffer_offset;
    size_t buffer_size;
    pthread_mutex_t *buffer_mutex;
    OutputFormat output_format;  // New: replaces debug_output
    int split_context_files;
    khash_t(pos) * pos_map;
    uint8_t *tnc_array; // New: TriNucleotideContexts[25] array
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
    {
        // Create parent directories recursively
        char *parent = strdup(path);
        char *slash = strrchr(parent, '/');
        if (slash && slash != parent)
        {
            *slash = '\0';
            if (make_directory(parent) != 0)
            {
                free(parent);
                return -1;
            }
        }
        free(parent);

        // Create the directory
        return mkdir(path, 0700);
    }
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

void initialize_buffer(MethylRecord *buffer, size_t site_count, const char *chr_seq, uint32_t chr_len, int keep_chg, int keep_chh)
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

size_t find_buffer_index(MethylRecord *buffer, size_t offset, size_t size, uint32_t pos)
{
    size_t left = offset, right = offset + size - 1;
    while (left <= right)
    {
        size_t mid = left + (right - left) / 2;
        if (buffer[mid].pos == pos + 1)
            return mid;
        if (buffer[mid].pos < pos + 1)
            left = mid + 1;
        else
            right = mid - 1;
    }
    return -1;
}

// Structure to hold statistics
typedef struct {
    size_t num_positions;
    uint64_t total_methylated;
    uint64_t total_unmethylated;
    double avg_methylation_level;
    double avg_coverage;
} MethylStats;

// Function to calculate statistics from filtered buffer
static MethylStats calculate_statistics(MethylRecord *filtered_buffer, size_t n_records)
{
    MethylStats stats = {0, 0, 0, 0.0, 0.0};
    
    if (n_records == 0)
        return stats;
    
    stats.num_positions = n_records;
    
    // Use uint64_t to avoid overflow for large datasets
    uint64_t total_mC = 0;
    uint64_t total_uC = 0;
    
    for (size_t i = 0; i < n_records; i++) 
    {
        total_mC += filtered_buffer[i].mC;
        total_uC += filtered_buffer[i].uC;
    }
    
    stats.total_methylated = total_mC;
    stats.total_unmethylated = total_uC;
    
    uint64_t total_coverage = total_mC + total_uC;
    
    if (total_coverage > 0) 
    {
        stats.avg_methylation_level = (double)total_mC / (double)total_coverage;
        stats.avg_coverage = (double)total_coverage / (double)n_records;
    }
    
    return stats;
}

// Function to write statistics JSON file
static int write_statistics_json(const char *base_filename, MethylStats stats)
{
    char json_filename[1024];
    
    // Create JSON filename by replacing the extension with .json
    const char *dot = strrchr(base_filename, '.');
    if (dot && (strcmp(dot, ".txt") == 0 || strcmp(dot, ".h5") == 0)) 
    {
        size_t base_len = dot - base_filename;
        strncpy(json_filename, base_filename, base_len);
        json_filename[base_len] = '\0';
        strcat(json_filename, ".json");
    } 
    else 
        snprintf(json_filename, sizeof(json_filename), "%s.json", base_filename);
    
    // Create JSON object
    cJSON *root = cJSON_CreateObject();
    if (!root) 
    {
        fprintf(stderr, "Failed to create JSON root object\n");
        return -1;
    }
    
    // Add statistics to JSON
    cJSON_AddNumberToObject(root, "num_positions", (double)stats.num_positions);
    cJSON_AddNumberToObject(root, "total_methylated", (double)stats.total_methylated);
    cJSON_AddNumberToObject(root, "total_unmethylated", (double)stats.total_unmethylated);
    cJSON_AddNumberToObject(root, "avg_methylation_level", stats.avg_methylation_level);
    cJSON_AddNumberToObject(root, "avg_coverage", stats.avg_coverage);
    
    // Write JSON to file
    char *json_string = cJSON_Print(root);
    if (!json_string) 
    {
        fprintf(stderr, "Failed to print JSON\n");
        cJSON_Delete(root);
        return -1;
    }
    
    FILE *json_fp = fopen(json_filename, "w");
    if (!json_fp) 
    {
        fprintf(stderr, "Failed to open JSON file for writing: %s\n", json_filename);
        free(json_string);
        cJSON_Delete(root);
        return -1;
    }
    
    fprintf(json_fp, "%s\n", json_string);
    fclose(json_fp);
    
    log_time("Wrote statistics to: %s\n", json_filename);
    
    // Cleanup
    free(json_string);
    cJSON_Delete(root);
    
    return 0;
}

size_t flush_buffer(const char *filename, MethylRecord *buffer, size_t n_records,
                   int compression, int chunk_size, int append_mode,
                   int min_cov, int cap_cov,
                   OutputFormat output_format)
{
    size_t records_written = 0;
    hsize_t dims[1] = {0};
    FILE *txt_fp = NULL;
    hid_t file = -1, dataset = -1, space = -1, type = -1, mem_type = -1, dcpl = -1;
    herr_t status = -1;

    // Calculate the number of records to write (filtered by coverage)
    // and the average coverage
    double avg_cov = 0.0;
    for (size_t i = 0; i < n_records; i++)
    {
        int total = buffer[i].mC + buffer[i].uC;
        if (total >= min_cov)
        {
            dims[0]++;
            // Online mean update using dims[0] as valid_count
            avg_cov += ((double)total - avg_cov) / dims[0];
        }
    }

    // If no records pass the filter, return early
    if (dims[0] == 0)
    {
        log_time("No records passed the coverage filter\n");
        return 0;
    }

    // Prepare filtered buffer
    MethylRecord *filtered_buffer = malloc(dims[0] * sizeof(MethylRecord));
    if (!filtered_buffer)
    {
        fprintf(stderr, "Failed to allocate memory for filtered_buffer\n");
        goto cleanup;
    }

    // fill filtered buffer
    size_t j = 0;
    for (size_t i = 0; i < n_records; i++)
    {
        int total = buffer[i].mC + buffer[i].uC;
        if (total >= min_cov)
        {
            // Only cap coverage if enabled
            if (cap_cov && total > avg_cov)
            {
                double prop = (double)buffer[i].mC / total;
                uint16_t old_mC = buffer[i].mC;
                buffer[i].mC = (uint16_t)round((avg_cov * prop));
                buffer[i].uC = (uint16_t)(avg_cov - buffer[i].mC);
                if (i == 41) {
                    debug_buffer_41_changed = 1;
                    debug_buffer_41_value = buffer[41].mC;
                    if (debug_file) {
                        fprintf(debug_file,
                                "RECALCULATE: buffer[41].mC changed from %u to %u "
                                "(capping: total=%d > avg_cov=%.1f, proportion=%.6f)\n",
                                old_mC, buffer[41].mC, total, avg_cov, prop);
                        fflush(debug_file);
                    }
                }
            }
            filtered_buffer[j++] = buffer[i];
        }
    }

    // Handle TXT output
    if (output_format == OUTPUT_TXT || output_format == OUTPUT_BOTH)
    {
        char txt_filename[1024];
        const char *chr_num = strrchr(filename, '/');
        const char *dir_end = chr_num;
        if (chr_num)
            chr_num++;
        else
            chr_num = filename;

        if (dir_end)
        {
            snprintf(
                txt_filename,
                sizeof(txt_filename),
                "%.*s%.*s.txt",
                (int)(dir_end - filename + 1),
                filename,
                (int)(strrchr(chr_num, '.') - chr_num),
                chr_num);
        }
        else
        {
            snprintf(
                txt_filename,
                sizeof(txt_filename),
                "%.*s.txt",
                (int)(strrchr(chr_num, '.') - chr_num),
                chr_num);
        }

        log_time("Starting to write text file: %s\n", txt_filename);
        txt_fp = fopen(txt_filename, "w");
        if (!txt_fp)
        {
            fprintf(stderr, "Failed to open text file for writing: %s\n", txt_filename);
            goto cleanup;
        }

        for (size_t i = 0; i < dims[0]; i++)
        {
            char tnc_str[4];
            decode_trinucleotide(filtered_buffer[i].tnc.tnc, tnc_str);
            char strand = filtered_buffer[i].tnc.strand ? '-' : '+';
            int context = filtered_buffer[i].tnc.context;
            fprintf(
                txt_fp,
                "%u\t%c\t%u\t%u\t%s\t%s\n",
                filtered_buffer[i].pos,
                strand,
                filtered_buffer[i].mC,
                filtered_buffer[i].uC,
                get_context_string(context),
                tnc_str);
        }

        log_time("Finished writing text file: %s\n", txt_filename);
        fclose(txt_fp);
        txt_fp = NULL;
    }

    // Handle HDF5 output
    if (output_format == OUTPUT_HDF5 || output_format == OUTPUT_BOTH)
    {
        log_time("Starting to write HDF5 file: %s\n", filename);
        type = H5Tcreate(H5T_COMPOUND, sizeof(MethylRecord));
        H5Tinsert(type, "pos", HOFFSET(MethylRecord, pos), H5T_NATIVE_UINT32);
        H5Tinsert(type, "mC", HOFFSET(MethylRecord, mC), H5T_NATIVE_UINT16);
        H5Tinsert(type, "uC", HOFFSET(MethylRecord, uC), H5T_NATIVE_UINT16);
        H5Tinsert(type, "tnc", HOFFSET(MethylRecord, tnc), H5T_NATIVE_UINT8);

        file = append_mode ? 
            H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT) : 
            H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (file < 0)
        {
            fprintf(stderr, "Failed to %s HDF5 file: %s\n", append_mode ? "open" : "create", filename);
            goto cleanup;
        }

        hsize_t maxdims[1] = {H5S_UNLIMITED};
        space = H5Screate_simple(1, dims, maxdims);
        if (space < 0)
        {
            fprintf(stderr, "Failed to create dataspace\n");
            goto cleanup;
        }

        dcpl = H5Pcreate(H5P_DATASET_CREATE);
        if (dcpl < 0)
        {
            fprintf(stderr, "Failed to create dataset creation property list\n");
            goto cleanup;
        }

        hsize_t chunk_dims[1] = {(hsize_t)chunk_size};
        if (H5Pset_chunk(dcpl, 1, chunk_dims) < 0)
        {
            fprintf(stderr, "Failed to set chunking\n");
            goto cleanup;
        }

        // Set Zstandard compression
        unsigned int cd_values[1] = {compression};
        status = H5Pset_filter(dcpl, ZSTD_FILTER, H5Z_FLAG_OPTIONAL, 1, cd_values);
        if (status < 0) 
        {
            fprintf(stderr, "Failed to set Zstandard filter\n");
            goto cleanup;    
        }

        mem_type = H5Tcopy(type);
        if (mem_type < 0)
        {
            fprintf(stderr, "Failed to create memory datatype\n");
            goto cleanup;
        }

        if (append_mode)
        {
            dataset = H5Dopen2(file, "methylation_data", H5P_DEFAULT);
            if (dataset >= 0)
            {
                hsize_t curr_size;
                hid_t file_space = H5Dget_space(dataset);
                H5Sget_simple_extent_dims(file_space, &curr_size, NULL);
                dims[0] += curr_size;
                H5Dset_extent(dataset, dims);
                file_space = H5Dget_space(dataset);
                hsize_t start[1] = {curr_size};
                hsize_t count[1] = {dims[0] - curr_size};
                H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);
                status = H5Dwrite(dataset, mem_type, space, file_space, H5P_DEFAULT, filtered_buffer);
                H5Sclose(file_space);
            }
            else
            {
                dataset = H5Dcreate2(file, "methylation_data", type, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
                if (dataset >= 0)
                    status = H5Dwrite(dataset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, filtered_buffer);
            }
        }
        else
        {
            dataset = H5Dcreate2(file, "methylation_data", type, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
            if (dataset >= 0)
                status = H5Dwrite(dataset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, filtered_buffer);
        }

        if (dataset < 0)
        {
            fprintf(stderr, "Failed to create or open dataset\n");
            goto cleanup;
        }

        if (status < 0)
        {
            fprintf(stderr, "Failed to write data to HDF5 file\n");
            goto cleanup;
        }

        if (dataset >= 0)
            H5Dflush(dataset);
        if (file >= 0)
            H5Fflush(file, H5F_SCOPE_GLOBAL);

        log_time("Finished writing HDF5 file: %s\n", filename);
    }

    // Calculate and write statistics (same for both formats since content is identical)
    MethylStats stats = calculate_statistics(filtered_buffer, dims[0]);
    write_statistics_json(filename, stats);

    records_written = dims[0];
    log_time("Finished writing output files, wrote %llu filtered records\n", (unsigned long long)dims[0]);

cleanup:
    if (filtered_buffer)
        free(filtered_buffer);
    if (txt_fp)
        fclose(txt_fp);
    if (mem_type >= 0)
        H5Tclose(mem_type);
    if (dataset >= 0)
        H5Dclose(dataset);
    if (dcpl >= 0)
        H5Pclose(dcpl);
    if (space >= 0)
        H5Sclose(space);
    if (file >= 0)
    {
        H5Fflush(file, H5F_SCOPE_GLOBAL);
        H5Fclose(file);
    }
    return records_written;
}

int getRealStrand(bam1_t *b)
{
    char *XG = (char *)bam_aux_get(b, "XG");
    if (XG != NULL && *(XG + 1) != 'C' && *(XG + 1) != 'G')
        XG = NULL;
    if (XG == NULL)
    {
        if (b->core.flag & BAM_FPAIRED)
        {
            if ((b->core.flag & (BAM_FREAD1 | BAM_FREVERSE)) == (BAM_FREAD1 | BAM_FREVERSE))
                return 2;
            else if (b->core.flag & BAM_FREAD1)
                return 1;
            else if ((b->core.flag & (BAM_FREAD2 | BAM_FREVERSE)) == (BAM_FREAD2 | BAM_FREVERSE))
                return 1;
            else if (b->core.flag & BAM_FREAD2)
                return 2;
            return 0;
        }
        else
        {
            if (b->core.flag & BAM_FREVERSE)
                return 2;
            return 1;
        }
    }
    else
    {
        if (*(XG + 1) == 'C')
        {
            if ((b->core.flag & (BAM_FREAD1 | BAM_FREVERSE)) == (BAM_FREAD1 | BAM_FREVERSE))
                return 1;
            else if ((b->core.flag & BAM_FREAD1) == BAM_FREAD1)
                return 3;
            else if ((b->core.flag & (BAM_FREAD2 | BAM_FREVERSE)) == (BAM_FREAD2 | BAM_FREVERSE))
                return 3;
            else if ((b->core.flag & BAM_FREAD2) == BAM_FREAD2)
                return 1;
            else if (b->core.flag & BAM_FREVERSE)
                return 3;
            else
                return 1;
        }
        else
        {
            if ((b->core.flag & (BAM_FREAD1 | BAM_FREVERSE)) == (BAM_FREAD1 | BAM_FREVERSE))
                return 4;
            else if ((b->core.flag & BAM_FREAD1) == BAM_FREAD1)
                return 2;
            else if ((b->core.flag & (BAM_FREAD2 | BAM_FREVERSE)) == (BAM_FREAD2 | BAM_FREVERSE))
                return 2;
            else if ((b->core.flag & BAM_FREAD2) == BAM_FREAD2)
                return 4;
            else if (b->core.flag & BAM_FREVERSE)
                return 2;
            else
                return 4;
        }
    }
}

static int mplp_fetch(void *data, bam1_t *b) 
{
    int rv;
    mplp_data_t *ldata = (mplp_data_t *)data;
    uint8_t *p;

    // Debug: Check for NULL pointers
    if (!ldata) {
        //fprintf(stderr, "[DEBUG] mplp_fetch: ldata is NULL!\n");
        return -1;
    }
    if (!ldata->in) {
        //fprintf(stderr, "[DEBUG] mplp_fetch: ldata->in is NULL!\n");
        return -1;
    }
    if (!ldata->hdr) {
        //fprintf(stderr, "[DEBUG] mplp_fetch: ldata->hdr is NULL!\n");
        return -1;
    }
    if (!b) {
        //fprintf(stderr, "[DEBUG] mplp_fetch: bam1_t *b is NULL!\n");
        return -1;
    }
    if (ldata->iter == NULL) {
        //fprintf(stderr, "[DEBUG] mplp_fetch: ldata->iter is NULL (using sam_read1 fallback)\n");
        return -1;
    }
    if (ldata->targ == NULL) {
        //fprintf(stderr, "[DEBUG] mplp_fetch: ldata->targ is NULL!\n");
        return -1;
    }
    //fprintf(stderr, "[DEBUG] mplp_fetch: ldata=%p, in=%p, hdr=%p, iter=%p, targ=%p, b=%p\n", ldata, ldata->in, ldata->hdr, ldata->iter, ldata->targ, b);

    while (1)
    {
        rv = ldata->iter ? sam_itr_next(ldata->in, ldata->iter, b) : sam_read1(ldata->in, ldata->hdr, b);

        if (rv < 0) {
            //fprintf(stderr, "[DEBUG] mplp_fetch: sam_itr_next/sam_read1 returned %d (EOF or error)\n", rv);
            return rv;
        }
        if (b->core.tid == -1 || b->core.flag & BAM_FUNMAP) {
            //fprintf(stderr, "[DEBUG] mplp_fetch: skipping unmapped read (tid == -1 or BAM_FUNMAP)\n");
            continue; // Unmapped
        }
        if (b->core.qual < ldata->targ->min_mapq) {
            //fprintf(stderr, "[DEBUG] mplp_fetch: skipping read with low mapping quality (%d < %d)\n", b->core.qual, ldata->targ->min_mapq);
            continue; //-q
        }
        if (b->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY)) {
            //fprintf(stderr, "[DEBUG] mplp_fetch: skipping read with flag 0xF00 (secondary, QC fail, duplicate, supplementary)\n");
            continue; // By default: secondary alignments, QC failed, PCR duplicates, and supplemental alignments
        }
        if (b->core.flag & BAM_FDUP) {
            //fprintf(stderr, "[DEBUG] mplp_fetch: skipping duplicate read (BAM_FDUP)\n");
            continue;
        }
        p = bam_aux_get(b, "NH");
        if (p) {
            int NH = bam_aux2i(p);
            if (NH > 1) {
                //fprintf(stderr, "[DEBUG] mplp_fetch: skipping multi-mapper (NH=%d)\n", NH);
                continue; // Ignore obvious multimappers
            }
        }
        if ((b->core.flag & (BAM_FPAIRED | BAM_FMUNMAP)) == (BAM_FPAIRED | BAM_FMUNMAP)) {
            //fprintf(stderr, "[DEBUG] mplp_fetch: skipping singleton (flag & (BAM_FPAIRED | BAM_FMUNMAP) == (BAM_FPAIRED | BAM_FMUNMAP))\n");
            continue; // Singleton
        }
        if ((b->core.flag & (BAM_FPAIRED | BAM_FPROPER_PAIR)) == BAM_FPAIRED) {
            //fprintf(stderr, "[DEBUG] mplp_fetch: skipping discordant (flag & (BAM_FPAIRED | BAM_FMUNMAP) == BAM_FPAIRED)\n");
            continue; // Discordant
        }
        if ((b->core.flag & (BAM_FPAIRED | BAM_FMUNMAP)) == BAM_FPAIRED) {
            //fprintf(stderr, "[DEBUG] mplp_fetch: marking discordant pair as proper (flag & (BAM_FPAIRED | BAM_FMUNMAP) == BAM_FPAIRED)\n");
            b->core.flag |= BAM_FPROPER_PAIR; // Discordant pairs can cause double counts
        }
        // If we reach here, the read passed all filters
        //fprintf(stderr, "[DEBUG] mplp_fetch: read passed all filters (qname=%s, tid=%d, pos=%ld, flag=0x%x)\n", bam_get_qname(b), b->core.tid, (long)b->core.pos, b->core.flag);
        break;
    }
    return rv;
}

void *process_chromosome_region(void *arg)
{
    ThreadArg *targ = (ThreadArg *)arg;
    samFile *in = sam_open(targ->bam_file, "r");
    if (!in)
    {
        fprintf(stderr, "Thread %s:%u-%u: Failed to open BAM file\n", targ->chr, targ->start_pos, targ->end_pos);
        return NULL;
    }
    bam_hdr_t *header = sam_hdr_read(in);
    hts_idx_t *idx = sam_index_load(in, targ->bam_file);
    hts_itr_t *iter = sam_itr_queryi(idx, targ->tid, targ->start_pos, targ->end_pos);
    if (!iter)
    {
        log_time("Thread %s:%u-%u: Failed to create iterator\n", targ->chr, targ->start_pos, targ->end_pos);
        sam_hdr_destroy(header);
        sam_close(in);
        return NULL;
    }

    mplp_data_t *mplp_data = malloc(sizeof(mplp_data_t));
    mplp_data->targ = targ;
    mplp_data->in = in;
    mplp_data->iter = iter;
    mplp_data->hdr = header;

    bam_mplp_t mplp = bam_mplp_init(1, mplp_fetch, (void **)&mplp_data);

    bam1_t *b = bam_init1();
    int tid, n_plp;
    hts_pos_t pos;
    const bam_pileup1_t *pileup;

    while (bam_mplp64_auto(mplp, &tid, &pos, &n_plp, &pileup) > 0)
    {
        if (tid != targ->tid) 
            continue;
        if (pos < targ->start_pos || pos >= targ->end_pos)
            continue;

        uint64_t key = ((uint64_t)targ->tid << 32) | pos;
        khint_t iter_kh = kh_get(pos, targ->pos_map, key);
        if (iter_kh == kh_end(targ->pos_map))
            continue;

        char ref_base = toupper(targ->chr_seq[pos]);
        if (ref_base != 'C' && ref_base != 'G')
            continue; // Only process C or G reference sites

        size_t idx = kh_val(targ->pos_map, iter_kh);

        // Debug: Log when we process buffer[41]
        if (idx == 41 && debug_file) {
            fprintf(debug_file, "DEBUG: Processing buffer[41] at genomic position %lu\n", pos);
            fflush(debug_file);
        }

        for (int i = 0; i < n_plp; i++)
        {
            const bam_pileup1_t *p = &pileup[i];
            if (p->is_del || p->is_refskip)
                continue;

            b = p->b;
            int read_qual = b->core.qual;
            // Ignore low mapping quality reads
            if (read_qual < targ->min_mapq)
                continue;

            int read_flag = b->core.flag;
            // Ignore reads with default flags
            if (read_flag & DEFAULT_FLAGS)
                continue;

            int strand = getRealStrand(b);
            if (strand == 0)
                continue;   // Strand not determined

            int seq_idx = p->qpos;
            uint8_t *seq = bam_get_seq(b);
            uint8_t *qual = bam_get_qual(b);

            // Ignore low quality bases
            int base_qual = qual[seq_idx];
            if (base_qual < targ->min_phred)
                continue;

            int base = bam_seqi(seq, seq_idx);

            pthread_mutex_lock(targ->buffer_mutex);

            // Now, apply MethylDackel's logic:
            // CpG
            if ((ref_base == 'C') && (strand == 1 || strand == 3))
            {
                if (base == 2) // G
                {
                    targ->buffer[idx].mC++;
                    if (idx == 41) 
                    {
                        debug_buffer_41_changed = 1;
                        debug_buffer_41_value = targ->buffer[41].mC;
                        if (debug_file) 
                        {
                            fprintf(debug_file,
                                    "INCREMENT: buffer[41].mC -> %u "
                                    "(read_qual=%d, read_flag=0x%x, base_qual=%d, "
                                    "base=%d, ref_base=%c, strand=%d, "
                                    "condition: ref_base=='C' && (strand==1||strand==3) && base==2 [G])\n",
                                    targ->buffer[41].mC,
                                    read_qual,
                                    read_flag,
                                    base_qual,
                                    base,
                                    ref_base,
                                    strand
                            );
                            fflush(debug_file);
                        }
                    }
                }
                else if (base == 8) // T
                    targ->buffer[idx].uC++;
            }
            else if ((ref_base == 'G') && (strand == 2 || strand == 4))
            {
                if (base == 4) // C
                {
                    targ->buffer[idx].mC++;
                    if (idx == 41) 
                    {
                        debug_buffer_41_changed = 1;
                        debug_buffer_41_value = targ->buffer[41].mC;
                        if (debug_file) 
                        {
                            fprintf(debug_file,
                                    "INCREMENT: buffer[41].mC -> %u "
                                    "(read_qual=%d, read_flag=0x%x, base_qual=%d, "
                                    "base=%d, ref_base=%c, strand=%d, "
                                    "condition: ref_base=='G' && (strand==2||strand==4) && base==4 [C])\n",
                                    targ->buffer[41].mC,
                                    read_qual,
                                    read_flag,
                                    base_qual,
                                    base,
                                    ref_base,
                                    strand
                            );
                            fflush(debug_file);
                        }
                    }
                }
                else if (base == 1) // A
                    targ->buffer[idx].uC++;
            }
            // Otherwise, ignore
            pthread_mutex_unlock(targ->buffer_mutex);
        }
    }
    bam_destroy1(b);
    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    sam_hdr_destroy(header);
    sam_close(in);
    return NULL;
}

void process_chromosome(ThreadArg *targ)
{
    size_t site_count = count_methylation_sites(targ->chr_seq, targ->chr_len, targ->keep_chg, targ->keep_chh);
    MethylRecord *buffer = calloc(site_count, sizeof(MethylRecord));
    if (!buffer)
    {
        fprintf(stderr, "Failed to allocate buffer for %s\n", targ->chr);
        return;
    }
    initialize_buffer(buffer, site_count, targ->chr_seq, targ->chr_len, targ->keep_chg, targ->keep_chh);

    khash_t(pos) *pos_map = kh_init(pos);
    for (size_t i = 0; i < site_count; i++)
    {
        uint64_t key = ((uint64_t)targ->tid << 32) | (buffer[i].pos - 1);
        int ret;
        khint_t iter = kh_put(pos, pos_map, key, &ret);
        kh_val(pos_map, iter) = i;
    }

    uint32_t chunk_size = targ->chunk_size;
    if (chunk_size > targ->chr_len)
        chunk_size = targ->chr_len;

    int n_regions = (int)ceil((double)targ->chr_len / chunk_size);
    if (n_regions > MAX_REGIONS_PER_CHR)
    {
        n_regions = MAX_REGIONS_PER_CHR;
        chunk_size = (targ->chr_len + n_regions - 1) / n_regions;
    }
    if (n_regions < 1)
        n_regions = 1;

    // Allocate thread arguments array
    ThreadArg *region_args = malloc(n_regions * sizeof(ThreadArg));
    if (!region_args) 
    {
        fprintf(stderr, "Failed to allocate region arguments\n");
        free(buffer);
        kh_destroy(pos, pos_map);
        return;
    }

    pthread_mutex_t buffer_mutex;
    pthread_mutex_init(&buffer_mutex, NULL);

    size_t sites_per_region = (site_count + n_regions - 1) / n_regions;

    // First, validate the chromosome name
    if (!targ->chr || strlen(targ->chr) == 0) 
    {
        fprintf(stderr, "Invalid chromosome name\n");
        free(region_args);
        pthread_mutex_destroy(&buffer_mutex);
        free(buffer);
        kh_destroy(pos, pos_map);
        return;
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
            free(buffer);
            kh_destroy(pos, pos_map);
            return;
        }

        region_args[i].chr_len = targ->chr_len;
        region_args[i].chr_seq = targ->chr_seq;
        region_args[i].min_mapq = targ->min_mapq;
        region_args[i].min_phred = targ->min_phred;
        region_args[i].min_cov = targ->min_cov;
        region_args[i].cap_cov = targ->cap_cov;
        region_args[i].keep_chg = targ->keep_chg;
        region_args[i].keep_chh = targ->keep_chh;
        region_args[i].hdf5_compression = targ->hdf5_compression;
        region_args[i].hdf5_chunk_size = targ->hdf5_chunk_size;
        region_args[i].chunk_size = targ->chunk_size;
        region_args[i].output_format = targ->output_format;
        region_args[i].split_context_files = targ->split_context_files;

        // Set shared resources
        region_args[i].buffer = buffer;
        region_args[i].buffer_mutex = &buffer_mutex;
        region_args[i].pos_map = pos_map;
        region_args[i].buffer_offset = i * sites_per_region;
        region_args[i].buffer_size = (i == n_regions - 1) ? site_count - i * sites_per_region : sites_per_region;

        uint32_t start_pos = i * chunk_size;
        uint32_t end_pos = (i == n_regions - 1) ? targ->chr_len : ((i + 1) * chunk_size);

        if (end_pos > targ->chr_len)
            end_pos = targ->chr_len;

        if (i > 0 && buffer[region_args[i].buffer_offset].pos > 0)
        {
            uint32_t buffer_start = buffer[region_args[i].buffer_offset].pos - 1;
            if (buffer_start < end_pos)
                start_pos = buffer_start;
        }
        if (start_pos >= end_pos)
            start_pos = (end_pos > chunk_size) ? end_pos - chunk_size : 0;

        region_args[i].start_pos = start_pos;
        region_args[i].end_pos = end_pos;
    }

    pthread_t threads[n_regions];
    int active_threads = 0;
    int *joined = calloc(n_regions, sizeof(int));
    if (!joined)
    {
        fprintf(stderr, "Failed to allocate memory for joined array\n");
        // Clean up region arguments
        for (int i = 0; i < n_regions; i++) 
            free(region_args[i].chr);
        free(region_args);
        pthread_mutex_destroy(&buffer_mutex);
        free(buffer);
        kh_destroy(pos, pos_map);
        return;
    }

    for (int i = 0; i < n_regions; i++)
    {
        if (pthread_create(&threads[i], NULL, process_chromosome_region, &region_args[i]) != 0)
        {
            fprintf(stderr, "Failed to create thread for %s\n", region_args[i].chr);
            continue;
        }
        active_threads++;
        for (int j = 0; j <= i; j++)
            if (!joined[j] && pthread_join(threads[j], NULL) == 0)
            {
                joined[j] = 1;
                active_threads--;
            }
        
        int max_region_threads = 8;
        while (active_threads >= max_region_threads)
            for (int j = 0; j <= i; j++)
                if (!joined[j] && pthread_join(threads[j], NULL) == 0)
                {
                    joined[j] = 1;
                    active_threads--;
                }

    }
    for (int i = 0; i < n_regions; i++)
    {
        if (!joined[i] && pthread_join(threads[i], NULL) == 0)
        {
            joined[i] = 1;
            if (active_threads > 0)
                active_threads--;
            fprintf(
                stderr,
                "Final join: Completed thread %d for chromosome region %s:%u-%u, active threads: %d\n",
                i,
                region_args[i].chr,
                region_args[i].start_pos,
                region_args[i].end_pos,
                active_threads);
        }
    }
    free(joined);
    pthread_mutex_destroy(&buffer_mutex);

    if (targ->split_context_files) 
    {
        for (int ctx = CONTEXT_CPG; ctx <= CONTEXT_CHH; ++ctx) 
        {
            if ((ctx == CONTEXT_CPG) ||
                (ctx == CONTEXT_CHG && targ->keep_chg) ||
                (ctx == CONTEXT_CHH && targ->keep_chh)) 
            {
                // Filter buffer for this context
                size_t n_ctx_records = 0;
                for (size_t i = 0; i < site_count; ++i)
                    if (buffer[i].tnc.context == ctx)
                        n_ctx_records++;

                if (n_ctx_records == 0)
                    continue;

                MethylRecord *ctx_buffer = malloc(n_ctx_records * sizeof(MethylRecord));
                size_t j = 0;
                for (size_t i = 0; i < site_count; ++i)
                    if (buffer[i].tnc.context == ctx)
                        ctx_buffer[j++] = buffer[i];

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
                    ctx_buffer,
                    n_ctx_records,
                    targ->hdf5_compression,
                    targ->hdf5_chunk_size,
                    0,
                    targ->min_cov,
                    targ->cap_cov,
                    targ->output_format);

                free(ctx_buffer);
            }
        }
    } 
    else 
    {
        // Output file name - use appropriate extension based on output format
        char out_path[1024];
        const char *ext = (targ->output_format == OUTPUT_TXT) ? ".txt" : ".h5";  // For OUTPUT_BOTH, use .h5
        snprintf(out_path, sizeof(out_path), "%s/%s%s", targ->out_dir, targ->chr, ext);
        flush_buffer(
            out_path,
            buffer,
            site_count,
            targ->hdf5_compression,
            targ->hdf5_chunk_size,
            0,
            targ->min_cov,
            targ->cap_cov,
            targ->output_format);
    }

    kh_destroy(pos, pos_map);
    free(buffer);

    // Clean up region arguments
    for (int i = 0; i < n_regions; i++)
        if (region_args[i].chr) 
            free(region_args[i].chr);
    free(region_args);
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
    {
        *reference_file = strdup(ref->valuestring);
        fprintf(stderr, "Found reference file: %s\n", *reference_file);
    }
    else
        fprintf(stderr, "Warning: No reference file specified in mapping\n");

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
        cJSON *extract = cJSON_GetObjectItem(item, "extract");
        if (!extract || !cJSON_IsBool(extract) || !cJSON_IsTrue(extract)) 
            continue;
        ChromMapEntry *e = &(*entries)[*n_entries];
        
        // Initialize all strings to empty
        e->fasta[0] = '\0';
        e->bam[0] = '\0';
        e->name[0] = '\0';
        
        cJSON *fasta = cJSON_GetObjectItem(item, "fasta");
        cJSON *bam = cJSON_GetObjectItem(item, "bam");
        cJSON *name = cJSON_GetObjectItem(item, "name");
        
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
        
        // Validate that we have all required fields
        if (e->fasta[0] == '\0' || e->bam[0] == '\0' || e->name[0] == '\0') 
        {
            fprintf(stderr, "Warning: Skipping chromosome entry with missing required fields\n");
            continue;
        }
        
        e->extract = 1;
        (*n_entries)++;
    }
    cJSON_Delete(json);
    return 0;
}

static void print_usage(const char *prog) 
{
    fprintf(stderr, "Usage: %s [options] <input.bam> <output_dir> [ref.fa]\n", prog);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h, --help                Show this help message\n");
    fprintf(stderr, "  -t, --threads INT         Number of threads [%d]\n", DEFAULT_THREADS);
    fprintf(stderr, "  -q, --min-mapq INT        Minimum mapping quality [%d]\n", DEFAULT_MIN_MAPQ);
    fprintf(stderr, "  -p, --min-phred INT       Minimum base quality [%d]\n", DEFAULT_MIN_PHRED);
    fprintf(stderr, "  -c, --min-cov INT         Minimum coverage [%d]\n", DEFAULT_MIN_COV);
    fprintf(stderr, "  -C, --cap-cov INT         Cap coverage (optional, default: 0)\n");
    fprintf(stderr, "  -G, --CHG                 Process CHG context\n");
    fprintf(stderr, "  -H, --CHH                 Process CHH context\n");
    fprintf(stderr, "  -m, --chrom-mapping FILE  Chromosome mapping file [chrom_mapping.json]\n");
    fprintf(stderr, "  -z, --compression INT     HDF5 compression level [%d]\n", DEFAULT_HDF5_COMPRESSION);
    fprintf(stderr, "  -k, --chunk-size INT      HDF5 chunk size [%d]\n", DEFAULT_HDF5_CHUNK_SIZE);
    fprintf(stderr, "  -f, --output-format STR   Output format (hdf5, txt, both) [hdf5]\n");
    fprintf(stderr, "  -s, --split               Split output by context\n");
    fprintf(stderr, "  -o, --output-dir DIR      Output directory\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Note: ref.fa is optional. If provided, it will override the reference in chrom_mapping.json\n");
    fprintf(stderr, "\n");
}

int main(int argc, char *argv[])
{
    fprintf(stderr, "Program: MethylExtractor\nParameters:\n");
    for (int i = 0; i < argc; i++)
        fprintf(stderr, "  Arg %d: %s\n", i, argv[i]);


    int hdf5_compression = DEFAULT_HDF5_COMPRESSION;
    int hdf5_chunk_size = DEFAULT_HDF5_CHUNK_SIZE;
    uint32_t chunk_size = DEFAULT_CHUNK_SIZE;
    int num_threads = DEFAULT_THREADS;
    int keep_chg = 0;
    int keep_chh = 0;
    int min_mapq = DEFAULT_MIN_MAPQ;
    int min_phred = DEFAULT_MIN_PHRED;
    int min_cov = DEFAULT_MIN_COV;
    int cap_cov = 0;  // Default to 0 (no capping)
    OutputFormat output_format = OUTPUT_HDF5;  // Default to HDF5 output
    const char *out_dir = NULL;
    int split_context_files = 0;
    const char *chrom_mapping_file = NULL;  // Default to NULL
    const char *ref_file = NULL;  // Will be set from chrom_mapping or command line

    struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"threads", required_argument, 0, 't'},
        {"min-mapq", required_argument, 0, 'q'},
        {"min-phred", required_argument, 0, 'p'},
        {"min-cov", required_argument, 0, 'c'},
        {"cap-cov", required_argument, 0, 'C'},
        {"CHG", no_argument, 0, 'G'},
        {"CHH", no_argument, 0, 'H'},
        {"chrom-mapping", required_argument, 0, 'm'},
        {"compression", required_argument, 0, 'z'},
        {"chunk-size", required_argument, 0, 'k'},
        {"output-format", required_argument, 0, 'f'},
        {"split", no_argument, 0, 's'},
        {"output-dir", required_argument, 0, 'o'},
        {0, 0, 0, 0}
    };
    int opt;
    while ((opt = getopt_long(argc, argv, "ht:q:p:c:C:GHm:z:k:f:so:", long_options, NULL)) != -1)
    {
        fprintf(stderr, "DEBUG: Processing option: %c, optarg: %s\n", opt, optarg ? optarg : "(null)");
        switch (opt)
        {
        case 'h':
            print_usage(argv[0]);
            return 0;
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
                fprintf(stderr, "Minimum coverage must be non-negative\n");
                return 1;
            }
            break;
        case 'C':
            fprintf(stderr, "DEBUG: Hit 'C' case with optarg: %s\n", optarg ? optarg : "(null)");
            cap_cov = atoi(optarg);
            fprintf(stderr, "DEBUG: cap_cov = %d\n", cap_cov);
            if (cap_cov < 0)
            {
                fprintf(stderr, "Cap coverage must be non-negative\n");
                return 1;
            }
            break;
        case 'G':
            keep_chg = 1;
            break;
        case 'H':
            keep_chh = 1;
            break;
        case 'm':
            chrom_mapping_file = optarg;
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
        case 'f':
            if (strcmp(optarg, "both") == 0)
                output_format = OUTPUT_BOTH;
            else if (strcmp(optarg, "hdf5") == 0)
                output_format = OUTPUT_HDF5;
            else if (strcmp(optarg, "txt") == 0)
                output_format = OUTPUT_TXT;
            else
            {
                fprintf(stderr, "Error: Invalid output format '%s'. Must be one of: both, hdf5, txt\n", optarg);
                return 1;
            }
            break;
        case 's':
            split_context_files = 1;
            break;
        case 'o':
            out_dir = optarg;
            break;
        case '?':
        default:
            print_usage(argv[0]);
            return 1;
        }
    }

    // Check if we have the required BAM file argument
    if (optind >= argc)
    {
        fprintf(stderr, "Error: No BAM file specified\n");
        return 1;
    }

    const char *bam_file = argv[optind];  // BAM file is the first non-option argument
    const char *cmd_ref_file = NULL;      // Will be set if reference file is provided

    if (optind + 1 < argc)  // If we have a reference file
    {
        cmd_ref_file = argv[optind + 1];  // It's the second non-option argument
    }

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

    // Open debug file for monitoring buffer[41].mC changes in output directory
    if (out_dir) {
        char debug_path[1024];
        snprintf(debug_path, sizeof(debug_path), "%s/debug_buffer_41.txt", out_dir);
        debug_file = fopen(debug_path, "w");
        if (debug_file) {
            fprintf(debug_file, "Debug log for buffer[41].mC changes\n");
            fprintf(debug_file, "==================================\n");
            fprintf(debug_file, "Program started at: %s\n", __TIME__ " " __DATE__);
            fprintf(debug_file, "Output directory: %s\n", out_dir);
            fprintf(debug_file, "BAM file: %s\n", bam_file);
            fprintf(debug_file, "Threads: %d\n\n", num_threads);
            fflush(debug_file);
        }
    }

    log_time("Starting processing...\n");

    log_time("Loading chromosome mapping...\n");
    int valid_chr_count = 0;
    ChromMapEntry *chroms = NULL;
    char chrom_mapping_path[1024];
    if (chrom_mapping_file)
    {
        // Use the provided mapping file
        strncpy(chrom_mapping_path, chrom_mapping_file, sizeof(chrom_mapping_path) - 1);
        chrom_mapping_path[sizeof(chrom_mapping_path) - 1] = '\0';
    }
    else
    {
        // Use chrom_mapping.json in the same directory as the BAM file
        char *bam_dir = strdup(bam_file);
        char *last_slash = strrchr(bam_dir, '/');
        if (last_slash) {
            *(last_slash + 1) = '\0';  // Keep the trailing slash
        } else {
            bam_dir[0] = '\0';  // No directory, use current directory
        }
        snprintf(chrom_mapping_path, sizeof(chrom_mapping_path), "%schrom_mapping.json", bam_dir);
        free(bam_dir);
    }
    
    int n_chroms = 0;
    char *mapping_ref_file = NULL;
    if (load_chrom_mapping(chrom_mapping_path, &chroms, &n_chroms, &mapping_ref_file) != 0)
    {
        fprintf(stderr, "Error: Failed to load chromosome mapping from %s\n", chrom_mapping_path);
        return 1;
    }

    // If reference file was provided on command line, use it
    if (cmd_ref_file)
    {
        ref_file = cmd_ref_file;
        log_time("Using reference from command line: %s\n", ref_file);
    }
    else if (mapping_ref_file)
    {
        // use the one from chrom_mapping
        ref_file = mapping_ref_file;
        log_time("Using reference from chrom_mapping: %s\n", ref_file);
    }
    else
    {
        fprintf(stderr, "No reference file specified in chrom_mapping.json or command line\n");
        free(chroms);
        return 1;
    }

    log_time("Loading reference...\n");
    faidx_t *fai = fai_load(ref_file);
    if (!fai)
    {
        fprintf(stderr, "Failed to load reference: %s\n", ref_file);
        free(chroms);
        return 1;
    }

    // Print chromosome names from FASTA index
    // fprintf(stderr, "\nChromosomes in FASTA file:\n");
    // for (int i = 0; i < faidx_nseq(fai); i++) {
    //     const char *name = faidx_iseq(fai, i);
    //     fprintf(stderr, "  %d: %s (length: %d)\n", i, name, faidx_seq_len(fai, name));
    // }
    // fprintf(stderr, "\n");

    log_time("Loading BAM...\n");
    samFile *in = sam_open(bam_file, "r");
    if (!in)
    {
        fprintf(stderr, "Failed to open BAM: %s\n", bam_file);
        fai_destroy(fai);
        return 1;
    }

    log_time("Reading BAM header...\n");
    bam_hdr_t *header = sam_hdr_read(in);
    if (!header)
    {
        fprintf(stderr, "Failed to read BAM header\n");
        sam_close(in);
        fai_destroy(fai);
        return 1;
    }

    // Print chromosome names from BAM header
    // fprintf(stderr, "\nChromosomes in BAM file:\n");
    // for (int i = 0; i < header->n_targets; i++) {
    //     fprintf(stderr, "  %d: %s (length: %d)\n", i, header->target_name[i], header->target_len[i]);
    // }
    // fprintf(stderr, "\n");

    sam_close(in);

    log_time("Allocating thread arguments...\n");
    ThreadArg *thread_args = malloc(header->n_targets * sizeof(ThreadArg));
    if (!thread_args)
    {
        fprintf(stderr, "Failed to allocate thread arguments\n");
        sam_hdr_destroy(header);
        fai_destroy(fai);
        return 1;
    }

    log_time("Processing chromosomes...\n");
    for (int i = 0; i < n_chroms; ++i) 
    {
        ChromMapEntry *entry = &chroms[i];

        int tid = bam_name2id(header, entry->bam);
        if (tid < 0) 
        {
            fprintf(stderr, "BAM does not contain chromosome %s\n", entry->bam);
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
        thread_args[valid_chr_count].keep_chg = keep_chg;
        thread_args[valid_chr_count].keep_chh = keep_chh;
        thread_args[valid_chr_count].hdf5_compression = hdf5_compression;
        thread_args[valid_chr_count].hdf5_chunk_size = hdf5_chunk_size;
        thread_args[valid_chr_count].chunk_size = chunk_size;
        thread_args[valid_chr_count].chr_seq = seq;
        thread_args[valid_chr_count].output_format = output_format;
        thread_args[valid_chr_count].split_context_files = split_context_files;
        valid_chr_count++;
    }
    free(chroms);

    log_time("Allocating threads...\n");
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

    log_time("Creating thread argument copies...\n");
    // Create a copy of thread arguments for each thread
    ThreadArg **thread_args_copies = malloc(valid_chr_count * sizeof(ThreadArg *));
    if (!thread_args_copies) 
    {
        fprintf(stderr, "Failed to allocate thread argument copies\n");
        free(threads);
        for (int i = 0; i < valid_chr_count; i++)
            free(thread_args[i].chr_seq);
        free(thread_args);
        sam_hdr_destroy(header);
        fai_destroy(fai);
        return 1;
    }

    log_time("Creating threads...\n");
    int active_threads = 0;
    for (int i = 0; i < valid_chr_count; i++)
    {
        // Create a deep copy of the thread arguments
        thread_args_copies[i] = malloc(sizeof(ThreadArg));
        if (!thread_args_copies[i]) 
        {
            fprintf(stderr, "Failed to allocate thread argument copy %d\n", i);
            continue;
        }

        // Initialize the structure to zero
        memset(thread_args_copies[i], 0, sizeof(ThreadArg));
        
        // Make deep copies of all string fields
        thread_args_copies[i]->bam_file = strdup(thread_args[i].bam_file);
        thread_args_copies[i]->out_dir = strdup(thread_args[i].out_dir);
        thread_args_copies[i]->chr = strdup(thread_args[i].chr);
        
        // Make a deep copy of the chromosome sequence
        thread_args_copies[i]->chr_seq = malloc(thread_args[i].chr_len + 1);
        memcpy(thread_args_copies[i]->chr_seq, thread_args[i].chr_seq, thread_args[i].chr_len + 1);

        // Copy all non-pointer fields
        thread_args_copies[i]->tid = thread_args[i].tid;
        thread_args_copies[i]->chr_len = thread_args[i].chr_len;
        thread_args_copies[i]->min_mapq = thread_args[i].min_mapq;
        thread_args_copies[i]->min_phred = thread_args[i].min_phred;
        thread_args_copies[i]->min_cov = thread_args[i].min_cov;
        thread_args_copies[i]->cap_cov = thread_args[i].cap_cov;
        thread_args_copies[i]->keep_chg = thread_args[i].keep_chg;
        thread_args_copies[i]->keep_chh = thread_args[i].keep_chh;
        thread_args_copies[i]->hdf5_compression = thread_args[i].hdf5_compression;
        thread_args_copies[i]->hdf5_chunk_size = thread_args[i].hdf5_chunk_size;
        thread_args_copies[i]->chunk_size = thread_args[i].chunk_size;
        thread_args_copies[i]->output_format = thread_args[i].output_format;
        thread_args_copies[i]->split_context_files = thread_args[i].split_context_files;

        if (pthread_create(&threads[i], NULL, (void *(*)(void *))process_chromosome, thread_args_copies[i]) != 0)
        {
            fprintf(stderr, "Failed to create thread for %s\n", thread_args[i].chr);
            free(thread_args_copies[i]->chr_seq);
            free(thread_args_copies[i]->chr);
            free(thread_args_copies[i]);
            continue;
        }
        active_threads++;
        for (int j = 0; j <= i; j++)
            if (pthread_join(threads[j], NULL) == 0)
                active_threads--;

        while (active_threads >= num_threads)
            for (int j = 0; j <= i; j++)
                if (pthread_join(threads[j], NULL) == 0)
                    active_threads--;
    }

    log_time("Joining threads...\n");
    for (int i = 0; i < valid_chr_count; i++)
        if (pthread_join(threads[i], NULL) == 0)
            if (active_threads > 0)
                active_threads--;

    log_time("Cleaning up HDF5...\n");
    cleanup_hdf5();

    log_time("Cleaning up thread argument copies...\n");
    // Clean up thread argument copies
    for (int i = 0; i < valid_chr_count; i++) 
        if (thread_args_copies[i]) 
        {
            free(thread_args_copies[i]->chr_seq);
            free(thread_args_copies[i]->chr);
            free(thread_args_copies[i]);
        }

    free(thread_args_copies);

    log_time("Cleaning up original thread arguments...\n");
    // Clean up original thread arguments
    for (int i = 0; i < valid_chr_count; i++)
    {
        free(thread_args[i].chr_seq);
        free(thread_args[i].chr);  // Free the copied chromosome name
    }
    free(threads);
    free(thread_args);
    sam_hdr_destroy(header);
    fai_destroy(fai);
    if (mapping_ref_file)
        free(mapping_ref_file);

    log_time("Processing complete. MethylExtractor has finished.\n");

    // Close debug file
    if (debug_file) {
        fprintf(debug_file, "\nDebug logging complete.\n");
        fclose(debug_file);
    }

    return 0;
}
