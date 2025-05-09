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

#define DEFAULT_MAX_CHR 24
#define DEFAULT_HDF5_COMPRESSION 6
#define DEFAULT_HDF5_CHUNK_SIZE 10000
#define DEFAULT_THREADS 16
#define DEFAULT_CHUNK_SIZE 1000000
#define DEFAULT_MIN_MAPQ 30
#define DEFAULT_MIN_PHRED 20
#define DEFAULT_MIN_COV 4
#define DEFAULT_MAX_COV 100
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
#define MAX_REGIONS_PER_CHR 8

// BAM flag constants for strand determination
#define BAM_FLAG_PAIRED 0x1
#define BAM_FLAG_FIRST_IN_PAIR 0x40
#define BAM_FLAG_SECOND_IN_PAIR 0x80
#define BAM_FLAG_REVERSE_STRAND 0x10
#define BAM_FLAG_READ1_REVERSE 0x50
#define BAM_FLAG_READ2_REVERSE 0x90

// Hash table for position-to-buffer-index mapping
KHASH_MAP_INIT_INT64(pos, size_t)

// Function prototypes
static inline char decode_nucleotide(uint8_t n);
static inline void decode_trinucleotide(uint8_t tnc, char *trinucl);
static inline const char *get_context_string(int context, char strand);

typedef struct
{
    uint32_t position;
    uint32_t methylated;
    uint32_t unmethylated;
    int8_t strand_ctx;
    uint8_t tnc;
} MethylRecord;

typedef struct
{
    const char *bam_file;
    const char *out_dir;
    int tid;
    const char *chr;
    uint32_t chr_len;
    char *chr_seq;
    int min_mapq;
    int min_phred;
    int min_cov;
    int max_cov;
    int min_meth;
    int max_meth;
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
    int debug_output;
    khash_t(pos) * pos_map;
} ThreadArg;

static const char *valid_chromosomes[] = {
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
    "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
    "21", "22", "X", "Y"};
static const int num_valid_chromosomes = 24;

const char *normalize_chromosome(const char *chr)
{
    if (strncmp(chr, "chr", 3) == 0)
        return chr + 3;
    return chr;
}

static const char *get_std_chr_name(const char *chr)
{
    const char *norm = normalize_chromosome(chr);
    for (int i = 0; i < num_valid_chromosomes; i++)
    {
        if (strcmp(norm, valid_chromosomes[i]) == 0)
            return valid_chromosomes[i];
    }
    return NULL;
}

static const char *get_bam_chr_name(const char *chr)
{
    static char bam_chr[16];
    const char *std_name = get_std_chr_name(chr);
    if (std_name == NULL)
        return NULL;
    snprintf(bam_chr, sizeof(bam_chr), "chr%s", std_name);
    return bam_chr;
}

int is_valid_chromosome(const char *chr)
{
    return get_std_chr_name(chr) != NULL;
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

static inline uint8_t encode_trinucleotide(const char *trinucl)
{
    uint8_t n1 = encode_nucleotide(trinucl[0]);
    uint8_t n2 = encode_nucleotide(trinucl[1]);
    uint8_t n3 = encode_nucleotide(trinucl[2]);
    return (n1 << 6) | (n2 << 3) | n3;
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

    if (pos < 1 || pos + 1 >= chr_len)
    {
        *tnc = encode_trinucleotide("NNN");
        *strand_ctx = encode_strand_context('+', context);
        return context;
    }

    char trinucl[4], c1, c2, c3;
    c1 = toupper(chr_seq[pos - 1]);
    c2 = toupper(chr_seq[pos]);
    c3 = toupper(chr_seq[pos + 1]);
    trinucl[0] = c1;
    trinucl[1] = c2;
    trinucl[2] = c3;
    trinucl[3] = '\0';
    *tnc = encode_trinucleotide(trinucl);

    if (c2 == 'C')
    {
        *strand_ctx = encode_strand_context('+', context);
        return context;
    }
    else if (c2 == 'G')
    {
        *strand_ctx = encode_strand_context('-', context);
        return context;
    }
    else
    {
        *strand_ctx = 0;
        return 0;
    }
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
        uint8_t tnc;
        int ctx = get_context(chr_seq, chr_len, pos, &strand_ctx, &tnc, keep_chg, keep_chh);
        if (ctx)
        {
            buffer[idx].position = pos + 1;
            buffer[idx].methylated = 0;
            buffer[idx].unmethylated = 0;
            buffer[idx].strand_ctx = strand_ctx;
            buffer[idx].tnc = tnc;
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
        if (buffer[mid].position == pos + 1)
            return mid;
        if (buffer[mid].position < pos + 1)
            left = mid + 1;
        else
            right = mid - 1;
    }
    return -1;
}

size_t flush_buffer_to_hdf5(const char *filename, MethylRecord *buffer, size_t n_records,
                            int compression, int chunk_size, int append_mode,
                            int min_cov, int max_cov, int min_meth, int max_meth,
                            int debug_output)
{
    fprintf(stderr, "Starting flush_buffer_to_hdf5 for file %s\n", filename);
    hid_t file = -1, dataset = -1, space = -1, type = -1, mem_type = -1, dcpl = -1;
    herr_t status = -1;
    size_t records_written = 0;
    type = H5Tcreate(H5T_COMPOUND, sizeof(MethylRecord));
    H5Tinsert(type, "position", HOFFSET(MethylRecord, position), H5T_NATIVE_UINT32);
    H5Tinsert(type, "methylated", HOFFSET(MethylRecord, methylated), H5T_NATIVE_UINT32);
    H5Tinsert(type, "unmethylated", HOFFSET(MethylRecord, unmethylated), H5T_NATIVE_UINT32);
    H5Tinsert(type, "strand_ctx", HOFFSET(MethylRecord, strand_ctx), H5T_NATIVE_INT8);
    H5Tinsert(type, "tnc", HOFFSET(MethylRecord, tnc), H5T_NATIVE_UINT8);
    file = append_mode ? H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT) : H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0)
    {
        fprintf(stderr, "Failed to %s HDF5 file: %s\n", append_mode ? "open" : "create", filename);
        goto cleanup;
    }
    fprintf(stderr, "HDF5 file %s opened/created\n", filename);
    hsize_t dims[1] = {0};
    FILE *debug_fp = NULL;

    for (size_t i = 0; i < n_records; i++)
    {
        int total = buffer[i].methylated + buffer[i].unmethylated;
        if (total >= min_cov && total <= max_cov)
        {
            double meth_level = total > 0 ? 100.0 * ((double)buffer[i].methylated / total) : 0.0;
            if (meth_level >= min_meth && meth_level <= max_meth)
                dims[0]++;
        }
    }

    if (debug_output)
    {
        char debug_filename[1024];
        const char *chr_num = strrchr(filename, '/');
        const char *dir_end = chr_num;
        if (chr_num)
            chr_num++;
        else
            chr_num = filename;

        if (dir_end)
        {
            snprintf(debug_filename, sizeof(debug_filename), "%.*s%.*s.txt",
                     (int)(dir_end - filename + 1), filename,
                     (int)(strrchr(chr_num, '.') - chr_num), chr_num);
        }
        else
        {
            snprintf(debug_filename, sizeof(debug_filename), "%.*s.txt",
                     (int)(strrchr(chr_num, '.') - chr_num), chr_num);
        }

        debug_fp = fopen(debug_filename, "w");
    }

    MethylRecord *filtered_buffer = malloc(dims[0] * sizeof(MethylRecord));
    size_t j = 0;

    for (size_t i = 0; i < n_records; i++)
    {
        int total = buffer[i].methylated + buffer[i].unmethylated;
        if (total >= min_cov && total <= max_cov)
        {
            double meth_level = total > 0 ? 100.0 * ((double)buffer[i].methylated / total) : 0.0;
            if (meth_level >= min_meth && meth_level <= max_meth)
            {
                if (j >= dims[0]) 
                {
                    fprintf(stderr, "Error: Attempting to write beyond allocated filtered_buffer size at index %zu\n", j);
                    break;
                }
                filtered_buffer[j++] = buffer[i];

                if (debug_fp)
                {
                    char tnc_str[4];
                    decode_trinucleotide(buffer[i].tnc, tnc_str);
                    char strand = buffer[i].strand_ctx > 0 ? '+' : '-';
                    fprintf(
                        debug_fp,
                        "%u\t%u\t%u\t%c\t%s\t%s\n",
                        buffer[i].position,
                        buffer[i].methylated,
                        buffer[i].unmethylated,
                        strand,
                        get_context_string(abs(buffer[i].strand_ctx), strand),
                        tnc_str);
                }
            }
        }
    }

    if (debug_fp)
        fclose(debug_fp);

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
    if (compression > 0)
    {
        if (H5Pset_deflate(dcpl, compression) < 0)
        {
            fprintf(stderr, "Failed to set compression\n");
            goto cleanup;
        }
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
    records_written = dims[0];
    fprintf(stderr, "Finished processing records for %s, wrote %llu filtered records\n", filename, (unsigned long long)dims[0]);
cleanup:
    if (filtered_buffer)
        free(filtered_buffer);
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
        if (b->core.flag & BAM_FLAG_PAIRED)
        {
            if ((b->core.flag & BAM_FLAG_READ1_REVERSE) == BAM_FLAG_READ1_REVERSE)
                return 2;
            else if (b->core.flag & BAM_FLAG_FIRST_IN_PAIR)
                return 1;
            else if ((b->core.flag & BAM_FLAG_READ2_REVERSE) == BAM_FLAG_READ2_REVERSE)
                return 1;
            else if (b->core.flag & BAM_FLAG_SECOND_IN_PAIR)
                return 2;
            return 0;
        }
        else
        {
            if (b->core.flag & BAM_FLAG_REVERSE_STRAND)
                return 2;
            return 1;
        }
    }
    else
    {
        if (*(XG + 1) == 'C')
        {
            if ((b->core.flag & (BAM_FLAG_FIRST_IN_PAIR | BAM_FLAG_REVERSE_STRAND)) == (BAM_FLAG_FIRST_IN_PAIR | BAM_FLAG_REVERSE_STRAND))
                return 1;
            else if ((b->core.flag & BAM_FLAG_READ1_REVERSE) == BAM_FLAG_READ1_REVERSE)
                return 3;
            else if ((b->core.flag & (BAM_FLAG_SECOND_IN_PAIR | BAM_FLAG_REVERSE_STRAND)) == (BAM_FLAG_SECOND_IN_PAIR | BAM_FLAG_REVERSE_STRAND))
                return 3;
            else if ((b->core.flag & BAM_FLAG_READ2_REVERSE) == BAM_FLAG_READ2_REVERSE)
                return 1;
            else if (b->core.flag & BAM_FLAG_REVERSE_STRAND)
                return 3;
            else
                return 1;
        }
        else
        {
            if ((b->core.flag & (BAM_FLAG_FIRST_IN_PAIR | BAM_FLAG_REVERSE_STRAND)) == (BAM_FLAG_FIRST_IN_PAIR | BAM_FLAG_REVERSE_STRAND))
                return 4;
            else if ((b->core.flag & BAM_FLAG_READ1_REVERSE) == BAM_FLAG_READ1_REVERSE)
                return 2;
            else if ((b->core.flag & (BAM_FLAG_SECOND_IN_PAIR | BAM_FLAG_REVERSE_STRAND)) == (BAM_FLAG_SECOND_IN_PAIR | BAM_FLAG_REVERSE_STRAND))
                return 2;
            else if ((b->core.flag & BAM_FLAG_READ2_REVERSE) == BAM_FLAG_READ2_REVERSE)
                return 4;
            else if (b->core.flag & BAM_FLAG_REVERSE_STRAND)
                return 2;
            else
                return 4;
        }
    }
}

void *process_chromosome_region(void *arg)
{
    ThreadArg *targ = (ThreadArg *)arg;
    fprintf(stderr, "Thread %s:%u-%u: Starting processing\n", targ->chr, targ->start_pos, targ->end_pos);
    samFile *in = sam_open(targ->bam_file, "r");
    if (!in)
    {
        fprintf(stderr, "Thread %s:%u-%u: Failed to open BAM file\n", targ->chr, targ->start_pos, targ->end_pos);
        return NULL;
    }
    fprintf(stderr, "Thread %s:%u-%u: BAM file opened\n", targ->chr, targ->start_pos, targ->end_pos);
    bam_hdr_t *header = sam_hdr_read(in);
    hts_idx_t *idx = sam_index_load(in, targ->bam_file);
    hts_itr_t *iter = sam_itr_queryi(idx, targ->tid, targ->start_pos, targ->end_pos);
    if (!iter)
    {
        fprintf(stderr, "Thread %s:%u-%u: Failed to create iterator\n", targ->chr, targ->start_pos, targ->end_pos);
        sam_hdr_destroy(header);
        sam_close(in);
        return NULL;
    }
    fprintf(stderr, "Thread %s:%u-%u: Iterator created, starting read loop\n", targ->chr, targ->start_pos, targ->end_pos);
    bam1_t *b = bam_init1();
    int read_count = 0;
    while (sam_itr_next(in, iter, b) >= 0)
    {
        read_count++;
        if (read_count % 10000 == 0)
            fprintf(stderr, "Thread %s:%u-%u: Processed %d reads\n", targ->chr, targ->start_pos, targ->end_pos, read_count);

        if (b->core.flag & DEFAULT_FLAGS || b->core.qual < targ->min_mapq)
            continue;

        uint8_t *seq = bam_get_seq(b);
        uint8_t *qual = bam_get_qual(b);
        int strand = getRealStrand(b);
        if (strand == 0)
            continue;

        uint32_t *cigar = bam_get_cigar(b);
        uint32_t pos = b->core.pos;
        int seq_idx = 0;

        for (int i = 0; i < b->core.n_cigar; i++)
        {
            int op = bam_cigar_op(cigar[i]);
            int len = bam_cigar_oplen(cigar[i]);

            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
            {
                for (int k = 0; k < len; k++)
                {
                    uint32_t refpos = pos + k;
                    if (refpos < targ->start_pos || refpos >= targ->end_pos)
                        continue;
                    if (qual[seq_idx] < targ->min_phred)
                    {
                        seq_idx++;
                        continue;
                    }

                    uint64_t key = ((uint64_t)targ->tid << 32) | refpos;
                    khint_t iter = kh_get(pos, targ->pos_map, key);
                    if (iter == kh_end(targ->pos_map))
                    {
                        seq_idx++;
                        continue;
                    }
                    size_t idx = kh_val(targ->pos_map, iter);

                    int base = bam_seqi(seq, seq_idx);
                    char cstrand = (strand == 1 || strand == 3) ? '+' : '-';
                    int ctx = abs(targ->buffer[idx].strand_ctx);
                    targ->buffer[idx].strand_ctx = encode_strand_context(cstrand, ctx);

                    pthread_mutex_lock(targ->buffer_mutex);

                    if (strand & 1)
                    {
                        if (base != 2)
                        {
                            seq_idx++;
                            pthread_mutex_unlock(targ->buffer_mutex);
                            continue;
                        }
                    }
                    else
                    {
                        if (base != 4)
                        {
                            seq_idx++;
                            pthread_mutex_unlock(targ->buffer_mutex);
                            continue;
                        }
                    }

                    if ((base == 8 && (strand & 1)) || (base == 1 && !(strand & 1)))
                        targ->buffer[idx].unmethylated++;
                    else if ((base == 2 && (strand & 1)) || (base == 4 && !(strand & 1)))
                        targ->buffer[idx].methylated++;

                    pthread_mutex_unlock(targ->buffer_mutex);
                    seq_idx++;
                }
                pos += len;
            }
            else if (op == BAM_CINS)
            {
                seq_idx += len;
            }
            else if (op == BAM_CDEL || op == BAM_CREF_SKIP)
            {
                pos += len;
            }
            else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)
            {
                if (op == BAM_CSOFT_CLIP)
                    seq_idx += len;
            }
        }
    }
    fprintf(stderr, "Thread %s:%u-%u: Finished processing %d reads\n", targ->chr, targ->start_pos, targ->end_pos, read_count);
    bam_destroy1(b);
    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    sam_hdr_destroy(header);
    sam_close(in);
    fprintf(stderr, "Thread %s:%u-%u: Completed and resources cleaned up\n", targ->chr, targ->start_pos, targ->end_pos);
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
        uint64_t key = ((uint64_t)targ->tid << 32) | (buffer[i].position - 1);
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

    ThreadArg region_args[n_regions];
    pthread_mutex_t buffer_mutex;
    pthread_mutex_init(&buffer_mutex, NULL);

    size_t sites_per_region = (site_count + n_regions - 1) / n_regions;

    for (int i = 0; i < n_regions; i++)
    {
        region_args[i] = *targ;
        region_args[i].buffer = buffer;
        region_args[i].buffer_mutex = &buffer_mutex;
        region_args[i].pos_map = pos_map;
        region_args[i].buffer_offset = i * sites_per_region;
        region_args[i].buffer_size = (i == n_regions - 1) ? site_count - i * sites_per_region : sites_per_region;

        uint32_t start_pos = i * chunk_size;
        uint32_t end_pos = (i == n_regions - 1) ? targ->chr_len : ((i + 1) * chunk_size);

        if (end_pos > targ->chr_len)
            end_pos = targ->chr_len;

        if (i > 0 && buffer[region_args[i].buffer_offset].position > 0)
        {
            uint32_t buffer_start = buffer[region_args[i].buffer_offset].position - 1;
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
    if (!joined) {
        fprintf(stderr, "Failed to allocate memory for joined array\n");
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
        fprintf(stderr, "Started thread %d for chromosome region %s:%u-%u, active threads: %d\n", i, region_args[i].chr, region_args[i].start_pos, region_args[i].end_pos, active_threads);
        // Check for completed threads to prevent hanging
        for (int j = 0; j <= i; j++)
        {
            if (!joined[j] && pthread_join(threads[j], NULL) == 0)
            {
                joined[j] = 1;
                active_threads--;
                fprintf(stderr, "Completed thread %d for chromosome region %s:%u-%u, active threads: %d\n", j, region_args[j].chr, region_args[j].start_pos, region_args[j].end_pos, active_threads);
            }
        }
        // Limit active threads if necessary
        int max_region_threads = 8; // Default limit for region threads per chromosome
        while (active_threads >= max_region_threads)
        {
            for (int j = 0; j <= i; j++)
            {
                if (!joined[j] && pthread_join(threads[j], NULL) == 0)
                {
                    joined[j] = 1;
                    active_threads--;
                    fprintf(stderr, "Completed thread %d for chromosome region %s:%u-%u, active threads: %d\n", j, region_args[j].chr, region_args[j].start_pos, region_args[j].end_pos, active_threads);
                }
            }
        }
    }
    // Ensure all threads are joined
    for (int i = 0; i < n_regions; i++)
    {
        if (!joined[i] && pthread_join(threads[i], NULL) == 0)
        {
            joined[i] = 1;
            if (active_threads > 0) active_threads--;
            fprintf(stderr, "Final join: Completed thread %d for chromosome region %s:%u-%u, active threads: %d\n", i, region_args[i].chr, region_args[i].start_pos, region_args[i].end_pos, active_threads);
        }
    }
    free(joined);
    pthread_mutex_destroy(&buffer_mutex);

    char out_path[1024];
    snprintf(out_path, sizeof(out_path), "%s/%s.h5", targ->out_dir, get_std_chr_name(targ->chr));
    flush_buffer_to_hdf5(
        out_path,
        buffer,
        site_count,
        targ->hdf5_compression,
        targ->hdf5_chunk_size,
        0,
        targ->min_cov,
        targ->max_cov,
        targ->min_meth,
        targ->max_meth,
        targ->debug_output);

    kh_destroy(pos, pos_map);
    free(buffer);
}

void cleanup_hdf5(void)
{
    H5close();
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
    uint8_t n1 = (tnc >> 6) & 0x07;
    uint8_t n2 = (tnc >> 3) & 0x07;
    uint8_t n3 = tnc & 0x07;

    trinucl[0] = decode_nucleotide(n1);
    trinucl[1] = decode_nucleotide(n2);
    trinucl[2] = decode_nucleotide(n3);
    trinucl[3] = '\0';
}

static inline const char *get_context_string(int context, char strand)
{
    if (context == CONTEXT_CPG)
        return strand == '+' ? "CpG" : "GpC";
    else if (context == CONTEXT_CHG)
        return strand == '+' ? "CHG" : "GHC";
    else if (context == CONTEXT_CHH)
        return strand == '+' ? "CHH" : "GHH";
    return "???";
}

int main(int argc, char *argv[])
{
    int max_chr = DEFAULT_MAX_CHR;
    int hdf5_compression = DEFAULT_HDF5_COMPRESSION;
    int hdf5_chunk_size = DEFAULT_HDF5_CHUNK_SIZE;
    uint32_t chunk_size = DEFAULT_CHUNK_SIZE;
    int num_threads = DEFAULT_THREADS;
    int keep_chg = 0;
    int keep_chh = 0;
    int min_mapq = DEFAULT_MIN_MAPQ;
    int min_phred = DEFAULT_MIN_PHRED;
    int min_cov = DEFAULT_MIN_COV;
    int max_cov = DEFAULT_MAX_COV;
    int min_meth = DEFAULT_MIN_METH;
    int max_meth = DEFAULT_MAX_METH;
    int debug_output = 0;
    const char *out_dir = NULL;
    struct option long_options[] = {
        {"max-chr", required_argument, 0, 'n'},
        {"o", required_argument, 0, 'o'},
        {"hdf5-compression", required_argument, 0, 'z'},
        {"hdf5-chunk-size", required_argument, 0, 'k'},
        {"@", required_argument, 0, 't'},
        {"chunk-size", required_argument, 0, 's'},
        {"CHG", no_argument, 0, 'G'},
        {"CHH", no_argument, 0, 'H'},
        {"q", required_argument, 0, 'q'},
        {"p", required_argument, 0, 'p'},
        {"c", required_argument, 0, 'c'},
        {"C", required_argument, 0, 'C'},
        {"l", required_argument, 0, 'l'},
        {"L", required_argument, 0, 'L'},
        {"debug", no_argument, 0, 'd'},
        {0, 0, 0, 0}};
    int opt;
    while ((opt = getopt_long(argc, argv, "n:o:z:k:t:GHq:p:c:C:l:L:d", long_options, NULL)) != -1)
    {
        switch (opt)
        {
        case 'n':
            max_chr = atoi(optarg);
            if (max_chr < 1)
            {
                fprintf(stderr, "Maximum chromosomes must be positive\n");
                return 1;
            }
            break;
        case 'G':
            keep_chg = 1;
            break;
        case 'H':
            keep_chh = 1;
            break;
        case 'o':
            out_dir = optarg;
            break;
        case 'z':
            hdf5_compression = atoi(optarg);
            if (hdf5_compression < 0 || hdf5_compression > 9)
            {
                fprintf(stderr, "HDF5 compression level must be 0-9\n");
                return 1;
            }
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
            max_cov = atoi(optarg);
            if (max_cov < 1)
            {
                fprintf(stderr, "Maximum coverage must be positive\n");
                return 1;
            }
            break;
        case 'l':
            min_meth = atoi(optarg);
            if (min_meth < 0 || min_meth > 100)
            {
                fprintf(stderr, "Minimum methylation level must be between 0 and 100\n");
                return 1;
            }
            break;
        case 'L':
            max_meth = atoi(optarg);
            if (max_meth < 0 || max_meth > 100)
            {
                fprintf(stderr, "Maximum methylation level must be between 0 and 100\n");
                return 1;
            }
            break;
        case 'd':
            debug_output = 1;
            break;
        case '?':
        default:
            fprintf(stderr, "Usage: %s [options] <ref.fa> <sorted_alignments.bam>\n", argv[0]);
            fprintf(stderr, "Options:\n");
            fprintf(stderr, "  --max-chr INT            Maximum number of chromosomes (default: %d)\n", DEFAULT_MAX_CHR);
            fprintf(stderr, "  --o DIR                  Output directory for HDF5 files\n");
            fprintf(stderr, "  --hdf5-compression INT   Compression level (0-9, default: %d)\n", DEFAULT_HDF5_COMPRESSION);
            fprintf(stderr, "  --hdf5-chunk-size INT    Chunk size for HDF5 datasets (default: %d)\n", DEFAULT_HDF5_CHUNK_SIZE);
            fprintf(stderr, "  --chunk-size INT         Genomic chunk size for threading (default: %d)\n", DEFAULT_CHUNK_SIZE);
            fprintf(stderr, "  --@ INT                  Number of threads to use (default: %d)\n", DEFAULT_THREADS);
            fprintf(stderr, "  --CHG                    Keep CHG context methylation data\n");
            fprintf(stderr, "  --CHH                    Keep CHH context methylation data\n");
            fprintf(stderr, "  --q INT                  Minimum mapping quality (default: %d)\n", DEFAULT_MIN_MAPQ);
            fprintf(stderr, "  --p INT                  Minimum Phred score (default: %d)\n", DEFAULT_MIN_PHRED);
            fprintf(stderr, "  --c INT                  Minimum coverage (default: %d)\n", DEFAULT_MIN_COV);
            fprintf(stderr, "  --C INT                  Maximum coverage (default: %d)\n", DEFAULT_MAX_COV);
            fprintf(stderr, "  --l INT                  Minimum methylation level (default: %d)\n", DEFAULT_MIN_METH);
            fprintf(stderr, "  --L INT                  Maximum methylation level (default: %d)\n", DEFAULT_MAX_METH);
            fprintf(stderr, "  --debug                  Enable debug output (.txt files)\n");
            return 1;
        }
    }
    if (argc - optind < 2)
    {
        fprintf(stderr, "Usage: %s [options] <ref.fa> <sorted_alignments.bam>\n", argv[0]);
        return 1;
    }
    if (min_cov > max_cov)
    {
        fprintf(stderr, "Minimum coverage (%d) must not exceed maximum coverage (%d)\n", min_cov, max_cov);
        return 1;
    }

    if (min_meth > max_meth)
    {
        fprintf(stderr, "Minimum methylation level (%d) must not exceed maximum methylation level (%d)\n", min_meth, max_meth);
        return 1;
    }

    const char *ref_file = argv[optind];
    const char *bam_file = argv[optind + 1];
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
    ThreadArg *thread_args = malloc(header->n_targets * sizeof(ThreadArg));
    if (!thread_args)
    {
        fprintf(stderr, "Failed to allocate thread arguments\n");
        sam_hdr_destroy(header);
        fai_destroy(fai);
        return 1;
    }
    int valid_chr_count = 0;
    for (int tid = 0; tid < header->n_targets && valid_chr_count < max_chr; tid++)
    {
        const char *chr = header->target_name[tid];
        const char *std_chr = get_std_chr_name(chr);
        if (!std_chr)
            continue;

        int seq_len;
        char *seq = NULL;

        seq = faidx_fetch_seq(fai, std_chr, 0, header->target_len[tid], &seq_len);

        if (!seq || seq_len <= 0)
        {
            const char *bam_chr = get_bam_chr_name(std_chr);
            if (bam_chr)
            {
                if (seq)
                    free(seq);
                seq = faidx_fetch_seq(fai, bam_chr, 0, header->target_len[tid], &seq_len);
            }
        }

        if (!seq || seq_len <= 0)
        {
            fprintf(stderr, "Failed to fetch sequence for %s (tried both %s and chr%s)\n",
                    chr, std_chr, std_chr);
            if (seq)
                free(seq);
            continue;
        }

        thread_args[valid_chr_count].bam_file = bam_file;
        thread_args[valid_chr_count].out_dir = out_dir;
        thread_args[valid_chr_count].tid = tid;
        thread_args[valid_chr_count].chr = chr;
        thread_args[valid_chr_count].chr_len = header->target_len[tid];
        thread_args[valid_chr_count].min_mapq = min_mapq;
        thread_args[valid_chr_count].min_phred = min_phred;
        thread_args[valid_chr_count].min_cov = min_cov;
        thread_args[valid_chr_count].max_cov = max_cov;
        thread_args[valid_chr_count].min_meth = min_meth;
        thread_args[valid_chr_count].max_meth = max_meth;
        thread_args[valid_chr_count].keep_chg = keep_chg;
        thread_args[valid_chr_count].keep_chh = keep_chh;
        thread_args[valid_chr_count].hdf5_compression = hdf5_compression;
        thread_args[valid_chr_count].hdf5_chunk_size = hdf5_chunk_size;
        thread_args[valid_chr_count].chunk_size = chunk_size;
        thread_args[valid_chr_count].chr_seq = seq;
        thread_args[valid_chr_count].debug_output = debug_output;
        valid_chr_count++;
    }
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
        fprintf(stderr, "Started thread %d for chromosome %s, active threads: %d\n", i, thread_args[i].chr, active_threads);
        // Check for completed threads to prevent hanging
        for (int j = 0; j <= i; j++)
        {
            if (pthread_join(threads[j], NULL) == 0)
            {
                active_threads--;
                fprintf(stderr, "Completed thread %d for chromosome %s, active threads: %d\n", j, thread_args[j].chr, active_threads);
            }
        }
        // Limit active threads if necessary
        while (active_threads >= num_threads)
        {
            for (int j = 0; j <= i; j++)
            {
                if (pthread_join(threads[j], NULL) == 0)
                {
                    active_threads--;
                    fprintf(stderr, "Completed thread %d for chromosome %s, active threads: %d\n", j, thread_args[j].chr, active_threads);
                }
            }
        }
    }
    // Ensure all threads are joined
    for (int i = 0; i < valid_chr_count; i++)
    {
        if (pthread_join(threads[i], NULL) == 0)
        {
            if (active_threads > 0) active_threads--;
            fprintf(stderr, "Final join: Completed thread %d for chromosome %s, active threads: %d\n", i, thread_args[i].chr, active_threads);
        }
    }
    cleanup_hdf5();
    for (int i = 0; i < valid_chr_count; i++)
        free(thread_args[i].chr_seq);
    free(threads);
    free(thread_args);
    sam_hdr_destroy(header);
    fai_destroy(fai);
    return 0;
}