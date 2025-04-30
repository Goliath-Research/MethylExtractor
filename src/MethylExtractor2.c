#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/stat.h>
#include <ctype.h>
#include <hdf5/serial/hdf5.h>
#include <pthread.h>
#include <errno.h>

#define MAX_PILEUP_ITERATIONS 10000000 // Safeguard against infinite pileup loops
#define DEFAULT_MAX_CHR 24
#define HDF5_COMPRESSION 6
#define HDF5_CHUNK_SIZE 1000
#define MIN_MAPQ 10
#define MIN_PHRED 5
#define MIN_COV 1
#define MAX_COV 1000
#define DEFAULT_FLAGS (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY)
#define BUFFER_SIZE 10000000 // Fixed buffer size (10M records, ~130 MB)

// Structure for methylation record
typedef struct
{
    uint32_t position;
    uint32_t methylated;
    uint32_t unmethylated;
    int8_t context; // +1,+2,+3 for '+' strand; -1,-2,-3 for '-' strand
    uint8_t tnc;    // Encoded trinucleotide
} MethylRecord;

// Structure to hold methylation counts
typedef struct
{
    uint32_t meth;
    uint32_t unmeth;
} meth_counter;

// Structure to hold pileup data
typedef struct
{
    int n;
    const bam_pileup1_t *pl;
} pileup_data;

// Structure for thread arguments
typedef struct
{
    const char *bam_file;
    const char *out_dir;
    int tid;
    const char *chr;
    uint32_t chr_len;
    char *chr_seq;
} ThreadArg;

// Valid chromosomes (without 'chr' prefix)
static const char *valid_chromosomes[] = {
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
    "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
    "21", "22", "X", "Y"};
static const int num_valid_chromosomes = 24;

// Utility functions for MethylRecord
int get_context_type(int8_t context)
{
    return abs(context); // 1=CpG, 2=CHG, 3=CHH
}

char get_strand(int8_t context)
{
    return context >= 0 ? '+' : '-';
}

void decode_trinucleotide(uint8_t tnc, char *trinucl)
{
    static const char bases[] = {'A', 'C', 'G', 'T'};
    trinucl[0] = bases[(tnc >> 4) & 3];
    trinucl[1] = bases[(tnc >> 2) & 3];
    trinucl[2] = bases[tnc & 3];
    trinucl[3] = '\0';
}

// Function to normalize chromosome name (strip 'chr' prefix)
const char *normalize_chromosome(const char *chr)
{
    if (strncmp(chr, "chr", 3) == 0)
    {
        return chr + 3;
    }
    return chr;
}

// Function to check if chromosome is valid
int is_valid_chromosome(const char *chr)
{
    const char *norm_chr = normalize_chromosome(chr);
    for (int i = 0; i < num_valid_chromosomes; i++)
    {
        if (strcmp(norm_chr, valid_chromosomes[i]) == 0)
        {
            return 1;
        }
    }
    return 0;
}

// Function to create directory
int make_directory(const char *path)
{
    struct stat st = {0};
    if (stat(path, &st) == -1)
    {
        if (mkdir(path, 0700) != 0)
        {
            fprintf(stderr, "Failed to create directory %s: %s\n", path, strerror(errno));
            return -1;
        }
    }
    return 0;
}

// Function to encode trinucleotide into uint8_t
uint8_t encode_trinucleotide(const char *trinucl)
{
    uint8_t tnc = 0;
    for (int i = 0; i < 3; i++)
    {
        uint8_t val;
        switch (toupper(trinucl[i]))
        {
        case 'A':
            val = 0;
            break;
        case 'C':
            val = 1;
            break;
        case 'G':
            val = 2;
            break;
        case 'T':
            val = 3;
            break;
        default:
            val = 0;
            break; // N or invalid → A
        }
        tnc |= (val << (2 * (2 - i))); // Pack: first<<4, second<<2, third
    }
    return tnc;
}

// Function to get sequence context from cached sequence
int get_context(const char *chr_seq, int chr_len, int pos, char *context_str, char *trinucl)
{
    if (pos < 1 || pos + 1 >= chr_len)
    {
        strcpy(context_str, "CHH");
        strcpy(trinucl, "NNN");
        return 3; // Always keep CHH
    }

    char c0 = toupper(chr_seq[pos - 1]), c1 = toupper(chr_seq[pos]), c2 = toupper(chr_seq[pos + 1]);

    trinucl[0] = c0;
    trinucl[1] = c1;
    trinucl[2] = c2;
    trinucl[3] = '\0';

    if (c1 != 'C')
    {
        strcpy(context_str, "");
        return 0;
    }

    if (c2 == 'G')
    {
        strcpy(context_str, "CpG");
        return 1;
    }
    else if (c2 == 'A' || c2 == 'C' || c2 == 'T')
    {
        if (c0 == 'C' || c0 == 'A' || c0 == 'T')
        {
            strcpy(context_str, "CHH");
            return 3;
        }
        else
        {
            strcpy(context_str, "CHG");
            return 2;
        }
    }
    else
    {
        if (c0 == 'C' || c0 == 'A' || c0 == 'T' || c0 == 'N')
        {
            strcpy(context_str, "CHH");
            return 3;
        }
        else
        {
            strcpy(context_str, "CHG");
            return 2;
        }
    }
}

// Function to process CIGAR for read position
int cigar2readpos(const bam1_t *b, int refpos, int *readpos)
{
    int rpos = 0, qpos = 0;
    uint32_t *cigar = bam_get_cigar(b);

    for (int k = 0; k < b->core.n_cigar; k++)
    {
        int op = bam_cigar_op(cigar[k]);
        int len = bam_cigar_oplen(cigar[k]);

        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
        {
            if (rpos + len > refpos)
            {
                *readpos = qpos + (refpos - rpos);
                return 0;
            }
            rpos += len;
            qpos += len;
        }
        else if (op == BAM_CDEL || op == BAM_CREF_SKIP)
        {
            rpos += len;
        }
        else if (op == BAM_CINS || op == BAM_CSOFT_CLIP)
        {
            qpos += len;
        }
    }
    return -1;
}

// Function to create or append to HDF5 file
int write_hdf5(const char *filename, MethylRecord *records, size_t n_records, const char *chr)
{
    hid_t file, dataset, space, type;
    hid_t mem_type, dcpl;
    herr_t status;

    fprintf(stderr, "Thread %s: Attempting to write HDF5 file: %s with %zu records\n", chr, filename, n_records);

    // Create compound datatype
    type = H5Tcreate(H5T_COMPOUND, sizeof(MethylRecord));
    if (type < 0)
    {
        fprintf(stderr, "Thread %s: Failed to create HDF5 compound datatype: %s\n", chr, strerror(errno));
        return 1;
    }

    H5Tinsert(type, "position", HOFFSET(MethylRecord, position), H5T_NATIVE_UINT32);
    H5Tinsert(type, "methylated", HOFFSET(MethylRecord, methylated), H5T_NATIVE_UINT32);
    H5Tinsert(type, "unmethylated", HOFFSET(MethylRecord, unmethylated), H5T_NATIVE_UINT32);
    H5Tinsert(type, "context", HOFFSET(MethylRecord, context), H5T_NATIVE_INT8);
    H5Tinsert(type, "tnc", HOFFSET(MethylRecord, tnc), H5T_NATIVE_UINT8);

    // Open or create file
    file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0)
    {
        fprintf(stderr, "Thread %s: Failed to create HDF5 file %s: %s\n", chr, filename, strerror(errno));
        H5Tclose(type);
        return 1;
    }

    // Create dataspace
    hsize_t dims[1] = {n_records};
    space = H5Screate_simple(1, dims, NULL);
    if (space < 0)
    {
        fprintf(stderr, "Thread %s: Failed to create HDF5 dataspace: %s\n", chr, strerror(errno));
        H5Fclose(file);
        H5Tclose(type);
        return 1;
    }

    // Create dataset with compression
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    if (dcpl < 0)
    {
        fprintf(stderr, "Thread %s: Failed to create HDF5 dataset creation property: %s\n", chr, strerror(errno));
        H5Sclose(space);
        H5Fclose(file);
        H5Tclose(type);
        return 1;
    }

    hsize_t chunk_dims[1] = {(hsize_t)(HDF5_CHUNK_SIZE > n_records ? n_records : HDF5_CHUNK_SIZE)};
    H5Pset_chunk(dcpl, 1, chunk_dims);
    H5Pset_deflate(dcpl, HDF5_COMPRESSION);

    dataset = H5Dcreate2(file, "methylation_data", type, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    if (dataset < 0)
    {
        fprintf(stderr, "Thread %s: Failed to create HDF5 dataset in %s: %s\n", chr, filename, strerror(errno));
        H5Pclose(dcpl);
        H5Sclose(space);
        H5Fclose(file);
        H5Tclose(type);
        return 1;
    }

    // Write data
    mem_type = H5Tcopy(type);
    status = H5Dwrite(dataset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, records);
    if (status < 0)
    {
        fprintf(stderr, "Thread %s: Failed to write data to %s: %s\n", chr, filename, strerror(errno));
    }

    // Cleanup
    H5Tclose(mem_type);
    H5Dclose(dataset);
    H5Pclose(dcpl);
    H5Sclose(space);
    H5Fclose(file);
    H5Tclose(type);

    if (status >= 0)
    {
        fprintf(stderr, "Thread %s: Successfully wrote HDF5 file: %s\n", chr, filename);
    }

    return status < 0 ? 1 : 0;
}

// Function to process pileup data
void process_pileup(pileup_data *p, meth_counter *c, int refpos, const char *chr_seq, int chr_len, MethylRecord *buffer, size_t *n_records, const char *chr)
{
    // Check if position is a cytosine
    if (refpos >= chr_len || toupper(chr_seq[refpos]) != 'C')
        return;

    char context_str[4], trinucl[4];
    int ctx = get_context(chr_seq, chr_len, refpos, context_str, trinucl);
    if (ctx == 0)
        return;

    c->meth = 0;
    c->unmeth = 0;

    for (int i = 0; i < p->n; i++)
    {
        const bam_pileup1_t *pl = p->pl + i;
        const bam1_t *b = pl->b;

        if (b->core.flag & DEFAULT_FLAGS)
            continue;
        if (b->core.qual < MIN_MAPQ)
            continue;
        if (pl->is_del || pl->is_refskip)
            continue;

        int readpos;
        if (cigar2readpos(b, refpos, &readpos) != 0)
            continue;

        uint8_t *seq = bam_get_seq(b);
        uint8_t base = bam_seqi(seq, readpos);
        uint8_t *qual = bam_get_qual(b);
        if (qual[readpos] < MIN_PHRED)
            continue;

        char bbase = seq_nt16_str[base];
        if (bbase == 'C' || bbase == 'c')
        {
            c->meth++;
        }
        else if (bbase == 'T' || bbase == 't')
        {
            c->unmeth++;
        }
    }

    int total_cov = c->meth + c->unmeth;
    if (total_cov >= MIN_COV && total_cov <= MAX_COV)
    {
        if (*n_records >= BUFFER_SIZE)
        {
            fprintf(stderr, "Thread %s: Buffer overflow at pos=%d, writing partial buffer\n", chr, refpos + 1);
            return; // Buffer full, will be handled in caller
        }

        MethylRecord *rec = &buffer[*n_records];
        rec->position = refpos + 1;
        rec->methylated = c->meth;
        rec->unmethylated = c->unmeth;
        rec->context = ctx; // Positive for '+' strand
        rec->tnc = encode_trinucleotide(trinucl);
        (*n_records)++;
        fprintf(stderr, "Thread %s: Recorded methylation at pos=%d, meth=%u, unmeth=%u, context=%d, tnc=%u\n",
                chr, refpos + 1, c->meth, c->unmeth, ctx, rec->tnc);
    }
}

// Thread function to process a single chromosome
void *process_chromosome(void *arg)
{
    ThreadArg *targ = (ThreadArg *)arg;
    fprintf(stderr, "Thread started for chromosome %s (tid=%d, len=%u)\n", targ->chr, targ->tid, targ->chr_len);

    // Open BAM
    samFile *in = sam_open(targ->bam_file, "r");
    if (!in)
    {
        fprintf(stderr, "Thread %s: Failed to open BAM: %s (%s)\n", targ->chr, targ->bam_file, strerror(errno));
        return NULL;
    }

    bam_hdr_t *header = sam_hdr_read(in);
    if (!header)
    {
        fprintf(stderr, "Thread %s: Failed to read BAM header\n", targ->chr);
        sam_close(in);
        return NULL;
    }

    hts_idx_t *idx = sam_index_load(in, targ->bam_file);
    if (!idx)
    {
        fprintf(stderr, "Thread %s: BAM index missing or corrupt. Try reindexing with 'samtools index %s'\n", targ->chr, targ->bam_file);
        sam_hdr_destroy(header);
        sam_close(in);
        return NULL;
    }

    // Create iterator
    hts_itr_t *iter = sam_itr_queryi(idx, targ->tid, 0, targ->chr_len);
    if (!iter)
    {
        fprintf(stderr, "Thread %s: Failed to create iterator\n", targ->chr);
        hts_idx_destroy(idx);
        sam_hdr_destroy(header);
        sam_close(in);
        return NULL;
    }

    // Initialize pileup
    bam_plp_t plp = bam_plp_init(NULL, NULL);
    if (!plp)
    {
        fprintf(stderr, "Thread %s: Failed to initialize pileup\n", targ->chr);
        hts_itr_destroy(iter);
        hts_idx_destroy(idx);
        sam_hdr_destroy(header);
        sam_close(in);
        return NULL;
    }
    bam_plp_set_maxcnt(plp, 8000);

    bam1_t *b = bam_init1();
    if (!b)
    {
        fprintf(stderr, "Thread %s: Failed to initialize BAM record\n", targ->chr);
        bam_plp_destroy(plp);
        hts_itr_destroy(iter);
        hts_idx_destroy(idx);
        sam_hdr_destroy(header);
        sam_close(in);
        return NULL;
    }

    // Initialize fixed buffer
    MethylRecord buffer[BUFFER_SIZE];
    size_t n_records = 0;

    // Process pileup
    pileup_data p;
    int pos;
    long iteration_count = 0;
    while (sam_itr_next(in, iter, b) >= 0)
    {
        bam_plp_push(plp, b);
        while ((p.pl = bam_plp_next(plp, &targ->tid, &pos, &p.n)) != NULL)
        {
            if (iteration_count++ >= MAX_PILEUP_ITERATIONS)
            {
                fprintf(stderr, "Thread %s: Exceeded maximum pileup iterations (%ld), breaking loop\n", targ->chr, MAX_PILEUP_ITERATIONS);
                break;
            }
            fprintf(stderr, "Thread %s: Pileup at tid=%d, pos=%d, n=%d\n", targ->chr, targ->tid, pos, p.n);
            meth_counter c;
            process_pileup(&p, &c, pos, targ->chr_seq, targ->chr_len, buffer, &n_records, targ->chr);

            // Flush buffer if full
            if (n_records >= BUFFER_SIZE)
            {
                char out_path[1024];
                snprintf(out_path, sizeof(out_path), "%s/%s.h5", targ->out_dir, normalize_chromosome(targ->chr));
                if (write_hdf5(out_path, buffer, n_records, targ->chr) != 0)
                {
                    fprintf(stderr, "Thread %s: Failed to write HDF5 file\n", targ->chr);
                }
                n_records = 0; // Reset buffer
            }
        }
        if (iteration_count >= MAX_PILEUP_ITERATIONS)
            break;
    }

    // Flush pileup buffer
    fprintf(stderr, "Thread %s: Flushing pileup buffer\n", targ->chr);
    bam_plp_push(plp, NULL);
    iteration_count = 0;
    while ((p.pl = bam_plp_next(plp, &targ->tid, &pos, &p.n)) != NULL)
    {
        if (iteration_count++ >= MAX_PILEUP_ITERATIONS)
        {
            fprintf(stderr, "Thread %s: Exceeded maximum pileup iterations (%ld) during flush, breaking loop\n", targ->chr, MAX_PILEUP_ITERATIONS);
            break;
        }
        fprintf(stderr, "Thread %s: Pileup flush at tid=%d, pos=%d, n=%d\n", targ->chr, targ->tid, pos, p.n);
        meth_counter c;
        process_pileup(&p, &c, pos, targ->chr_seq, targ->chr_len, buffer, &n_records, targ->chr);

        // Flush buffer if full
        if (n_records >= BUFFER_SIZE)
        {
            char out_path[1024];
            snprintf(out_path, sizeof(out_path), "%s/%s.h5", targ->out_dir, normalize_chromosome(targ->chr));
            if (write_hdf5(out_path, buffer, n_records, targ->chr) != 0)
            {
                fprintf(stderr, "Thread %s: Failed to write HDF5 file\n", targ->chr);
            }
            n_records = 0; // Reset buffer
        }
    }

    // Reset pileup to clear internal state
    bam_plp_reset(plp);

    // Write remaining records
    if (n_records > 0)
    {
        char out_path[1024];
        snprintf(out_path, sizeof(out_path), "%s/%s.h5", targ->out_dir, normalize_chromosome(targ->chr));
        if (write_hdf5(out_path, buffer, n_records, targ->chr) != 0)
        {
            fprintf(stderr, "Thread %s: Failed to write HDF5 file\n", targ->chr);
        }
    }
    else
    {
        fprintf(stderr, "Thread %s: No valid methylation records generated, skipping HDF5 write\n", targ->chr);
    }

    // Cleanup
    fprintf(stderr, "Thread %s: Completed processing\n", targ->chr);
    bam_destroy1(b);
    bam_plp_destroy(plp);
    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    sam_hdr_destroy(header);
    sam_close(in);
    return NULL;
}

// Main function
int main(int argc, char *argv[])
{
    int max_chr = DEFAULT_MAX_CHR;
    const char *out_dir = NULL;

    struct option long_options[] = {
        {"max-chr", required_argument, 0, 'n'},
        {"o", required_argument, 0, 'o'},
        {0, 0, 0, 0}};

    int opt;
    while ((opt = getopt_long(argc, argv, "n:o:", long_options, NULL)) != -1)
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
        case 'o':
            out_dir = optarg;
            break;
        case '?':
        default:
            fprintf(stderr, "Usage: %s [options] <ref.fa> <sorted_alignments.bam>\n", argv[0]);
            fprintf(stderr, "Options:\n");
            fprintf(stderr, "  --max-chr INT  Maximum number of chromosomes (default: %d)\n", DEFAULT_MAX_CHR);
            fprintf(stderr, "  --o DIR        Output directory for HDF5 files\n");
            return 1;
        }
    }

    if (argc - optind != 2 || !out_dir)
    {
        fprintf(stderr, "Usage: %s [options] <ref.fa> <sorted_alignments.bam>\n", argv[0]);
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  --max-chr INT  Maximum number of chromosomes (default: %d)\n", DEFAULT_MAX_CHR);
        fprintf(stderr, "  --o DIR        Output directory for HDF5 files\n");
        return 1;
    }

    const char *ref_file = argv[optind];
    const char *bam_file = argv[optind + 1];

    // Create output directory
    if (make_directory(out_dir) != 0)
    {
        fprintf(stderr, "Failed to create output directory: %s\n", out_dir);
        return 1;
    }

    // Load reference
    faidx_t *fai = fai_load(ref_file);
    if (!fai)
    {
        fprintf(stderr, "Failed to load reference: %s (%s)\n", ref_file, strerror(errno));
        return 1;
    }

    // Open BAM to get header and validate index
    samFile *in = sam_open(bam_file, "r");
    if (!in)
    {
        fprintf(stderr, "Failed to open BAM: %s (%s)\n", bam_file, strerror(errno));
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

    hts_idx_t *idx = sam_index_load(in, bam_file);
    if (!idx)
    {
        fprintf(stderr, "BAM index missing or corrupt. Try reindexing with 'samtools index %s'\n", bam_file);
        sam_close(in);
        fai_destroy(fai);
        return 1;
    }
    hts_idx_destroy(idx);
    sam_close(in);

    // Collect valid chromosomes
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
        if (!is_valid_chromosome(chr))
        {
            fprintf(stderr, "Skipping invalid chromosome: %s\n", chr);
            continue;
        }

        thread_args[valid_chr_count].bam_file = bam_file;
        thread_args[valid_chr_count].out_dir = out_dir;
        thread_args[valid_chr_count].tid = tid;
        thread_args[valid_chr_count].chr = chr;
        thread_args[valid_chr_count].chr_len = header->target_len[tid];

        // Cache chromosome sequence
        int seq_len;
        char *seq = faidx_fetch_seq(fai, chr, 0, header->target_len[tid], &seq_len);
        if (!seq || seq_len <= 0)
        {
            fprintf(stderr, "Failed to fetch sequence for %s (%s)\n", chr, strerror(errno));
            if (seq)
                free(seq);
            continue;
        }
        thread_args[valid_chr_count].chr_seq = seq;

        valid_chr_count++;
    }
    fprintf(stderr, "Found %d valid chromosomes for processing\n", valid_chr_count);

    // Create threads (one per chromosome)
    pthread_t *threads = malloc(valid_chr_count * sizeof(pthread_t));
    if (!threads)
    {
        fprintf(stderr, "Failed to allocate threads\n");
        for (int i = 0; i < valid_chr_count; i++)
        {
            free(thread_args[i].chr_seq);
        }
        free(thread_args);
        sam_hdr_destroy(header);
        fai_destroy(fai);
        return 1;
    }

    for (int i = 0; i < valid_chr_count; i++)
    {
        if (pthread_create(&threads[i], NULL, process_chromosome, &thread_args[i]) != 0)
        {
            fprintf(stderr, "Failed to create thread for %s: %s\n", thread_args[i].chr, strerror(errno));
            continue;
        }
    }

    // Join all threads
    for (int i = 0; i < valid_chr_count; i++)
    {
        pthread_join(threads[i], NULL);
    }

    // Cleanup
    for (int i = 0; i < valid_chr_count; i++)
    {
        free(thread_args[i].chr_seq);
    }
    free(threads);
    free(thread_args);
    sam_hdr_destroy(header);
    fai_destroy(fai);
    fprintf(stderr, "Program completed\n");
    return 0;
}