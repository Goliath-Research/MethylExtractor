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
#include <sys/sysinfo.h>

#define DEFAULT_MAX_CHR 24
#define DEFAULT_HDF5_COMPRESSION 9
#define DEFAULT_HDF5_CHUNK_SIZE 1000
#define DEFAULT_THREADS 4
#define DEFAULT_MIN_MAPQ 30
#define DEFAULT_MIN_PHRED 20
#define DEFAULT_MIN_COV 4
#define DEFAULT_MAX_COV 100
#define DEFAULT_FLAGS (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY)
#define TNC_A 0
#define TNC_C 1
#define TNC_G 2
#define TNC_T 3
#define TNC_N 4
#define CONTEXT_CPG 1
#define CONTEXT_CHG 2
#define CONTEXT_CHH 3
#define STRAND_MASK 0x80   // bit 7
#define CONTEXT_MASK 0x03  // bits 0-1
#define MAX_CHR_NAME 2     // Enough for "22" or "X" or "Y"
#define BUFFER_SIZE 100000 // Size of the methylation record buffer

// Structure for methylation record
typedef struct
{
    uint32_t position;
    uint32_t methylated;
    uint32_t unmethylated;
    int8_t strand_ctx; // bit 7: strand (0=+, 1=-), bits 0-1: context
    uint8_t tnc;       // trinucleotide context packed into a single byte
} MethylRecord;

// Structure for thread arguments
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
    int keep_chg;
    int keep_chh;
    int hdf5_compression;
    int hdf5_chunk_size;
} ThreadArg;

// Valid chromosomes (without 'chr' prefix)
static const char *valid_chromosomes[] = {
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
    "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
    "21", "22", "X", "Y"};
static const int num_valid_chromosomes = 24;

// Function to normalize chromosome name (strip 'chr' prefix)
const char *normalize_chromosome(const char *chr)
{
    if (strncmp(chr, "chr", 3) == 0)
    {
        return chr + 3;
    }
    return chr;
}

// Function to get standardized chromosome number/name
static const char *get_std_chr_name(const char *chr)
{
    // Strip any prefix (chr, Chr, etc.)
    const char *norm = normalize_chromosome(chr);

    // Return NULL if not a valid chromosome
    for (int i = 0; i < num_valid_chromosomes; i++)
    {
        if (strcmp(norm, valid_chromosomes[i]) == 0)
        {
            return valid_chromosomes[i];
        }
    }
    return NULL;
}

// Function to check if chromosome is valid
int is_valid_chromosome(const char *chr)
{
    return get_std_chr_name(chr) != NULL;
}

// Function to create directory
int make_directory(const char *path)
{
    struct stat st = {0};
    if (stat(path, &st) == -1)
    {
        return mkdir(path, 0700);
    }
    return 0;
}
// Function to encode a nucleotide to a 3-bit number
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

// Function to encode trinucleotide context into a single byte
static inline uint8_t encode_trinucleotide(const char *trinucl)
{
    uint8_t n1 = encode_nucleotide(trinucl[0]);
    uint8_t n2 = encode_nucleotide(trinucl[1]);
    uint8_t n3 = encode_nucleotide(trinucl[2]);
    return (n1 << 6) | (n2 << 3) | n3;
}

// Function to encode strand and context into a single byte
static inline int8_t encode_strand_context(char strand, int context)
{
    return (strand == '-' ? STRAND_MASK : 0) | (context & CONTEXT_MASK);
}

// Function to get sequence context from cached sequence
int get_context(const char *chr_seq, int chr_len, int pos, int8_t *strand_ctx, uint8_t *tnc, int keep_chg, int keep_chh)
{
    if (pos < 1 || pos + 1 >= chr_len)
    {
        *tnc = encode_trinucleotide("NNN");
        *strand_ctx = encode_strand_context('+', CONTEXT_CHH);
        return keep_chh ? CONTEXT_CHH : 0;
    }

    char trinucl[4];
    trinucl[0] = toupper(chr_seq[pos - 1]);
    trinucl[1] = toupper(chr_seq[pos]);
    trinucl[2] = toupper(chr_seq[pos + 1]);
    trinucl[3] = '\0';
    *tnc = encode_trinucleotide(trinucl);

    char c0 = trinucl[0], c1 = trinucl[1], c2 = trinucl[2];

    if (c1 != 'C')
    {
        *strand_ctx = 0;
        return 0;
    }

    if (c2 == 'G')
    {
        *strand_ctx = encode_strand_context('+', CONTEXT_CPG);
        return CONTEXT_CPG;
    }
    else if (c2 == 'A' || c2 == 'C' || c2 == 'T')
    {
        if (c0 == 'C' || c0 == 'A' || c0 == 'T')
        {
            *strand_ctx = encode_strand_context('+', CONTEXT_CHH);
            return keep_chh ? CONTEXT_CHH : 0;
        }
        else
        {
            *strand_ctx = encode_strand_context('+', CONTEXT_CHG);
            return keep_chg ? CONTEXT_CHG : 0;
        }
    }
    else
    {
        if (c0 == 'C' || c0 == 'A' || c0 == 'T' || c0 == 'N')
        {
            *strand_ctx = encode_strand_context('+', CONTEXT_CHH);
            return keep_chh ? CONTEXT_CHH : 0;
        }
        else
        {
            *strand_ctx = encode_strand_context('+', CONTEXT_CHG);
            return keep_chg ? CONTEXT_CHG : 0;
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

// Function to write buffer to HDF5 and return number of records written
static size_t flush_buffer_to_hdf5(const char *filename, MethylRecord *buffer,
                                   size_t n_records, int compression, int chunk_size,
                                   int append_mode)
{
    hid_t file = -1, dataset = -1, space = -1, type = -1;
    hid_t mem_type = -1, dcpl = -1;
    herr_t status = -1;
    size_t records_written = 0;

    // Create compound datatype
    type = H5Tcreate(H5T_COMPOUND, sizeof(MethylRecord));
    if (type < 0)
    {
        fprintf(stderr, "Failed to create compound datatype\n");
        goto cleanup;
    }

    H5Tinsert(type, "position", HOFFSET(MethylRecord, position), H5T_NATIVE_UINT32);
    H5Tinsert(type, "methylated", HOFFSET(MethylRecord, methylated), H5T_NATIVE_UINT32);
    H5Tinsert(type, "unmethylated", HOFFSET(MethylRecord, unmethylated), H5T_NATIVE_UINT32);
    H5Tinsert(type, "strand_ctx", HOFFSET(MethylRecord, strand_ctx), H5T_NATIVE_INT8);
    H5Tinsert(type, "tnc", HOFFSET(MethylRecord, tnc), H5T_NATIVE_UINT8);

    // Open or create file
    file = append_mode ? H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT)
                       : H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0)
    {
        fprintf(stderr, "Failed to %s HDF5 file: %s\n",
                append_mode ? "open" : "create", filename);
        goto cleanup;
    }

    // Create dataspace
    hsize_t dims[1] = {n_records};
    hsize_t maxdims[1] = {H5S_UNLIMITED};
    space = H5Screate_simple(1, dims, maxdims);
    if (space < 0)
    {
        fprintf(stderr, "Failed to create dataspace\n");
        goto cleanup;
    }

    // Create property list for dataset creation
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    if (dcpl < 0)
    {
        fprintf(stderr, "Failed to create dataset creation property list\n");
        goto cleanup;
    }

    // Set chunking and compression
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

    // Create memory datatype
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
            // Get current size and extend
            hsize_t curr_size;
            hid_t file_space = H5Dget_space(dataset);
            if (file_space < 0)
            {
                fprintf(stderr, "Failed to get dataset space\n");
                goto cleanup;
            }

            if (H5Sget_simple_extent_dims(file_space, &curr_size, NULL) < 0)
            {
                fprintf(stderr, "Failed to get dataset dimensions\n");
                H5Sclose(file_space);
                goto cleanup;
            }

            dims[0] = curr_size + n_records;
            if (H5Dset_extent(dataset, dims) < 0)
            {
                fprintf(stderr, "Failed to extend dataset\n");
                H5Sclose(file_space);
                goto cleanup;
            }

            // Write at the end
            file_space = H5Dget_space(dataset);
            if (file_space < 0)
            {
                fprintf(stderr, "Failed to get dataset space after extension\n");
                goto cleanup;
            }

            hsize_t start[1] = {curr_size};
            hsize_t count[1] = {n_records};
            if (H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL) < 0)
            {
                fprintf(stderr, "Failed to select hyperslab\n");
                H5Sclose(file_space);
                goto cleanup;
            }

            status = H5Dwrite(dataset, mem_type, space, file_space, H5P_DEFAULT, buffer);
            H5Sclose(file_space);
        }
        else
        {
            // Dataset doesn't exist, create it
            dataset = H5Dcreate2(
                file, "methylation_data",
                type,
                space,
                H5P_DEFAULT,
                dcpl,
                H5P_DEFAULT);
            if (dataset >= 0)
                status = H5Dwrite(dataset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
        }
    }
    else
    {
        // Create new dataset
        dataset = H5Dcreate2(
            file,
            "methylation_data",
            type,
            space,
            H5P_DEFAULT,
            dcpl,
            H5P_DEFAULT);
        if (dataset >= 0)
            status = H5Dwrite(dataset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
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

    // Explicitly flush the dataset and file
    if (dataset >= 0)
        H5Dflush(dataset);

    if (file >= 0)
        H5Fflush(file, H5F_SCOPE_GLOBAL);

    records_written = n_records;

cleanup:
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
        H5Fflush(file, H5F_SCOPE_GLOBAL); // One final flush before closing
        H5Fclose(file);
    }
    if (type >= 0)
        H5Tclose(type);

    return records_written;
}

// Thread function to process a single chromosome
void *process_chromosome(void *arg)
{
    ThreadArg *targ = (ThreadArg *)arg;

    // Get standardized chromosome name
    const char *std_chr = get_std_chr_name(targ->chr);
    if (!std_chr)
    {
        fprintf(stderr, "Invalid chromosome name: %s\n", targ->chr);
        return NULL;
    }

    char out_path[1024];
    snprintf(out_path, sizeof(out_path), "%s/%s.h5", targ->out_dir, std_chr);

    // Open BAM file and get index
    samFile *in = sam_open(targ->bam_file, "r");
    if (!in)
    {
        fprintf(stderr, "Thread %s: Failed to open BAM file\n", targ->chr);
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
        fprintf(stderr, "Thread %s: BAM index missing\n", targ->chr);
        sam_hdr_destroy(header);
        sam_close(in);
        return NULL;
    }

    // Create iterator for this chromosome
    hts_itr_t *iter = sam_itr_queryi(idx, targ->tid, 0, header->target_len[targ->tid]);
    if (!iter)
    {
        fprintf(stderr, "Thread %s: Failed to create iterator\n", targ->chr);
        hts_idx_destroy(idx);
        sam_hdr_destroy(header);
        sam_close(in);
        return NULL;
    }

    // Initialize fixed-size buffers
    MethylRecord current = {0};       // Current position being processed
    MethylRecord buffer[BUFFER_SIZE]; // Buffer for valid records
    size_t buffer_count = 0;          // Number of records in buffer
    int first_write = 1;              // Flag for first write to file

    // Process alignments
    bam1_t *b = bam_init1();
    uint32_t last_pos = 0;

    while (sam_itr_next(in, iter, b) >= 0)
    {
        // Skip if not meeting quality criteria
        if (b->core.flag & DEFAULT_FLAGS)
            continue;
        if (b->core.qual < targ->min_mapq)
            continue;

        // Get read sequence and quality
        uint8_t *seq = bam_get_seq(b);
        uint8_t *qual = bam_get_qual(b);

        // Process each base in the read
        for (int i = 0; i < b->core.l_qseq; i++)
        {
            uint32_t refpos = b->core.pos + i;

            // Skip if moving backwards (shouldn't happen in sorted BAM)
            if (refpos < last_pos)
                continue;

            // Check base quality
            if (qual[i] < targ->min_phred)
                continue;

            // Get context and check if it's a cytosine we care about
            int8_t strand_ctx;
            uint8_t tnc;
            int ctx = get_context(targ->chr_seq, targ->chr_len, refpos,
                                  &strand_ctx, &tnc, targ->keep_chg, targ->keep_chh);
            if (ctx == 0)
                continue;

            // If this is a new position, process the current one
            if (current.position != refpos + 1)
            {
                // Check if current position meets coverage criteria
                if (current.position > 0)
                {
                    int total = current.methylated + current.unmethylated;
                    if (total >= targ->min_cov && total <= targ->max_cov)
                    {
                        // Add to buffer
                        buffer[buffer_count++] = current;

                        // If buffer is full, write to file
                        if (buffer_count == BUFFER_SIZE)
                        {
                            if (!flush_buffer_to_hdf5(out_path, buffer, buffer_count,
                                                      targ->hdf5_compression,
                                                      targ->hdf5_chunk_size,
                                                      !first_write))
                            {
                                fprintf(stderr, "Thread %s: Failed to write buffer to HDF5\n",
                                        targ->chr);
                                break;
                            }
                            first_write = 0;
                            buffer_count = 0;
                        }
                    }
                }

                // Start new position
                current.position = refpos + 1; // 1-based
                current.methylated = 0;
                current.unmethylated = 0;
                current.strand_ctx = strand_ctx;
                current.tnc = tnc;
            }

            // Count methylation
            char base = seq_nt16_str[bam_seqi(seq, i)];
            if (base == 'C' || base == 'c')
                current.methylated++;
            else if (base == 'T' || base == 't')
                current.unmethylated++;

            last_pos = refpos;
        }
    }

    // Process last position
    if (current.position > 0)
    {
        int total = current.methylated + current.unmethylated;
        if (total >= targ->min_cov && total <= targ->max_cov)
            buffer[buffer_count++] = current;
    }

    // Flush remaining records
    if (buffer_count > 0)
    {
        if (!flush_buffer_to_hdf5(out_path, buffer, buffer_count,
                                  targ->hdf5_compression, targ->hdf5_chunk_size,
                                  !first_write))
            fprintf(stderr, "Thread %s: Failed to write final buffer to HDF5\n", targ->chr);
    }

    // Cleanup
    bam_destroy1(b);
    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    sam_hdr_destroy(header);
    sam_close(in);
    return NULL;
}

// Function to cleanup HDF5 resources
void cleanup_hdf5(void)
{
    H5close(); // Close all remaining open HDF5 identifiers
}

// Main function
int main(int argc, char *argv[])
{
    int max_chr = DEFAULT_MAX_CHR;
    int hdf5_compression = DEFAULT_HDF5_COMPRESSION;
    int hdf5_chunk_size = DEFAULT_HDF5_CHUNK_SIZE;
    int num_threads = DEFAULT_THREADS;
    int keep_chg = 0;
    int keep_chh = 0;
    int min_mapq = DEFAULT_MIN_MAPQ;
    int min_phred = DEFAULT_MIN_PHRED;
    int min_cov = DEFAULT_MIN_COV;
    int max_cov = DEFAULT_MAX_COV;
    const char *out_dir = NULL;

    struct option long_options[] = {
        {"max-chr", required_argument, 0, 'n'},
        {"o", required_argument, 0, 'o'},
        {"hdf5-compression", required_argument, 0, 'z'},
        {"hdf5-chunk-size", required_argument, 0, 'k'},
        {"@", required_argument, 0, 't'},
        {"CHG", no_argument, 0, 'G'},
        {"CHH", no_argument, 0, 'H'},
        {"q", required_argument, 0, 'q'},
        {"p", required_argument, 0, 'p'},
        {"c", required_argument, 0, 'c'},
        {"C", required_argument, 0, 'C'},
        {0, 0, 0, 0}};

    int opt;
    while ((opt = getopt_long(argc, argv, "n:o:z:k:t:GHq:p:c:C:", long_options, NULL)) != -1)
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
            if (min_cov < 1)
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
        case '?':
        default:
            fprintf(stderr, "Usage: %s [options] <ref.fa> <sorted_alignments.bam>\n", argv[0]);
            fprintf(stderr, "Options:\n");
            fprintf(stderr, "  --max-chr INT            Maximum number of chromosomes (default: %d)\n", DEFAULT_MAX_CHR);
            fprintf(stderr, "  --o DIR                  Output directory for HDF5 files\n");
            fprintf(stderr, "  --hdf5-compression INT   Compression level (0-9, default: %d)\n", DEFAULT_HDF5_COMPRESSION);
            fprintf(stderr, "  --hdf5-chunk-size INT    Chunk size for HDF5 datasets (default: %d)\n", DEFAULT_HDF5_CHUNK_SIZE);
            fprintf(stderr, "  --@ INT                  Number of threads to use (default: %d)\n", DEFAULT_THREADS);
            fprintf(stderr, "  --CHG                    Keep CHG context methylation data (default: false)\n");
            fprintf(stderr, "  --CHH                    Keep CHH context methylation data (default: false)\n");
            fprintf(stderr, "  --q INT                  Minimum mapping quality (default: %d)\n", DEFAULT_MIN_MAPQ);
            fprintf(stderr, "  --p INT                  Minimum Phred score (default: %d)\n", DEFAULT_MIN_PHRED);
            fprintf(stderr, "  --c INT                  Minimum coverage (default: %d)\n", DEFAULT_MIN_COV);
            fprintf(stderr, "  --C INT                  Maximum coverage (default: %d)\n", DEFAULT_MAX_COV);
            return 1;
        }
    }

    if (min_cov > max_cov)
    {
        fprintf(stderr, "Minimum coverage (%d) must not exceed maximum coverage (%d)\n", min_cov, max_cov);
        return 1;
    }

    const char *ref_file = argv[optind];
    const char *bam_file = argv[optind + 1];

    // Auto-detect threads if not specified
    if (num_threads == DEFAULT_THREADS)
    {
        num_threads = sysconf(_SC_NPROCESSORS_ONLN);
        if (num_threads < 1)
            num_threads = DEFAULT_THREADS;
    }

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
        fprintf(stderr, "Failed to load reference: %s\n", ref_file);
        return 1;
    }

    // Open BAM to get header
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

    // Collect valid chromosomes and preallocate storage
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
            continue;

        thread_args[valid_chr_count].bam_file = bam_file;
        thread_args[valid_chr_count].out_dir = out_dir;
        thread_args[valid_chr_count].tid = tid;
        thread_args[valid_chr_count].chr = chr;
        thread_args[valid_chr_count].chr_len = header->target_len[tid];
        thread_args[valid_chr_count].min_mapq = min_mapq;
        thread_args[valid_chr_count].min_phred = min_phred;
        thread_args[valid_chr_count].min_cov = min_cov;
        thread_args[valid_chr_count].max_cov = max_cov;
        thread_args[valid_chr_count].keep_chg = keep_chg;
        thread_args[valid_chr_count].keep_chh = keep_chh;
        thread_args[valid_chr_count].hdf5_compression = hdf5_compression;
        thread_args[valid_chr_count].hdf5_chunk_size = hdf5_chunk_size;

        // Cache chromosome sequence
        int seq_len;
        char *seq = faidx_fetch_seq(fai, chr, 0, header->target_len[tid], &seq_len);
        if (!seq || seq_len <= 0)
        {
            fprintf(stderr, "Failed to fetch sequence for %s\n", chr);
            if (seq)
                free(seq);
            continue;
        }
        thread_args[valid_chr_count].chr_seq = seq;

        valid_chr_count++;
    }

    // Process chromosomes in parallel
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

    int active_threads = 0;
    for (int i = 0; i < valid_chr_count; i++)
    {
        if (pthread_create(&threads[i], NULL, process_chromosome, &thread_args[i]) != 0)
        {
            fprintf(stderr, "Failed to create thread for %s\n", thread_args[i].chr);
            continue;
        }
        active_threads++;
        while (active_threads >= num_threads)
        {
            for (int j = 0; j < i; j++)
            {
                if (pthread_join(threads[j], NULL) == 0)
                {
                    active_threads--;
                }
            }
        }
    }

    // Join remaining threads
    for (int i = 0; i < valid_chr_count; i++)
        pthread_join(threads[i], NULL);

    // Cleanup HDF5 resources
    cleanup_hdf5();

    // Cleanup
    for (int i = 0; i < valid_chr_count; i++)
    {
        free(thread_args[i].chr_seq);
    }
    free(threads);
    free(thread_args);
    sam_hdr_destroy(header);
    fai_destroy(fai);
    return 0;
}