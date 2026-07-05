#include "methyl_extractor.h"

static void print_usage(const char *prog)
{
    fprintf(stderr, "Usage: %s [options] <input.bam> [output_dir] [ref.fa]\n",
            prog);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h, --help                Show this help message\n");
    fprintf(stderr, "  -t, --threads INT         Number of worker threads "
                    "[auto: CPU count]\n");
    fprintf(stderr, "  -q, --min-mapq INT        Minimum mapping quality [%d]\n",
            DEFAULT_MIN_MAPQ);
    fprintf(stderr, "  -p, --min-phred INT       Minimum base quality [%d]\n",
            DEFAULT_MIN_PHRED);
    fprintf(stderr, "  -c, --min-cov INT         Minimum coverage [%d]\n",
            DEFAULT_MIN_COV);
    fprintf(stderr,
            "  -C, --cap-cov INT         Cap coverage (optional, default: 0)\n");
    fprintf(stderr, "  -G, --CHG                 Process CHG context\n");
    fprintf(stderr, "  -H, --CHH                 Process CHH context\n");
    fprintf(stderr, "  -m, --chrom-mapping FILE  Chromosome mapping file "
                    "[chrom_mapping.json]\n");
    fprintf(
        stderr,
        "  -z, --compression INT     Compression level for HDF5 [%d]\n",
        DEFAULT_HDF5_COMPRESSION);
    fprintf(stderr, "                            0=none, 1-19=Zstd (gzip fallback if Zstd unavailable)\n");
    fprintf(stderr, "  -k, --chunk-size INT      HDF5 chunk size [%d]\n",
            DEFAULT_HDF5_CHUNK_SIZE);
    fprintf(stderr, "  -f, --output-format STR   Output format (hdf5, txt, both) "
                    "[hdf5]\n");
    fprintf(stderr, "  -s, --split               Split output by context\n");
    fprintf(stderr, "  -o, --output-dir DIR      Output directory\n");
    fprintf(stderr, "  -R, --read-level          Emit read-level pattern sidecars "
                    "({chrom}-{ctx}.patterns.h5)\n");
    fprintf(stderr, "  -T, --tile-size INT       CpG sites per read-level tile "
                    "[%d]\n", DEFAULT_TILE_SIZE);
    fprintf(stderr, "\n");
    fprintf(stderr, "Note: output_dir may be given positionally (2nd argument) "
                    "or via -o/--output-dir.\n");
    fprintf(stderr, "Note: ref.fa is optional. If provided, it will override the "
                    "reference in chrom_mapping.json\n");
    fprintf(stderr, "\n");
}

int main(int argc, char *argv[])
{
    fprintf(stderr, "Program: MethylExtractor\nParameters:\n");
    for (int i = 0; i < argc; i++)
        fprintf(stderr, "  Arg %d: %s\n", i, argv[i]);

    int compression = DEFAULT_HDF5_COMPRESSION;
    int hdf5_chunk_size = DEFAULT_HDF5_CHUNK_SIZE;
    uint32_t chunk_size = DEFAULT_CHUNK_SIZE;
    int num_threads = 0; // 0 = auto-detect from CPU count
    int keep_chg = 0;
    int keep_chh = 0;
    int min_mapq = DEFAULT_MIN_MAPQ;
    int min_phred = DEFAULT_MIN_PHRED;
    int min_cov = DEFAULT_MIN_COV;
    int cap_cov = DEFAULT_CAP_COVERAGE;
    OutputFormat output_format = OUTPUT_HDF5; // Default to HDF5 output
    const char *out_dir = NULL;
    int split_context_files = 0;
    int read_level = 0;
    int tile_size = DEFAULT_TILE_SIZE;
    const char *chrom_mapping_file = NULL; // Default to NULL
    const char *ref_file = NULL;           // Will be set from chrom_mapping or command line

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
        {"read-level", no_argument, 0, 'R'},
        {"tile-size", required_argument, 0, 'T'},
        {0, 0, 0, 0}};
    int opt;
    while ((opt = getopt_long(argc, argv, "ht:q:p:c:C:GHm:z:k:f:so:RT:",
                              long_options, NULL)) != -1)
    {
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
            cap_cov = atoi(optarg);
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
            compression = atoi(optarg);
            if (compression < 0 || compression > 19)
            {
                fprintf(stderr, "Compression level must be between 0 and 19\n");
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
        case 'f':
            if (strcmp(optarg, "both") == 0)
                output_format = OUTPUT_BOTH;
            else if (strcmp(optarg, "hdf5") == 0)
                output_format = OUTPUT_HDF5;
            else if (strcmp(optarg, "txt") == 0)
                output_format = OUTPUT_TXT;
            else
            {
                fprintf(stderr,
                        "Error: Invalid output format '%s'. Must be one of: both, "
                        "hdf5, txt\n",
                        optarg);
                return 1;
            }
            break;
        case 's':
            split_context_files = 1;
            break;
        case 'o':
            out_dir = optarg;
            break;
        case 'R':
            read_level = 1;
            break;
        case 'T':
            tile_size = atoi(optarg);
            if (tile_size < MIN_TILE_SIZE || tile_size > MAX_TILE_SIZE)
            {
                fprintf(stderr, "Tile size must be between %d and %d\n",
                        MIN_TILE_SIZE, MAX_TILE_SIZE);
                return 1;
            }
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

    const char *bam_file =
        argv[optind];                // BAM file is the first non-option argument
    const char *cmd_ref_file = NULL; // Will be set if reference file is provided

    // Positional layout after options: <input.bam> [output_dir] [ref.fa].
    // If -o/--output-dir was not given, the next positional is the output dir;
    // the following positional (if any) is the reference.
    int next_pos = optind + 1;
    if (!out_dir && next_pos < argc)
        out_dir = argv[next_pos++];
    if (next_pos < argc)
        cmd_ref_file = argv[next_pos];

    if (num_threads <= 0)
    {
        num_threads = sysconf(_SC_NPROCESSORS_ONLN);
        if (num_threads < 1)
            num_threads = DEFAULT_THREADS;
    }

    if (!out_dir)
    {
        fprintf(stderr, "Error: no output directory specified (provide it as the "
                        "2nd argument or via -o/--output-dir)\n");
        return 1;
    }

    if (make_directory(out_dir) != 0)
    {
        fprintf(stderr, "Failed to create output directory: %s\n", out_dir);
        return 1;
    }

    log_time("Starting processing...\n");

    log_time("Loading chromosome mapping...\n");
    int valid_chr_count = 0;
    ChromMapEntry *chroms = NULL;
    char chrom_mapping_path[1024];
    if (chrom_mapping_file)
    {
        // Use the provided mapping file
        strncpy(chrom_mapping_path, chrom_mapping_file,
                sizeof(chrom_mapping_path) - 1);
        chrom_mapping_path[sizeof(chrom_mapping_path) - 1] = '\0';
    }
    else
    {
        // Use chrom_mapping.json in the same directory as the BAM file
        char *bam_dir = strdup(bam_file);
        char *last_slash = strrchr(bam_dir, '/');
        if (last_slash)
            *(last_slash + 1) = '\0'; // Keep the trailing slash
        else
            bam_dir[0] = '\0'; // No directory, use current directory

        snprintf(chrom_mapping_path, sizeof(chrom_mapping_path),
                 "%schrom_mapping.json", bam_dir);
        free(bam_dir);
    }

    int n_chroms = 0;
    char *mapping_ref_file = NULL;
    if (load_chrom_mapping(chrom_mapping_path, &chroms, &n_chroms,
                           &mapping_ref_file) != 0)
    {
        fprintf(stderr, "Error: Failed to load chromosome mapping from %s\n",
                chrom_mapping_path);
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
        fprintf(
            stderr,
            "No reference file specified in chrom_mapping.json or command line\n");
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

    sam_close(in);

    log_time("Allocating thread arguments...\n");
    ThreadArg *thread_args = malloc(header->n_targets * sizeof(ThreadArg));
    ChromosomeExport *chr_exports = calloc(header->n_targets, sizeof(ChromosomeExport));
    if (!thread_args || !chr_exports)
    {
        fprintf(stderr, "Failed to allocate thread arguments\n");
        free(thread_args);
        free(chr_exports);
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
        char *seq = faidx_fetch_seq(fai, entry->fasta, 0, header->target_len[tid],
                                    &seq_len);
        if (!seq || seq_len <= 0)
        {
            fprintf(stderr, "Failed to fetch sequence for %s\n", entry->fasta);
            if (seq)
                free(seq);
            continue;
        }
        // Set up ThreadArg as before, but use entry->name for output
        thread_args[valid_chr_count].bam_file = bam_file;
        thread_args[valid_chr_count].out_dir = out_dir;
        thread_args[valid_chr_count].tid = tid;
        thread_args[valid_chr_count].chr =
            strdup(entry->name); // Make a copy of the name
        thread_args[valid_chr_count].chr_len = header->target_len[tid];
        thread_args[valid_chr_count].min_mapq = min_mapq;
        thread_args[valid_chr_count].min_phred = min_phred;
        thread_args[valid_chr_count].min_cov = min_cov;
        thread_args[valid_chr_count].cap_cov = cap_cov;
        thread_args[valid_chr_count].keep_chg = keep_chg;
        thread_args[valid_chr_count].keep_chh = keep_chh;
        thread_args[valid_chr_count].compression = compression;
        thread_args[valid_chr_count].hdf5_chunk_size = hdf5_chunk_size;
        thread_args[valid_chr_count].chunk_size = chunk_size;
        thread_args[valid_chr_count].chr_seq = seq;
        thread_args[valid_chr_count].output_format = output_format;
        thread_args[valid_chr_count].split_context_files = split_context_files;
        thread_args[valid_chr_count].num_threads = num_threads;
        thread_args[valid_chr_count].read_level = read_level;
        thread_args[valid_chr_count].tile_size = tile_size;
        thread_args[valid_chr_count].chr_export = &chr_exports[valid_chr_count];
        valid_chr_count++;
    }
    free(chroms);

    // Chromosomes are processed one at a time; parallelism happens inside
    // process_chromosome(), which splits each chromosome into num_threads
    // regions. This bounds peak memory to a single chromosome's working set
    // while letting the user's --threads setting govern the actual work.
    log_time("Processing %d chromosome(s), %d worker thread(s) each...\n",
             valid_chr_count, num_threads);
    for (int i = 0; i < valid_chr_count; i++)
        process_chromosome(&thread_args[i]);

    SampleRunInfo run_info = {
        .bam_file = bam_file,
        .out_dir = out_dir,
        .reference = ref_file,
        .min_mapq = min_mapq,
        .min_phred = min_phred,
        .min_cov = min_cov,
        .cap_cov = cap_cov,
        .keep_chg = keep_chg,
        .keep_chh = keep_chh,
        .split_context_files = split_context_files,
        .read_level = read_level,
        .tile_size = tile_size,
    };
    if (write_extraction_manifest(&run_info, chr_exports, valid_chr_count) != 0)
        fprintf(stderr, "Warning: failed to write extraction manifest\n");

    log_time("Cleaning up HDF5...\n");
    cleanup_hdf5();

    log_time("Cleaning up thread arguments...\n");
    for (int i = 0; i < valid_chr_count; i++)
    {
        free(thread_args[i].chr_seq);
        free(thread_args[i].chr); // Free the copied chromosome name
    }
    free(thread_args);
    free(chr_exports);
    sam_hdr_destroy(header);
    fai_destroy(fai);
    if (mapping_ref_file)
        free(mapping_ref_file);

    log_time("Processing complete. MethylExtractor has finished.\n");

    return 0;
}
