#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>
#include <pthread.h>

#define RUNOFFSET 99
#define BBM_VERSION 1

typedef struct
{
    int minMapq, minPhred, minDepth;
    FILE** output_fp;
    char* BAMName;
    htsFile* fp;
    hts_idx_t* bai;
    char* FastaName;
    int nThreads;
    unsigned long chunkSize;
} Config;

struct lastCall
{
    int32_t tid, pos;
    uint32_t nmethyl, nunmethyl;
};

const char* TriNucleotideContexts[25] = { "CAA", "CAC", "CAG", "CAT", "CAN",
                                         "CCA", "CCC", "CCG", "CCT", "CCN",
                                         "CGA", "CGC", "CGG", "CGT", "CGN",
                                         "CTA", "CTC", "CTG", "CTT", "CTN",
                                         "CNA", "CNC", "CNG", "CNT", "CNN" };

void writeCall(kstring_t* ks, Config* config, char* chrom, int32_t pos, uint32_t nmethyl, uint32_t nunmethyl, char base, char* context, const char* tnc)
{
    char str[10000];
    char strand = (base == 'C' || base == 'c') ? '+' : '-';
    snprintf(str, 10000, "%s\t%i\t%c\t%" PRIu32 "\t%" PRIu32 "\tC%s\t%s\n",
        chrom, pos + 1, strand, nmethyl, nunmethyl, context, tnc);
    kputs(str, ks);
}

char revcomp(char b)
{
    switch (b)
    {
    case 'A':
    case 'a':
        return 'T';
    case 'C':
    case 'c':
        return 'G';
    case 'G':
    case 'g':
        return 'C';
    case 'T':
    case 't':
        return 'A';
    default:
        return 'N';
    }
}

int getTriNucContext(char* seq, uint32_t offset, int seqlen, int direction)
{
    int rv = 0;
    char base;
    if ((direction > 0 && offset + 2 >= seqlen) || (direction < 0 && offset <= 1))
        rv = 4;
    else
    {
        base = *(seq + offset + 2 * direction);
        if (direction < 0)
            base = revcomp(base);
        switch (base)
        {
        case 'A':
        case 'a':
            rv = 0;
            break;
        case 'C':
        case 'c':
            rv = 1;
            break;
        case 'G':
        case 'g':
            rv = 2;
            break;
        case 'T':
        case 't':
            rv = 3;
            break;
        default:
            rv = 4;
            break;
        }
    }
    if ((direction > 0 && offset + 1 >= seqlen) || (direction < 0 && offset == 0))
        rv += 20;
    else
    {
        base = *(seq + offset + direction);
        if (direction < 0)
            base = revcomp(base);
        switch (base)
        {
        case 'A':
        case 'a':
            rv += 0;
            break;
        case 'C':
        case 'c':
            rv += 5;
            break;
        case 'G':
        case 'g':
            rv += 10;
            break;
        case 'T':
        case 't':
            rv += 15;
            break;
        default:
            rv += 20;
            break;
        }
    }
    return rv;
}

int isCpG(char* seq, uint32_t offset, int seqlen)
{
    if (offset + 1 >= seqlen)
        return 0;
    char b1 = seq[offset], b2 = seq[offset + 1];
    if ((b1 == 'C' || b1 == 'c') && (b2 == 'G' || b2 == 'g'))
        return 1;
    if (offset == 0)
        return 0;
    b1 = seq[offset - 1];
    b2 = seq[offset];
    if ((b1 == 'C' || b1 == 'c') && (b2 == 'G' || b2 == 'g'))
        return -1;
    return 0;
}

int isCHG(char* seq, uint32_t offset, int seqlen)
{
    if (offset + 2 >= seqlen || offset < 1)
        return 0;
    char b1 = seq[offset], b2 = seq[offset + 1], b3 = seq[offset + 2];
    char b0 = seq[offset - 1];
    if ((b1 == 'C' || b1 == 'c') && (b3 == 'G' || b3 == 'g') &&
        (b2 == 'A' || b2 == 'a' || b2 == 'C' || b2 == 'c' || b2 == 'T' || b2 == 't'))
        return 1;
    if ((b0 == 'C' || b0 == 'c') && (b2 == 'G' || b2 == 'g') &&
        (b1 == 'A' || b1 == 'a' || b1 == 'C' || b1 == 'c' || b1 == 'T' || b1 == 't'))
        return -1;
    return 0;
}

int isCHH(char* seq, uint32_t offset, int seqlen)
{
    if (offset + 2 >= seqlen || offset < 1)
        return 0;
    char b1 = seq[offset], b2 = seq[offset + 1], b3 = seq[offset + 2];
    char b0 = seq[offset - 1];
    if ((b1 == 'C' || b1 == 'c') &&
        (b2 == 'A' || b2 == 'a' || b2 == 'C' || b2 == 'c' || b2 == 'T' || b2 == 't') &&
        (b3 == 'A' || b3 == 'a' || b3 == 'C' || b3 == 'c' || b3 == 'T' || b3 == 't'))
        return 1;
    if ((b0 == 'C' || b0 == 'c') &&
        (b1 == 'A' || b1 == 'a' || b1 == 'C' || b1 == 'c' || b1 == 'T' || b1 == 't') &&
        (b2 == 'A' || b2 == 'a' || b2 == 'C' || b2 == 'c' || b2 == 'T' || b2 == 't'))
        return -1;
    return 0;
}

void writeBlank(kstring_t* ks, Config* config, char* chrom, int32_t pos, uint32_t localPos2, uint32_t* lastPos, char* seq, int seqlen)
{
    int triNucContext, direction;
    char context[3] = "HG";
    if (pos == -1)
        return;
    for (; *lastPos < pos; (*lastPos)++)
    {
        if ((direction = isCpG(seq, *lastPos - localPos2, seqlen)) != 0)
        {
            triNucContext = getTriNucContext(seq, *lastPos - localPos2, seqlen, direction);
            context[0] = 'G';
            context[1] = 0;
        }
        else if ((direction = isCHG(seq, *lastPos - localPos2, seqlen)) != 0)
        {
            triNucContext = getTriNucContext(seq, *lastPos - localPos2, seqlen, direction);
            context[0] = 'H';
            context[1] = 'G';
        }
        else if ((direction = isCHH(seq, *lastPos - localPos2, seqlen)) != 0)
        {
            triNucContext = getTriNucContext(seq, *lastPos - localPos2, seqlen, direction);
            context[0] = 'H';
            context[1] = 'H';
        }
        else
        {
            continue;
        }
        writeCall(ks, config, chrom, *lastPos, 0, 0, (direction > 0) ? 'C' : 'G', context, TriNucleotideContexts[triNucContext]);
    }
}

void* extractCalls(void* foo)
{
    Config* config = (Config*)foo;
    bam_hdr_t* hdr;
    bam_mplp_t iter;
    int ret, tid, i, seqlen, rv;
    hts_pos_t pos;
    int n_plp, strand, direction, tnc;
    uint32_t nmethyl, nunmethyl, localPos, localEnd, localTid, localPos2, lastPos;
    const bam_pileup1_t** plp;
    char* seq, base, context[3] = "HG";
    kstring_t* os;
    faidx_t* fai;
    hts_idx_t* bai;
    htsFile* fp;
    mplp_data* data;

    os = calloc(1, sizeof(kstring_t));

    if ((fai = fai_load(config->FastaName)) == NULL)
    {
        fprintf(stderr, "Couldn't open the index for %s!\n", config->FastaName);
        return NULL;
    }
    if ((fp = hts_open(config->BAMName, "rb")) == NULL)
    {
        fprintf(stderr, "Couldn't open %s for reading!\n", config->BAMName);
        return NULL;
    }
    if ((bai = sam_index_load(fp, config->BAMName)) == NULL)
    {
        fprintf(stderr, "Couldn't load the index for %s\n", config->BAMName);
        return NULL;
    }
    hdr = sam_hdr_read(fp);

    data = calloc(1, sizeof(mplp_data));
    data->config = config;
    data->hdr = hdr;
    data->fp = fp;

    plp = calloc(1, sizeof(bam_pileup1_t*));

    while (1)
    {
        pthread_mutex_lock(&positionMutex);
        localTid = globalTid;
        localPos = globalPos;
        localEnd = localPos + config->chunkSize;
        if (localTid >= hdr->n_targets)
        {
            pthread_mutex_unlock(&positionMutex);
            break;
        }
        if (globalEnd && localEnd > globalEnd)
            localEnd = globalEnd;
        globalPos = localEnd;
        if (globalEnd > 0 && globalPos >= globalEnd)
            globalTid = (uint32_t)-1;
        if (localTid < hdr->n_targets&& globalTid != (uint32_t)-1)
        {
            if (globalPos >= hdr->target_len[localTid])
            {
                localEnd = hdr->target_len[localTid];
                globalTid++;
                globalPos = 0;
            }
        }
        pthread_mutex_unlock(&positionMutex);

        localPos2 = localPos > 1 ? localPos - 2 : 0;
        lastPos = localPos;

        if (localTid >= hdr->n_targets || (globalEnd && localPos >= globalEnd))
            break;
        data->iter = sam_itr_queryi(bai, localTid, localPos, localEnd);

        seq = faidx_fetch_seq(fai, hdr->target_name[localTid], localPos2, localEnd + 10, &seqlen);
        if (seqlen < 0)
        {
            fprintf(stderr, "faidx_fetch_seq failed for %s:%" PRIu32 "-%" PRIu32 "\n",
                hdr->target_name[localTid], localPos2, localEnd);
            continue;
        }
        data->seq = seq;
        data->offset = localPos2;
        data->lseq = seqlen;

        data->ohash = initOlapHash();
        iter = bam_mplp_init(1, filter_func, (void**)&data);
        bam_mplp_set_maxcnt(iter, INT_MAX);
        bam_mplp_constructor(iter, custom_overlap_constructor);
        bam_mplp_destructor(iter, custom_overlap_destructor);

        while ((ret = bam_mplp64_auto(iter, &tid, &pos, &n_plp, plp)) > 0)
        {
            if (pos < localPos || pos >= localEnd)
                continue;

            if ((direction = isCpG(seq, pos - localPos2, seqlen)))
            {
                context[0] = 'G';
                context[1] = 0;
            }
            else if ((direction = isCHG(seq, pos - localPos2, seqlen)))
            {
                context[0] = 'H';
                context[1] = 'G';
            }
            else if ((direction = isCHH(seq, pos - localPos2, seqlen)))
            {
                context[0] = 'H';
                context[1] = 'H';
            }
            else
            {
                continue;
            }

            nmethyl = nunmethyl = 0;
            base = *(seq + pos - localPos2);
            for (i = 0; i < n_plp; i++)
            {
                if (plp[0][i].is_del || plp[0][i].is_refskip)
                    continue;
                strand = getStrand((plp[0] + i)->b);
                if (strand & 1)
                {
                    if (base != 'C' && base != 'c')
                        continue;
                }
                else
                {
                    if (base != 'G' && base != 'g')
                        continue;
                }
                rv = updateMetrics(config, plp[0] + i);
                if (rv > 0)
                    nmethyl++;
                else if (rv < 0)
                    nunmethyl++;
            }

            writeBlank(os, config, hdr->target_name[localTid], pos, localPos2, &lastPos, seq, seqlen);
            tnc = getTriNucContext(seq, pos - localPos2, seqlen, direction);
            writeCall(os, config, hdr->target_name[tid], pos, nmethyl, nunmethyl, base, context, TriNucleotideContexts[tnc]);
            lastPos = pos + 1;
        }
        bam_mplp_destroy(iter);
        writeBlank(os, config, hdr->target_name[localTid], localEnd, localPos2, &lastPos, seq, seqlen);
        hts_itr_destroy(data->iter);
        free(seq);
        destroyOlapHash(data->ohash);
    }

    pthread_mutex_lock(&outputMutex);
    if (os->l)
    {
        fputs(os->s, config->output_fp[0]);
        os->l = 0;
    }
    pthread_mutex_unlock(&outputMutex);

    free(os->s);
    free(os);
    bam_hdr_destroy(hdr);
    fai_destroy(fai);
    hts_close(fp);
    hts_idx_destroy(bai);
    free(data);
    free(plp);
    return NULL;
}

void extract_usage()
{
    fprintf(stderr, "\nUsage: MethylDackel extract [OPTIONS] <ref.fa> <sorted_alignments.bam>\n");
    fprintf(stderr,
        "\nOptions:\n"
        " -q INT           Minimum MAPQ threshold (default 10)\n"
        " -p INT           Minimum Phred threshold (default 5)\n"
        " -d INT           Minimum per-base depth (default 1)\n"
        " -@ INT           Number of threads (default 1)\n"
        " --chunkSize INT  Genome chunk size (default 1000000)\n"
        " -o STR           Output filename prefix\n"
        " -h               Show this help\n");
}

int extract_main(int argc, char* argv[])
{
    char* opref = NULL, * oname, * p;
    Config config;
    int c, i;

    globalTid = globalPos = globalEnd = 0;

    config.minMapq = 10;
    config.minPhred = 5;
    config.minDepth = 1;
    config.fp = NULL;
    config.bai = NULL;
    config.nThreads = 1;
    config.chunkSize = 1000000;

    static struct option lopts[] = {
        {"opref", 1, NULL, 'o'},
        {"minDepth", 1, NULL, 'd'},
        {"chunkSize", 1, NULL, 19},
        {"help", 0, NULL, 'h'},
        {0, 0, NULL, 0} };
    while ((c = getopt_long(argc, argv, "hq:p:o:d:@:", lopts, NULL)) >= 0)
    {
        switch (c)
        {
        case 'h':
            extract_usage();
            return 0;
        case 'o':
            opref = strdup(optarg);
            break;
        case 'd':
            config.minDepth = atoi(optarg);
            if (config.minDepth < 1)
            {
                fprintf(stderr, "Minimum depth must be >= 1!\n");
                return 1;
            }
            break;
        case 19:
            config.chunkSize = strtoul(optarg, NULL, 10);
            if (config.chunkSize < 1)
            {
                fprintf(stderr, "Chunk size must be >= 1!\n");
                return 1;
            }
            break;
        case 'q':
            config.minMapq = atoi(optarg);
            break;
        case 'p':
            config.minPhred = atoi(optarg);
            break;
        case '@':
            config.nThreads = atoi(optarg);
            break;
        default:
            fprintf(stderr, "Invalid option '%c'\n", c);
            extract_usage();
            return 1;
        }
    }

    if (argc - optind < 2)
    {
        fprintf(stderr, "Must supply reference genome and BAM file!\n");
        extract_usage();
        return -1;
    }

    if (config.minPhred < 1)
    {
        fprintf(stderr, "-p %i invalid, resetting to 1.\n", config.minPhred);
        config.minPhred = 1;
    }
    if (config.minMapq < 0)
    {
        fprintf(stderr, "-q %i invalid, resetting to 0.\n", config.minMapq);
        config.minMapq = 0;
    }

    config.FastaName = argv[optind];
    config.BAMName = argv[optind + 1];

    if ((config.fp = hts_open(argv[optind + 1], "rb")) == NULL)
    {
        fprintf(stderr, "Couldn't open %s!\n", argv[optind + 1]);
        return -4;
    }
    if ((config.bai = sam_index_load(config.fp, argv[optind + 1])) == NULL)
    {
        fprintf(stderr, "Couldn't load index for %s, attempting to build.\n", argv[optind + 1]);
        if (bam_index_build(argv[optind + 1], 0) < 0 || (config.bai = sam_index_load(config.fp, argv[optind + 1])) == NULL)
        {
            fprintf(stderr, "Couldn't build/load index for %s!\n", argv[optind + 1]);
            return -5;
        }
    }

    config.output_fp = malloc(sizeof(FILE*) * 3);
    if (opref == NULL)
    {
        opref = strdup(argv[optind + 1]);
        p = strrchr(opref, '.');
        if (p != NULL)
            *p = '\0';
    }
    oname = malloc(strlen(opref) + 21);
    sprintf(oname, "%s.cytosine_report.txt", opref);
    config.output_fp[0] = fopen(oname, "w");
    config.output_fp[1] = config.output_fp[0];
    config.output_fp[2] = config.output_fp[0];
    if (config.output_fp[0] == NULL)
    {
        fprintf(stderr, "Couldn't open output file!\n");
        return -3;
    }

    pthread_mutex_init(&positionMutex, NULL);
    pthread_mutex_init(&outputMutex, NULL);
    pthread_t* threads = calloc(config.nThreads, sizeof(pthread_t));
    for (i = 0; i < config.nThreads; i++)
        pthread_create(threads + i, NULL, &extractCalls, &config);
    for (i = 0; i < config.nThreads; i++)
        pthread_join(threads[i], NULL);
    free(threads);
    pthread_mutex_destroy(&outputMutex);
    pthread_mutex_destroy(&positionMutex);

    hts_close(config.fp);
    fclose(config.output_fp[0]);
    hts_idx_destroy(config.bai);
    free(opref);
    free(oname);
    free(config.output_fp);

    return 0;
}