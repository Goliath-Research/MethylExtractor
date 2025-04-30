#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <unistd.h>
#include <sys/stat.h>
#include <ctype.h>

#define MIN_MAPQ 10
#define MIN_PHRED 5
#define DEFAULT_FLAGS (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY)

// Structure to hold methylation counts
typedef struct {
    uint32_t meth;
    uint32_t unmeth;
} meth_counter;

// Structure to hold pileup data
typedef struct {
    int n; // Number of reads covering position
    bam_pileup1_t *pl; // Pileup data
} pileup_data;

// Global variables
static char *context_str[] = {"", "CpG", "CHG", "CHH"};

// Function to create directory
int make_directory(const char *path) {
    struct stat st = {0};
    if (stat(path, &st) == -1) {
        return mkdir(path, 0700);
    }
    return 0;
}

// Function to get sequence context
int get_context(faidx_t *fai, const char *chr, int pos, char *context, char *trinucl) {
    int len;
    char *seq = fai_fetch(fai, chr, pos-1, pos+1, &len);
    if (!seq || len < 3) {
        if (seq) free(seq);
        strcpy(context, "CHH");
        strcpy(trinucl, "NNN");
        return 3;
    }

    char c0 = toupper(seq[0]), c1 = toupper(seq[1]), c2 = toupper(seq[2]);
    free(seq);

    trinucl[0] = c0; trinucl[1] = c1; trinucl[2] = c2; trinucl[3] = '\0';

    if (c1 != 'C') {
        strcpy(context, "");
        return 0;
    }

    if (c2 == 'G') {
        strcpy(context, "CpG");
        return 1;
    } else if (c2 == 'A' || c2 == 'C' || c2 == 'T') {
        if (c0 == 'C' || c0 == 'A' || c0 == 'T') {
            strcpy(context, "CHH");
            return 3;
        } else {
            strcpy(context, "CHG");
            return 2;
        }
    } else {
        if (c0 == 'C' || c0 == 'A' || c0 == 'T' || c0 == 'N') {
            strcpy(context, "CHH");
            return 3;
        } else {
            strcpy(context, "CHG");
            return 2;
        }
    }
}

// Function to process CIGAR for read position
int cigar2readpos(const bam1_t *b, int refpos, int *readpos) {
    int i, rpos = 0, qpos = 0;
    uint32_t *cigar = bam_get_cigar(b);
    int k;

    for (k = 0; k < b->core.n_cigar; k++) {
        int op = bam_cigar_op(cigar[k]);
        int len = bam_cigar_oplen(cigar[k]);

        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            if (rpos + len > refpos) {
                *readpos = qpos + (refpos - rpos);
                return 0;
            }
            rpos += len;
            qpos += len;
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            rpos += len;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            qpos += len;
        }
    }
    return -1;
}

// Function to process pileup data
void process_pileup(pileup_data *p, meth_counter *c, int refpos, faidx_t *fai, const char *chr, FILE *out) {
    char context[4], trinucl[4];
    int ctx = get_context(fai, chr, refpos, context, trinucl);
    if (ctx == 0) return;

    c->meth = 0;
    c->unmeth = 0;

    for (int i = 0; i < p->n; i++) {
        bam_pileup1_t *pl = p->pl + i;
        const bam1_t *b = pl->b;

        if (b->core.flag & DEFAULT_FLAGS) continue;
        if (b->core.qual < MIN_MAPQ) continue;
        if (pl->is_del || pl->is_refskip) continue;

        int readpos;
        if (cigar2readpos(b, refpos, &readpos) != 0) continue;

        uint8_t *seq = bam_get_seq(b);
        uint8_t base = bam_seqi(seq, readpos);
        uint8_t *qual = bam_get_qual(b);
        if (qual[readpos] < MIN_PHRED) continue;

        char bbase = seq_nt16_str[base];
        if (bbase == 'C' || bbase == 'c') {
            c->meth++;
        } else if (bbase == 'T' || bbase == 't') {
            c->unmeth++;
        }
    }

    if (c->meth + c->unmeth > 0) {
        char strand = '+';
        fprintf(out, "%s\t%d\t%c\t%u\t%u\t%s\t%s\n",
                chr, refpos, strand, c->meth, c->unmeth, context, trinucl);
    }
}

// Main function
int main(int argc, char *argv[]) 
{
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <reference.fa> <input.bam> <output_dir>\n", argv[0]);
        return 1;
    }

    const char *ref_file = argv[1];
    const char *bam_file = argv[2];
    const char *out_dir = argv[3];

    // Create output directory
    if (make_directory(out_dir) != 0) {
        fprintf(stderr, "Failed to create output directory: %s\n", out_dir);
        return 1;
    }

    // Load reference
    faidx_t *fai = fai_load(ref_file);
    if (!fai) {
        fprintf(stderr, "Failed to load reference: %s\n", ref_file);
        return 1;
    }

    // Open BAM
    samFile *in = sam_open(bam_file, "r");
    if (!in) {
        fprintf(stderr, "Failed to open BAM: %s\n", bam_file);
        fai_destroy(fai);
        return 1;
    }

    bam_hdr_t *header = sam_hdr_read(in);
    if (!header) {
        fprintf(stderr, "Failed to read BAM header\n");
        sam_close(in);
        fai_destroy(fai);
        return 1;
    }

    hts_idx_t *idx = sam_index_load(in, bam_file);
    if (!idx) {
        fprintf(stderr, "BAM index file missing or corrupt\n");
        sam_hdr_destroy(header);
        sam_close(in);
        fai_destroy(fai);
        return 1;
    }

    // Process each chromosome
    for (int tid = 0; tid < header->n_targets; tid++) {
        const char *chr = header->target_name[tid];
        char out_path[1024];
        snprintf(out_path, sizeof(out_path), "%s/%s_cytosine_report.txt", out_dir, chr);
        FILE *out = fopen(out_path, "w");
        if (!out) {
            fprintf(stderr, "Failed to open output file: %s\n", out_path);
            continue;
        }

        hts_itr_t *iter = sam_itr_queryi(idx, tid, 0, header->target_len[tid]);
        if (!iter) {
            fprintf(stderr, "Failed to create iterator for %s\n", chr);
            fclose(out);
            continue;
        }

        bam_plp_t plp = bam_plp_init(NULL, NULL);
        bam_plp_set_maxcnt(plp, 8000);
        pileup_data p;

        bam1_t *b = bam_init1();
        while (sam_itr_next(in, iter, b) >= 0) {
            bam_plp_push(plp, b);
            while ((p.pl = bam_plp_next(plp, &tid, &p.n)) != NULL) {
                meth_counter c;
                process_pileup(&p, &c, p.pl[0].b->core.pos + 1, fai, chr, out);
            }
        }

        bam_plp_push(plp, NULL);
        while ((p.pl = bam_plp_next(plp, &tid, &p.n)) != NULL) {
            meth_counter c;
            process_pileup(&p, &c, p.pl[0].b->core.pos + 1, fai, chr, out);
        }

        bam_destroy1(b);
        bam_plp_destroy(plp);
        hts_itr_destroy(iter);
        fclose(out);
    }

    // Cleanup
    hts_idx_destroy(idx);
    sam_hdr_destroy(header);
    sam_close(in);
    fai_destroy(fai);
    return 0;
}