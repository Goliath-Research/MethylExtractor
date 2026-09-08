#include "methyl_extractor.h"

int getRealStrand(bam1_t *b)
{
    char *XG = (char *)bam_aux_get(b, "XG");
    if (XG && (XG[1] == 'C' || XG[1] == 'G'))
    {
        if (XG[1] == 'C')
        {
            if (b->core.flag & BAM_FREVERSE)
                return (b->core.flag & BAM_FREAD1) ? 1 : 3;
            return (b->core.flag & BAM_FREAD1) ? 3 : 1;
        }
        else
        {
            if (b->core.flag & BAM_FREVERSE)
                return (b->core.flag & BAM_FREAD1) ? 4 : 2;
            return (b->core.flag & BAM_FREAD1) ? 2 : 4;
        }
    }
    if (b->core.flag & BAM_FPAIRED)
    {
        if ((b->core.flag & (BAM_FREAD1 | BAM_FREVERSE)) ==
            (BAM_FREAD1 | BAM_FREVERSE))
            return 2;
        else if (b->core.flag & BAM_FREAD1)
            return 1;
        else if ((b->core.flag & (BAM_FREAD2 | BAM_FREVERSE)) ==
                 (BAM_FREAD2 | BAM_FREVERSE))
            return 1;
        else
            return 2;
    }
    return (b->core.flag & BAM_FREVERSE) ? 2 : 1;
}

void *process_region_direct(void *arg)
{
    RegionArg *r = (RegionArg *)arg;
    ThreadArg *t = &r->base;

    samFile *in = sam_open(t->bam_file, "r");
    if (!in)
        return NULL;
    /* One handle per region thread (htslib is not multi-reader on one fp).
     * Index is shared read-only. Header is loaded once in main and not
     * re-read here — iterators seek via the BAM index. */
    if (t->bgzf_threads > 1)
        hts_set_threads(in, t->bgzf_threads);
    hts_itr_t *iter = sam_itr_queryi(r->idx, t->tid, t->start_pos, t->end_pos);
    if (!iter)
    {
        sam_close(in);
        return NULL;
    }

    bam1_t *b = bam_init1();

    while (sam_itr_next(in, iter, b) >= 0)
    {
        bam1_core_t *c = &b->core;
        int used_this_read = 0;

        r->filter_stats.reads_seen++;

        if (c->flag & BAM_FUNMAP)
        {
            r->filter_stats.reads_dropped_unmapped++;
            continue;
        }
        if (c->flag & BAM_FSECONDARY)
        {
            r->filter_stats.reads_dropped_secondary++;
            continue;
        }
        if (c->flag & BAM_FQCFAIL)
        {
            r->filter_stats.reads_dropped_qc_fail++;
            continue;
        }
        if (c->flag & BAM_FDUP)
        {
            r->filter_stats.reads_dropped_duplicate++;
            continue;
        }
        if (c->flag & BAM_FSUPPLEMENTARY)
        {
            r->filter_stats.reads_dropped_supplementary++;
            continue;
        }
        if (c->qual < t->min_mapq)
        {
            r->filter_stats.reads_dropped_low_mapq++;
            continue;
        }
        uint8_t *nh = bam_aux_get(b, "NH");
        if (nh && bam_aux2i(nh) > 1)
        {
            r->filter_stats.reads_dropped_multimap++;
            continue;
        }

        int strand = getRealStrand(b);
        if (strand == 0)
        {
            r->filter_stats.reads_dropped_no_strand++;
            continue;
        }

        // Overlapping mate de-duplication (coordinate-based clip). For an FR
        // read pair the two mates overlap; to count each reference position
        // exactly once, the left mate (smaller start) yields the overlapping
        // region to the right mate by skipping reference positions >= the
        // mate's start. On equal starts (dovetailed pairs) READ2 yields to
        // READ1. The decision uses only this read's own fields, so it stays
        // consistent even when the two mates fall in different region threads.
        int32_t clip_from = INT32_MAX;
        if ((c->flag & BAM_FPAIRED) && !(c->flag & BAM_FMUNMAP) &&
            c->mtid == c->tid)
        {
            if (c->pos < c->mpos)
                clip_from = c->mpos;
            else if (c->pos == c->mpos && (c->flag & BAM_FREAD2))
                clip_from = c->mpos;
        }

        uint8_t *seq = bam_get_seq(b);
        uint8_t *qual = bam_get_qual(b);
        int32_t pos = c->pos;
        uint32_t *cigar = bam_get_cigar(b);
        int qpos = 0;
        uint8_t *xm_aux = bam_aux_get(b, "XM");
        const char *xmz = NULL;
        if (xm_aux && xm_aux[0] == 'Z')
            xmz = (const char *)(xm_aux + 1);

        int32_t mhap_pos[MHAP_MAX_CPG_PER_READ];
        uint8_t mhap_meth[MHAP_MAX_CPG_PER_READ];
        int mhap_n = 0;

        for (int k = 0; k < c->n_cigar; ++k)
        {
            int op = bam_cigar_op(cigar[k]);
            int len = bam_cigar_oplen(cigar[k]);

            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
            {
                for (int i = 0; i < len; ++i)
                {
                    if (pos >= clip_from || pos < t->start_pos ||
                        pos >= t->end_pos)
                    {
                        if (pos >= clip_from && pos < t->end_pos &&
                            pos >= t->start_pos)
                            r->filter_stats.bases_skipped_overlap_clip++;
                        pos++;
                        qpos++;
                        continue;
                    }
                    char ref = toupper(t->chr_seq[pos]);
                    if (ref != 'C' && ref != 'G')
                    {
                        r->filter_stats.bases_skipped_non_cytosine++;
                        pos++;
                        qpos++;
                        continue;
                    }
                    if (qual[qpos] < t->min_phred)
                    {
                        r->filter_stats.bases_skipped_low_phred++;
                        pos++;
                        qpos++;
                        continue;
                    }

                    int base = bam_seqi(seq, qpos);
                    uint32_t rel = (uint32_t)(pos - t->start_pos);
                    int called = -1;
                    if (xmz)
                    {
                        char xc = xmz[qpos];
                        if (xc == 'Z' || xc == 'X' || xc == 'H' || xc == 'U')
                            called = 1;
                        else if (xc == 'z' || xc == 'x' || xc == 'h' || xc == 'u')
                            called = 0;
                    }
                    if (called < 0)
                    {
                        if (ref == 'C' && (strand == 1 || strand == 3))
                        {
                            if (base == 2)
                                called = 1;
                            else if (base == 8)
                                called = 0;
                        }
                        else if (ref == 'G' && (strand == 2 || strand == 4))
                        {
                            if (base == 4)
                                called = 1;
                            else if (base == 1)
                                called = 0;
                        }
                    }

                    if (called == 1)
                    {
                        r->counts.mC[rel]++;
                        r->filter_stats.bases_counted_methylated++;
                        used_this_read = 1;
                    }
                    else if (called == 0)
                    {
                        r->counts.uC[rel]++;
                        r->filter_stats.bases_counted_unmethylated++;
                        used_this_read = 1;
                    }

                    if (t->mhap && r->mhap && called >= 0 &&
                        mhap_n < MHAP_MAX_CPG_PER_READ)
                    {
                        int cpg = isCpG(t->chr_seq, pos, (int)t->chr_len);
                        if (cpg != 0)
                        {
                            int32_t canonical = (cpg == -1) ? (int32_t)(pos - 1)
                                                            : (int32_t)pos;
                            if (mhap_n == 0 ||
                                mhap_pos[mhap_n - 1] != canonical)
                            {
                                mhap_pos[mhap_n] = canonical;
                                mhap_meth[mhap_n] = (uint8_t)called;
                                mhap_n++;
                            }
                        }
                    }
                    pos++;
                    qpos++;
                }
            }
            else if (op == BAM_CDEL || op == BAM_CREF_SKIP)
                pos += len;
            else if (op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)
                qpos += len;
        }

        if (t->mhap && r->mhap && mhap_n > 0)
            mhap_store_add(r->mhap, c->pos, (int8_t)strand, mhap_pos, mhap_meth,
                           mhap_n);

        if (used_this_read)
            r->filter_stats.reads_used++;

        if (r->rl && r->rl->n_active > 0 && r->rl_tile_start && r->rl_tile_end)
        {
            for (int ci = 0; ci < r->rl->n_active; ++ci)
            {
                if (r->rl_tile_end[ci] <= r->rl_tile_start[ci])
                    continue;
                ReadLevelContext *ctx = &r->rl->ctx[ci];
                read_level_accumulate(
                    &ctx->tiles, ctx->hist, r->rl_tile_start[ci],
                    r->rl_tile_end[ci], t->chr_seq, seq, qual, cigar, c->n_cigar,
                    c->pos, clip_from, strand, t->min_phred, t->chr_len);
            }
        }
    }

    bam_destroy1(b);
    hts_itr_destroy(iter);
    sam_close(in);
    return NULL;
}

static void record_context_export(ChromosomeExport *export_out, const char *context,
                                  const char *output_path, const MethylStats *stats)
{
    if (!export_out || export_out->n_contexts >= MAX_CONTEXTS_PER_CHR)
        return;

    ContextExport *ctx = &export_out->contexts[export_out->n_contexts++];
    strncpy(ctx->context, context, sizeof(ctx->context) - 1);
    ctx->context[sizeof(ctx->context) - 1] = '\0';
    strncpy(ctx->output_path, output_path, sizeof(ctx->output_path) - 1);
    ctx->output_path[sizeof(ctx->output_path) - 1] = '\0';
    if (stats)
        ctx->stats = *stats;
}

void process_chromosome(ThreadArg *targ)
{
    uint64_t chrom_t0 = now_ms();
    ChromosomeTiming *timing = targ->timing;
    if (timing)
    {
        memset(timing, 0, sizeof(*timing));
        strncpy(timing->chromosome, targ->chr, sizeof(timing->chromosome) - 1);
        timing->chromosome[sizeof(timing->chromosome) - 1] = '\0';
        timing->rss_est_bytes = estimate_chromosome_rss_bytes(
            targ->chr_len, targ->keep_chg, targ->keep_chh, targ->read_level,
            targ->split_context_files);
    }

    ChromosomeExport *export_out = targ->chr_export;
    if (export_out)
    {
        memset(export_out, 0, sizeof(*export_out));
        strncpy(export_out->chromosome, targ->chr, sizeof(export_out->chromosome) - 1);
        export_out->chromosome[sizeof(export_out->chromosome) - 1] = '\0';
    }

    uint64_t t_enum = now_ms();
    size_t site_count = count_methylation_sites(targ->chr_seq, targ->chr_len,
                                                targ->keep_chg, targ->keep_chh);
    if (site_count == 0)
    {
        if (timing)
            timing->total_ms = elapsed_ms(chrom_t0);
        return;
    }

    MethylRecord *buffer = calloc(site_count, sizeof(MethylRecord));
    initialize_buffer(buffer, site_count, targ->chr_seq, targ->chr_len,
                      targ->keep_chg, targ->keep_chh);
    if (timing)
        timing->site_enum_ms = elapsed_ms(t_enum);

    uint32_t *site_positions = malloc(site_count * sizeof(uint32_t));
    for (size_t i = 0; i < site_count; i++)
        site_positions[i] = buffer[i].pos - 1;

    int n_threads = targ->num_threads;
    if (n_threads < 1)
        n_threads = 1;
    if ((size_t)n_threads > site_count)
        n_threads = (int)site_count;

    // Load the BAM index once and share it (read-only) across region threads.
    samFile *idx_fp = sam_open(targ->bam_file, "r");
    hts_idx_t *idx = idx_fp ? sam_index_load(idx_fp, targ->bam_file) : NULL;
    if (!idx)
    {
        fprintf(stderr,
                "Failed to load BAM index for %s (is the BAM sorted and "
                "indexed?)\n",
                targ->chr);
        if (idx_fp)
            sam_close(idx_fp);
        free(site_positions);
        free(buffer);
        return;
    }

    size_t sites_per_thread = (site_count + n_threads - 1) / n_threads;

    ReadLevelChromData rl_data;
    memset(&rl_data, 0, sizeof(rl_data));
    if (targ->read_level)
    {
        int ctx_ids[] = {CONTEXT_CPG, CONTEXT_CHG, CONTEXT_CHH};
        int ctx_keep[] = {1, targ->keep_chg, targ->keep_chh};
        for (int ci = 0; ci < 3; ++ci)
        {
            if (!ctx_keep[ci])
                continue;
            if (!targ->split_context_files && ctx_ids[ci] != CONTEXT_CPG)
                continue;

            ReadLevelContext *rc = &rl_data.ctx[rl_data.n_active];
            if (build_read_level_tiles(targ->chr_seq, targ->chr_len, ctx_ids[ci],
                                       targ->tile_size, &rc->tiles) != 0)
            {
                fprintf(stderr, "Failed to build read-level tiles for %s %s\n",
                        targ->chr, get_context_string(ctx_ids[ci]));
                hts_idx_destroy(idx);
                if (idx_fp)
                    sam_close(idx_fp);
                free(site_positions);
                free(buffer);
                free_read_level_chrom_data(&rl_data);
                return;
            }
            if (rc->tiles.n_tiles == 0)
            {
                free_read_level_tiles(&rc->tiles);
                continue;
            }

            size_t hist_size =
                rc->tiles.n_tiles * (size_t)(1U << (unsigned)targ->tile_size);
            rc->hist = calloc(hist_size, sizeof(uint32_t));
            if (!rc->hist)
            {
                free_read_level_tiles(&rc->tiles);
                hts_idx_destroy(idx);
                if (idx_fp)
                    sam_close(idx_fp);
                free(site_positions);
                free(buffer);
                free_read_level_chrom_data(&rl_data);
                return;
            }
            rc->active = 1;
            rl_data.n_active++;
        }
    }

    RegionArg *regions = calloc(n_threads, sizeof(RegionArg));
    pthread_t *threads = malloc(n_threads * sizeof(pthread_t));

    // Only create threads for non-empty site ranges; track how many were
    // actually launched so join/merge never touch uninitialized entries.
    int n_created = 0;
    for (int i = 0; i < n_threads; i++)
    {
        size_t start_idx = (size_t)i * sites_per_thread;
        if (start_idx >= site_count)
            break;

        RegionArg *r = &regions[i];
        memset(&r->filter_stats, 0, sizeof(r->filter_stats));
        r->base = *targ;
        r->idx = idx;
        r->start_site_idx = start_idx;
        r->end_site_idx = start_idx + sites_per_thread;
        if (r->end_site_idx > site_count)
            r->end_site_idx = site_count;

        r->base.start_pos = site_positions[r->start_site_idx];
        r->base.end_pos = (r->end_site_idx < site_count)
                              ? site_positions[r->end_site_idx - 1] + 10
                              : targ->chr_len;

        uint32_t range = r->base.end_pos - r->base.start_pos;
        r->counts.mC = calloc(range, sizeof(uint32_t));
        r->counts.uC = calloc(range, sizeof(uint32_t));
        r->counts.size = range;

        if (rl_data.n_active > 0)
        {
            r->rl = &rl_data;
            r->rl_tile_start = calloc(rl_data.n_active, sizeof(size_t));
            r->rl_tile_end = calloc(rl_data.n_active, sizeof(size_t));
            if (r->rl_tile_start && r->rl_tile_end)
            {
                for (int ci = 0; ci < rl_data.n_active; ++ci)
                    read_level_tile_range_for_region(
                        &rl_data.ctx[ci].tiles, r->base.start_pos,
                        r->base.end_pos, &r->rl_tile_start[ci],
                        &r->rl_tile_end[ci]);
            }
        }
        r->mhap = NULL;
        if (targ->mhap)
        {
            r->mhap = malloc(sizeof(MhapStore));
            if (r->mhap)
                mhap_store_init(r->mhap);
        }

        pthread_create(&threads[i], NULL, process_region_direct, r);
        n_created++;
    }

    uint64_t t_bam = now_ms();
    for (int i = 0; i < n_created; i++)
        pthread_join(threads[i], NULL);
    if (timing)
        timing->bam_scan_ms = elapsed_ms(t_bam);

    FilterStats chr_filters;
    memset(&chr_filters, 0, sizeof(chr_filters));
    for (int i = 0; i < n_created; i++)
        merge_filter_stats(&chr_filters, &regions[i].filter_stats);
    if (export_out)
        export_out->filter_stats = chr_filters;

    uint64_t t_merge = now_ms();
    MhapStore mhap_merged;
    mhap_store_init(&mhap_merged);
    for (int i = 0; i < n_created; i++)
    {
        uint32_t offset = regions[i].base.start_pos;
        for (size_t j = regions[i].start_site_idx; j < regions[i].end_site_idx;
             j++)
        {
            uint32_t gpos = site_positions[j];
            uint32_t rel = gpos - offset;
            buffer[j].mC += regions[i].counts.mC[rel];
            buffer[j].uC += regions[i].counts.uC[rel];
        }
        free(regions[i].counts.mC);
        free(regions[i].counts.uC);
        free(regions[i].rl_tile_start);
        free(regions[i].rl_tile_end);
        if (regions[i].mhap)
        {
            mhap_store_merge(&mhap_merged, regions[i].mhap);
            mhap_store_free(regions[i].mhap);
            free(regions[i].mhap);
        }
    }
    if (timing)
        timing->merge_ms = elapsed_ms(t_merge);

    uint64_t t_write = now_ms();

    if (targ->mhap && mhap_merged.n_reads > 0)
    {
        char mhap_path[1024];
        snprintf(mhap_path, sizeof(mhap_path), "%s/%s-CG.mhap.h5", targ->out_dir,
                 targ->chr);
        if (write_mhap_h5(mhap_path, targ->chr, "CG", targ->compression,
                          &mhap_merged) != 0)
        {
            fprintf(stderr, "Warning: failed to write haplotype sidecar %s\n",
                    mhap_path);
        }
    }
    mhap_store_free(&mhap_merged);

    if (rl_data.n_active > 0)
    {
        for (int ci = 0; ci < rl_data.n_active; ++ci)
        {
            ReadLevelContext *rc = &rl_data.ctx[ci];
            char pattern_path[1024];
            snprintf(pattern_path, sizeof(pattern_path), "%s/%s-%s.patterns.h5",
                     targ->out_dir, targ->chr,
                     get_context_string(rc->tiles.context));
            if (write_read_level_patterns_h5(
                    pattern_path, get_context_string(rc->tiles.context),
                    targ->tile_size, 1, targ->compression, &rc->tiles,
                    rc->hist) != 0)
            {
                fprintf(stderr, "Warning: failed to write read-level sidecar %s\n",
                        pattern_path);
            }
        }
    }

    if (targ->split_context_files)
    {
        for (int ctx = CONTEXT_CPG; ctx <= CONTEXT_CHH; ++ctx)
        {
            if ((ctx == CONTEXT_CPG) || (ctx == CONTEXT_CHG && targ->keep_chg) ||
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
                const char *ext = (targ->output_format == OUTPUT_TXT) ? ".txt" : ".h5";
                snprintf(out_path, sizeof(out_path), "%s/%s-%s%s", targ->out_dir,
                         targ->chr, get_context_string(ctx), ext);

                ExtractionMeta meta = {
                    .chromosome = targ->chr,
                    .context = get_context_string(ctx),
                    .output_path = out_path,
                    .min_mapq = targ->min_mapq,
                    .min_phred = targ->min_phred,
                    .min_cov = targ->min_cov,
                    .cap_cov = targ->cap_cov,
                    .timing = timing,
                };
                MethylStats ctx_stats;
                flush_buffer(out_path, ctx_buffer, n_ctx_records, targ->compression,
                             targ->hdf5_chunk_size, 0, targ->min_cov, targ->cap_cov,
                             targ->output_format, &meta, &chr_filters, &ctx_stats);
                record_context_export(export_out, get_context_string(ctx), out_path,
                                      &ctx_stats);

                free(ctx_buffer);

                // Note: External compression removed to maintain HDF5 compatibility
            }
        }
    }
    else
    {
        // Output file name - use appropriate extension based on output format
        char out_path[1024];
        const char *ext = (targ->output_format == OUTPUT_TXT) ? ".txt"
                          : ".h5"; // For OUTPUT_BOTH, use .h5
        snprintf(out_path, sizeof(out_path), "%s/%s%s", targ->out_dir, targ->chr,
                 ext);
        ExtractionMeta meta = {
            .chromosome = targ->chr,
            .context = "ALL",
            .output_path = out_path,
            .min_mapq = targ->min_mapq,
            .min_phred = targ->min_phred,
            .min_cov = targ->min_cov,
            .cap_cov = targ->cap_cov,
            .timing = timing,
        };
        MethylStats chr_stats;
        flush_buffer(out_path, buffer, site_count, targ->compression,
                     targ->hdf5_chunk_size, 0, targ->min_cov, targ->cap_cov,
                     targ->output_format, &meta, &chr_filters, &chr_stats);
        record_context_export(export_out, "ALL", out_path, &chr_stats);

        // Note: External compression removed to maintain HDF5 compatibility
    }

    if (timing)
        timing->hdf5_write_ms = elapsed_ms(t_write);

    hts_idx_destroy(idx);
    if (idx_fp)
        sam_close(idx_fp);
    free(site_positions);
    free(buffer);
    free(threads);
    free(regions);
    free_read_level_chrom_data(&rl_data);

    if (timing)
    {
        timing->total_ms = elapsed_ms(chrom_t0);
        log_time("%s site_enum=%lums bam_scan=%lums merge=%lums write=%lums "
                 "qc=%lums total=%lums\n",
                 targ->chr, (unsigned long)timing->site_enum_ms,
                 (unsigned long)timing->bam_scan_ms,
                 (unsigned long)timing->merge_ms,
                 (unsigned long)timing->hdf5_write_ms,
                 (unsigned long)timing->qc_json_ms,
                 (unsigned long)timing->total_ms);
    }
}

typedef struct
{
    ThreadArg *jobs;
    uint64_t *rss_est;
    int n;
    int next;
    int in_flight;
    int max_inflight;
    uint64_t in_flight_rss;
    uint64_t max_rss;
    pthread_mutex_t mu;
    pthread_cond_t cv;
} ChromPool;

static void *chrom_pool_worker(void *arg)
{
    ChromPool *pool = (ChromPool *)arg;
    for (;;)
    {
        ThreadArg *job = NULL;
        uint64_t est = 0;
        pthread_mutex_lock(&pool->mu);
        for (;;)
        {
            if (pool->next >= pool->n)
            {
                pthread_mutex_unlock(&pool->mu);
                return NULL;
            }
            est = pool->rss_est[pool->next];
            int at_cap = pool->in_flight >= pool->max_inflight;
            int would_exceed =
                (pool->in_flight > 0 && pool->in_flight_rss + est > pool->max_rss);
            if (at_cap || would_exceed)
            {
                pthread_cond_wait(&pool->cv, &pool->mu);
                continue;
            }
            job = &pool->jobs[pool->next++];
            pool->in_flight++;
            pool->in_flight_rss += est;
            break;
        }
        pthread_mutex_unlock(&pool->mu);

        process_chromosome(job);

        pthread_mutex_lock(&pool->mu);
        pool->in_flight--;
        pool->in_flight_rss -= est;
        pthread_cond_broadcast(&pool->cv);
        pthread_mutex_unlock(&pool->mu);
    }
}

void process_chromosomes_parallel(ThreadArg *args, int n_chroms, int chrom_parallel,
                                  uint64_t max_rss_bytes)
{
    if (n_chroms <= 0)
        return;
    if (chrom_parallel < 1)
        chrom_parallel = 1;
    if (chrom_parallel > n_chroms)
        chrom_parallel = n_chroms;

    if (chrom_parallel == 1)
    {
        for (int i = 0; i < n_chroms; i++)
            process_chromosome(&args[i]);
        return;
    }

    uint64_t *rss_est = malloc((size_t)n_chroms * sizeof(uint64_t));
    if (!rss_est)
    {
        for (int i = 0; i < n_chroms; i++)
            process_chromosome(&args[i]);
        return;
    }
    for (int i = 0; i < n_chroms; i++)
        rss_est[i] = estimate_chromosome_rss_bytes(
            args[i].chr_len, args[i].keep_chg, args[i].keep_chh,
            args[i].read_level, args[i].split_context_files);

    ChromPool pool;
    memset(&pool, 0, sizeof(pool));
    pool.jobs = args;
    pool.rss_est = rss_est;
    pool.n = n_chroms;
    pool.max_inflight = chrom_parallel;
    pool.max_rss = max_rss_bytes ? max_rss_bytes : UINT64_MAX;
    pthread_mutex_init(&pool.mu, NULL);
    pthread_cond_init(&pool.cv, NULL);

    pthread_t *workers = malloc((size_t)chrom_parallel * sizeof(pthread_t));
    int n_workers = 0;
    for (int i = 0; i < chrom_parallel; i++)
    {
        if (pthread_create(&workers[i], NULL, chrom_pool_worker, &pool) == 0)
            n_workers++;
        else
            break;
    }
    if (n_workers == 0)
    {
        for (int i = 0; i < n_chroms; i++)
            process_chromosome(&args[i]);
    }
    else
    {
        for (int i = 0; i < n_workers; i++)
            pthread_join(workers[i], NULL);
    }

    pthread_mutex_destroy(&pool.mu);
    pthread_cond_destroy(&pool.cv);
    free(workers);
    free(rss_est);
}

int load_chrom_mapping(const char *filename, ChromMapEntry **entries,
                       int *n_entries, char **reference_file)
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
            fprintf(
                stderr,
                "Warning: Skipping chromosome entry with missing required fields\n");
            continue;
        }

        e->extract = 1;
        (*n_entries)++;
    }
    cJSON_Delete(json);
    return 0;
}
