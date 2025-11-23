#include "methyl_extractor.h"

int getRealStrand(bam1_t *b) {
  char *XG = (char *)bam_aux_get(b, "XG");
  if (XG && (XG[1] == 'C' || XG[1] == 'G')) {
    if (XG[1] == 'C') {
      if (b->core.flag & BAM_FREVERSE)
        return (b->core.flag & BAM_FREAD1) ? 1 : 3;
      return (b->core.flag & BAM_FREAD1) ? 3 : 1;
    } else {
      if (b->core.flag & BAM_FREVERSE)
        return (b->core.flag & BAM_FREAD1) ? 4 : 2;
      return (b->core.flag & BAM_FREAD1) ? 2 : 4;
    }
  }
  if (b->core.flag & BAM_FPAIRED) {
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

void *process_region_direct(void *arg) {
  RegionArg *r = (RegionArg *)arg;
  ThreadArg *t = &r->base;

  samFile *in = sam_open(t->bam_file, "r");
  if (!in)
    return NULL;
  hts_idx_t *idx = sam_index_load(in, t->bam_file);
  if (!idx) {
    sam_close(in);
    return NULL;
  }
  hts_itr_t *iter = sam_itr_queryi(idx, t->tid, t->start_pos, t->end_pos);
  if (!iter) {
    hts_idx_destroy(idx);
    sam_close(in);
    return NULL;
  }

  bam_hdr_t *hdr = sam_hdr_read(in);
  bam1_t *b = bam_init1();

  while (sam_itr_next(in, iter, b) >= 0) {
    bam1_core_t *c = &b->core;
    if (c->flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP |
                   BAM_FSUPPLEMENTARY))
      continue;
    if (c->qual < t->min_mapq)
      continue;
    uint8_t *nh = bam_aux_get(b, "NH");
    if (nh && bam_aux2i(nh) > 1)
      continue;

    int strand = getRealStrand(b);
    if (strand == 0)
      continue;

    uint8_t *seq = bam_get_seq(b);
    uint8_t *qual = bam_get_qual(b);
    int32_t pos = c->pos;
    uint32_t *cigar = bam_get_cigar(b);
    int qpos = 0;

    for (int k = 0; k < c->n_cigar; ++k) {
      int op = bam_cigar_op(cigar[k]);
      int len = bam_cigar_oplen(cigar[k]);

      if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
        for (int i = 0; i < len; ++i) {
          if (pos < t->start_pos || pos >= t->end_pos) {
            pos++;
            qpos++;
            continue;
          }
          char ref = toupper(t->chr_seq[pos]);
          if (ref != 'C' && ref != 'G') {
            pos++;
            qpos++;
            continue;
          }
          if (qual[qpos] < t->min_phred) {
            pos++;
            qpos++;
            continue;
          }

          int base = bam_seqi(seq, qpos);
          uint32_t rel = (uint32_t)(pos - t->start_pos);

          if (ref == 'C' && (strand == 1 || strand == 3)) {
            if (base == 2)
              r->counts.mC[rel]++; // C->C methylated
            else if (base == 8)
              r->counts.uC[rel]++; // C->T unmethylated
          } else if (ref == 'G' && (strand == 2 || strand == 4)) {
            if (base == 4)
              r->counts.mC[rel]++; // G->G methylated
            else if (base == 1)
              r->counts.uC[rel]++; // G->A unmethylated
          }
          pos++;
          qpos++;
        }
      } else if (op == BAM_CDEL || op == BAM_CREF_SKIP)
        pos += len;
      else if (op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)
        qpos += len;
    }
  }

  bam_destroy1(b);
  sam_hdr_destroy(hdr);
  hts_itr_destroy(iter);
  hts_idx_destroy(idx);
  sam_close(in);
  return NULL;
}

void process_chromosome(ThreadArg *targ) {
  size_t site_count = count_methylation_sites(targ->chr_seq, targ->chr_len,
                                              targ->keep_chg, targ->keep_chh);
  if (site_count == 0)
    return;

  MethylRecord *buffer = calloc(site_count, sizeof(MethylRecord));
  initialize_buffer(buffer, site_count, targ->chr_seq, targ->chr_len,
                    targ->keep_chg, targ->keep_chh);

  uint32_t *site_positions = malloc(site_count * sizeof(uint32_t));
  for (size_t i = 0; i < site_count; i++)
    site_positions[i] = buffer[i].pos - 1;

  int n_threads = get_nprocs();
  if (n_threads > 32)
    n_threads = 32;
  if (n_threads < 1)
    n_threads = 8;

  size_t sites_per_thread = (site_count + n_threads - 1) / n_threads;

  RegionArg *regions = calloc(n_threads, sizeof(RegionArg));
  pthread_t *threads = malloc(n_threads * sizeof(pthread_t));

  for (int i = 0; i < n_threads; i++) {
    RegionArg *r = &regions[i];
    r->base = *targ;
    r->start_site_idx = i * sites_per_thread;
    r->end_site_idx = (i + 1) * sites_per_thread;
    if (r->end_site_idx > site_count)
      r->end_site_idx = site_count;
    if (r->start_site_idx >= site_count)
      break;

    r->base.start_pos = site_positions[r->start_site_idx];
    r->base.end_pos = (r->end_site_idx < site_count)
                          ? site_positions[r->end_site_idx - 1] + 10
                          : targ->chr_len;

    uint32_t range = r->base.end_pos - r->base.start_pos;
    r->counts.mC = calloc(range, sizeof(uint32_t));
    r->counts.uC = calloc(range, sizeof(uint32_t));
    r->counts.size = range;

    pthread_create(&threads[i], NULL, process_region_direct, r);
  }

  for (int i = 0; i < n_threads; i++) {
    if (regions[i].start_site_idx < site_count)
      pthread_join(threads[i], NULL);
  }

  for (int i = 0; i < n_threads; i++) {
    if (regions[i].start_site_idx >= site_count)
      continue;
    uint32_t offset = regions[i].base.start_pos;
    for (size_t j = regions[i].start_site_idx; j < regions[i].end_site_idx;
         j++) {
      uint32_t gpos = site_positions[j];
      uint32_t rel = gpos - offset;
      buffer[j].mC += regions[i].counts.mC[rel];
      buffer[j].uC += regions[i].counts.uC[rel];
    }
    free(regions[i].counts.mC);
    free(regions[i].counts.uC);
  }

  // Output (unchanged)
  char final_out_path[1024] = {
      0}; // Store the last written file path for compression

  if (targ->split_context_files) {
    for (int ctx = CONTEXT_CPG; ctx <= CONTEXT_CHH; ++ctx) {
      if ((ctx == CONTEXT_CPG) || (ctx == CONTEXT_CHG && targ->keep_chg) ||
          (ctx == CONTEXT_CHH && targ->keep_chh)) {
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
        const char *ext = (targ->output_format == OUTPUT_TXT) ? ".txt"
                          : (targ->output_format == OUTPUT_PARQUET)
                              ? ".parquet"
                              : ".h5"; // For OUTPUT_BOTH, use .h5
        snprintf(out_path, sizeof(out_path), "%s/%s-%s%s", targ->out_dir,
                 targ->chr, get_context_string(ctx), ext);

        flush_buffer(out_path, ctx_buffer, n_ctx_records, targ->compression,
                     targ->hdf5_chunk_size, 0, targ->min_cov, targ->cap_cov,
                     targ->output_format);

        free(ctx_buffer);

        // OPTIMIZATION: External Zstd compression
        if (targ->compression >= 9 && (targ->output_format == OUTPUT_HDF5 ||
                                       targ->output_format == OUTPUT_BOTH)) {
          char cmd[2048];
          snprintf(cmd, sizeof(cmd),
                   "zstd --ultra -9 -T0 --rm \"%s\" -o \"%s.zst.h5\" && "
                   "mv \"%s.zst.h5\" \"%s\" 2>/dev/null",
                   out_path, out_path, out_path, out_path);

          log_time("Launching multi-threaded Zstd-9 compression for %s ...\n",
                   out_path);
          int ret = system(cmd);
          if (ret == 0)
            log_time(
                "Finished ultra-compression of %s (Zstd-9, multi-threaded)\n",
                out_path);
          else
            log_time("Warning: external Zstd failed for %s (you can compress "
                     "manually)\n",
                     out_path);
        }
      }
    }
  } else {
    // Output file name - use appropriate extension based on output format
    char out_path[1024];
    const char *ext = (targ->output_format == OUTPUT_TXT) ? ".txt"
                      : (targ->output_format == OUTPUT_PARQUET)
                          ? ".parquet"
                          : ".h5"; // For OUTPUT_BOTH, use .h5
    snprintf(out_path, sizeof(out_path), "%s/%s%s", targ->out_dir, targ->chr,
             ext);
    flush_buffer(out_path, buffer, site_count, targ->compression,
                 targ->hdf5_chunk_size, 0, targ->min_cov, targ->cap_cov,
                 targ->output_format);

    // OPTIMIZATION: External Zstd compression
    if (targ->compression >= 9 && (targ->output_format == OUTPUT_HDF5 ||
                                   targ->output_format == OUTPUT_BOTH)) {
      char cmd[2048];
      snprintf(cmd, sizeof(cmd),
               "zstd --ultra -9 -T0 --rm \"%s\" -o \"%s.zst.h5\" && "
               "mv \"%s.zst.h5\" \"%s\" 2>/dev/null",
               out_path, out_path, out_path, out_path);

      log_time("Launching multi-threaded Zstd-9 compression for %s ...\n",
               out_path);
      int ret = system(cmd);
      if (ret == 0)
        log_time("Finished ultra-compression of %s (Zstd-9, multi-threaded)\n",
                 out_path);
      else
        log_time("Warning: external Zstd failed for %s (you can compress "
                 "manually)\n",
                 out_path);
    }
  }

  free(site_positions);
  free(buffer);
  free(threads);
  free(regions);
}

int load_chrom_mapping(const char *filename, ChromMapEntry **entries,
                       int *n_entries, char **reference_file) {
  FILE *fp = fopen(filename, "r");
  if (!fp)
    return -1;
  fseek(fp, 0, SEEK_END);
  long len = ftell(fp);
  fseek(fp, 0, SEEK_SET);
  char *data = malloc(len + 1);
  size_t nread = fread(data, 1, len, fp);
  if (nread != len) {
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
  if (ref && cJSON_IsString(ref) && ref->valuestring) {
    *reference_file = strdup(ref->valuestring);
    fprintf(stderr, "Found reference file: %s\n", *reference_file);
  } else
    fprintf(stderr, "Warning: No reference file specified in mapping\n");

  cJSON *chroms = cJSON_GetObjectItem(json, "chromosomes");
  if (!chroms || !cJSON_IsArray(chroms)) {
    cJSON_Delete(json);
    return -4;
  }
  int count = cJSON_GetArraySize(chroms);
  *entries = calloc(count, sizeof(ChromMapEntry));
  *n_entries = 0;
  for (int i = 0; i < count; ++i) {
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
    if (fasta && cJSON_IsString(fasta) && fasta->valuestring) {
      strncpy(e->fasta, fasta->valuestring, sizeof(e->fasta) - 1);
      e->fasta[sizeof(e->fasta) - 1] = '\0';
    }

    if (bam && cJSON_IsString(bam) && bam->valuestring) {
      strncpy(e->bam, bam->valuestring, sizeof(e->bam) - 1);
      e->bam[sizeof(e->bam) - 1] = '\0';
    }

    if (name && cJSON_IsString(name) && name->valuestring) {
      strncpy(e->name, name->valuestring, sizeof(e->name) - 1);
      e->name[sizeof(e->name) - 1] = '\0';
    }

    // Validate that we have all required fields
    if (e->fasta[0] == '\0' || e->bam[0] == '\0' || e->name[0] == '\0') {
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
