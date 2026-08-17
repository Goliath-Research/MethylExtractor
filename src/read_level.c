#include "methyl_extractor.h"

static int site_matches_context(const char *chr_seq, uint32_t chr_len, int pos,
                                int context, int *canonical_pos)
{
    if (pos < 0 || (uint32_t)pos >= chr_len)
        return 0;

    if (context == CONTEXT_CPG || context == CONTEXT_CHG)
    {
        int cpg = (context == CONTEXT_CPG) ? isCpG(chr_seq, pos, (int)chr_len)
                                           : isCHG(chr_seq, pos, (int)chr_len);
        if (cpg == 1)
        {
            if (canonical_pos)
                *canonical_pos = pos;
            return 1;
        }
        if (cpg == -1)
        {
            if (canonical_pos)
                *canonical_pos = pos - 1;
            return 1;
        }
        return 0;
    }

    if (context == CONTEXT_CHH)
    {
        int chh = isCHH(chr_seq, pos, (int)chr_len);
        if (chh == 1 || chh == -1)
        {
            if (canonical_pos)
                *canonical_pos = pos;
            return 1;
        }
    }
    return 0;
}

static int collect_context_sites(const char *chr_seq, uint32_t chr_len, int context,
                                 uint32_t **out_pos, size_t *out_n)
{
    size_t cap = 1024;
    size_t n = 0;
    uint32_t *pos = malloc(cap * sizeof(uint32_t));
    if (!pos)
        return -1;

    for (uint32_t p = 0; p < chr_len; ++p)
    {
        int canonical = 0;
        if (!site_matches_context(chr_seq, chr_len, (int)p, context, &canonical))
            continue;
        if (canonical < 0)
            continue;

        if (context == CONTEXT_CPG || context == CONTEXT_CHG)
        {
            if (n > 0 && pos[n - 1] == (uint32_t)canonical)
                continue;
        }

        if (n >= cap)
        {
            cap *= 2;
            uint32_t *tmp = realloc(pos, cap * sizeof(uint32_t));
            if (!tmp)
            {
                free(pos);
                return -1;
            }
            pos = tmp;
        }
        pos[n++] = (uint32_t)canonical;
    }

    *out_pos = pos;
    *out_n = n;
    return 0;
}

void free_read_level_tiles(ReadLevelTiles *tiles)
{
    if (!tiles)
        return;
    free(tiles->cpg_pos);
    free(tiles->pos_to_cpg);
    free(tiles->tile_start_pos);
    free(tiles->tile_cpg_positions);
    memset(tiles, 0, sizeof(*tiles));
}

int build_read_level_tiles(const char *chr_seq, uint32_t chr_len, int context,
                           int tile_size, ReadLevelTiles *out)
{
    if (!chr_seq || !out || tile_size < MIN_TILE_SIZE || tile_size > MAX_TILE_SIZE)
        return -1;

    memset(out, 0, sizeof(*out));
    out->context = context;
    out->tile_size = tile_size;

    uint32_t *cpg_pos = NULL;
    size_t n_cpg = 0;
    if (collect_context_sites(chr_seq, chr_len, context, &cpg_pos, &n_cpg) != 0)
        return -1;

    out->n_cpg = n_cpg;
    out->cpg_pos = cpg_pos;
    if (n_cpg == 0)
        return 0;

    out->pos_to_cpg = malloc(chr_len * sizeof(int));
    if (!out->pos_to_cpg)
    {
        free_read_level_tiles(out);
        return -1;
    }
    for (uint32_t i = 0; i < chr_len; ++i)
        out->pos_to_cpg[i] = -1;
    for (size_t i = 0; i < n_cpg; ++i)
        out->pos_to_cpg[cpg_pos[i]] = (int)i;

    out->n_tiles = n_cpg / (size_t)tile_size;
    if (out->n_tiles == 0)
        return 0;

    out->tile_start_pos = malloc(out->n_tiles * sizeof(uint32_t));
    out->tile_cpg_positions =
        malloc(out->n_tiles * (size_t)tile_size * sizeof(uint32_t));
    if (!out->tile_start_pos || !out->tile_cpg_positions)
    {
        free_read_level_tiles(out);
        return -1;
    }

    for (size_t t = 0; t < out->n_tiles; ++t)
    {
        size_t base = t * (size_t)tile_size;
        out->tile_start_pos[t] = cpg_pos[base] + 1;
        for (int j = 0; j < tile_size; ++j)
            out->tile_cpg_positions[t * (size_t)tile_size + (size_t)j] =
                cpg_pos[base + (size_t)j] + 1;
    }

    return 0;
}

size_t read_level_tile_range_for_region(const ReadLevelTiles *tiles,
                                        uint32_t region_start,
                                        uint32_t region_end,
                                        size_t *out_start, size_t *out_end)
{
    if (!tiles || tiles->n_tiles == 0)
    {
        *out_start = 0;
        *out_end = 0;
        return 0;
    }

    size_t first = tiles->n_tiles;
    size_t last = 0;
    for (size_t t = 0; t < tiles->n_tiles; ++t)
    {
        uint32_t tile_pos0 = tiles->cpg_pos[t * (size_t)tiles->tile_size];
        if (tile_pos0 >= region_start && tile_pos0 < region_end)
        {
            if (t < first)
                first = t;
            last = t + 1;
        }
    }

    if (first >= tiles->n_tiles)
    {
        *out_start = 0;
        *out_end = 0;
        return 0;
    }

    *out_start = first;
    *out_end = last;
    return last - first;
}

static int record_site_call(int8_t *state, size_t state_len, size_t state_off,
                            const ReadLevelTiles *tiles, int cpg_idx, int methylated)
{
    if (cpg_idx < 0 || (size_t)cpg_idx >= tiles->n_cpg)
        return 0;
    size_t rel = (size_t)cpg_idx - state_off;
    if (rel >= state_len)
        return 0;
    if (state[rel] >= 0 && state[rel] != methylated)
        return -1;
    state[rel] = (int8_t)methylated;
    return 0;
}

static int map_read_call_to_cpg(const char *chr_seq, uint32_t chr_len, int pos,
                                int context, int strand, int base,
                                const ReadLevelTiles *tiles, int *cpg_idx,
                                int *methylated)
{
    int canonical = pos;
    if (context == CONTEXT_CPG || context == CONTEXT_CHG)
    {
        char ref = toupper(chr_seq[pos]);
        if (ref == 'C' && (strand == 1 || strand == 3))
        {
            canonical = pos;
            if (base == 2)
                *methylated = 1;
            else if (base == 8)
                *methylated = 0;
            else
                return 0;
        }
        else if (ref == 'G' && (strand == 2 || strand == 4))
        {
            canonical = pos - 1;
            if (canonical < 0)
                return 0;
            if (base == 4)
                *methylated = 1;
            else if (base == 1)
                *methylated = 0;
            else
                return 0;
        }
        else
            return 0;
    }
    else if (context == CONTEXT_CHH)
    {
        char ref = toupper(chr_seq[pos]);
        if (ref == 'C' && (strand == 1 || strand == 3))
        {
            if (base == 2)
                *methylated = 1;
            else if (base == 8)
                *methylated = 0;
            else
                return 0;
        }
        else if (ref == 'G' && (strand == 2 || strand == 4))
        {
            if (base == 4)
                *methylated = 1;
            else if (base == 1)
                *methylated = 0;
            else
                return 0;
        }
        else
            return 0;
    }
    else
        return 0;

    if (canonical < 0 || (uint32_t)canonical >= chr_len)
        return 0;
    int idx = tiles->pos_to_cpg[canonical];
    if (idx < 0)
        return 0;
    *cpg_idx = idx;
    return 1;
}

void read_level_accumulate(const ReadLevelTiles *tiles, uint32_t *hist,
                           size_t tile_start, size_t tile_end,
                           const char *chr_seq, uint8_t *seq, uint8_t *qual,
                           uint32_t *cigar, int n_cigar, int32_t read_start,
                           int32_t clip_from, int strand, int min_phred,
                           uint32_t chr_len)
{
    if (!tiles || !hist || tile_start >= tile_end || tiles->n_tiles == 0)
        return;

    int k = tiles->tile_size;
    int n_patterns = 1 << k;
    size_t cpg_start = tile_start * (size_t)k;
    size_t cpg_end = tile_end * (size_t)k;
    if (cpg_end > tiles->n_cpg)
        cpg_end = tiles->n_cpg;
    size_t state_len = cpg_end - cpg_start;
    if (state_len == 0)
        return;

    int8_t *state = calloc(state_len, sizeof(int8_t));
    if (!state)
        return;
    for (size_t i = 0; i < state_len; ++i)
        state[i] = -1;

    int32_t pos = read_start;
    int qpos = 0;
    int conflict = 0;

    for (int c = 0; c < n_cigar && !conflict; ++c)
    {
        int op = bam_cigar_op(cigar[c]);
        int len = bam_cigar_oplen(cigar[c]);

        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
        {
            for (int i = 0; i < len; ++i)
            {
                if (pos >= clip_from)
                {
                    pos++;
                    qpos++;
                    continue;
                }
                if (pos < 0 || (uint32_t)pos >= chr_len)
                {
                    pos++;
                    qpos++;
                    continue;
                }
                if (qual[qpos] < min_phred)
                {
                    pos++;
                    qpos++;
                    continue;
                }

                char ref = toupper(chr_seq[pos]);
                if (ref != 'C' && ref != 'G')
                {
                    pos++;
                    qpos++;
                    continue;
                }

                int base = bam_seqi(seq, qpos);
                int cpg_idx = 0;
                int methylated = 0;
                if (map_read_call_to_cpg(chr_seq, chr_len, pos, tiles->context,
                                         strand, base, tiles, &cpg_idx,
                                         &methylated))
                {
                    if (cpg_idx < (int)cpg_start || cpg_idx >= (int)cpg_end)
                    {
                        pos++;
                        qpos++;
                        continue;
                    }
                    if (record_site_call(state, state_len, cpg_start, tiles,
                                         cpg_idx, methylated) != 0)
                    {
                        conflict = 1;
                        break;
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

    if (!conflict)
    {
        for (size_t t = tile_start; t < tile_end; ++t)
        {
            size_t base = t * (size_t)k;
            if (base + (size_t)k > tiles->n_cpg)
                break;

            int pattern = 0;
            int complete = 1;
            for (int j = 0; j < k; ++j)
            {
                size_t rel = base + (size_t)j - cpg_start;
                if (state[rel] < 0)
                {
                    complete = 0;
                    break;
                }
                if (state[rel] > 0)
                    pattern |= (1 << (k - 1 - j));
            }
            if (complete)
                hist[t * (size_t)n_patterns + (size_t)pattern]++;
        }
    }

    free(state);
}

static int set_hdf5_compression(hid_t dcpl, int compression)
{
    if (compression <= 0)
        return 0;

    htri_t zstd_avail = H5Zfilter_avail(ZSTD_FILTER);
    if (zstd_avail > 0)
    {
        unsigned int cd_values[1] = {(unsigned int)compression};
        if (H5Pset_filter(dcpl, ZSTD_FILTER, H5Z_FLAG_MANDATORY, 1, cd_values) >= 0)
            return 0;
    }

    int gzip_level = compression > 9 ? 9 : compression;
    if (H5Pset_deflate(dcpl, gzip_level) < 0)
        return -1;
    return 0;
}

static int write_str_attr(hid_t obj, const char *name, const char *value)
{
    hid_t type = H5Tcopy(H5T_C_S1);
    size_t n = strlen(value) + 1;
    H5Tset_size(type, n);
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate2(obj, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr >= 0)
    {
        H5Awrite(attr, type, value);
        H5Aclose(attr);
    }
    H5Tclose(type);
    H5Sclose(space);
    return attr >= 0 ? 0 : -1;
}

static hid_t make_chunked_dcpl(int compression, int rank, const hsize_t *dims)
{
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    if (dcpl < 0)
        return dcpl;
    if (compression > 0 && rank > 0 && dims)
    {
        hsize_t chunk[3] = {1, 1, 1};
        for (int i = 0; i < rank; ++i)
            chunk[i] = dims[i] > 0 ? dims[i] : 1;
        if (H5Pset_chunk(dcpl, rank, chunk) >= 0)
            set_hdf5_compression(dcpl, compression);
    }
    return dcpl;
}

static int write_u32_dataset(hid_t grp, const char *name, const uint32_t *data,
                             hsize_t len, int compression)
{
    hsize_t dims[1] = {len > 0 ? len : 1};
    hid_t space = H5Screate_simple(1, dims, NULL);
    if (space < 0)
        return -1;
    hid_t dcpl = make_chunked_dcpl(compression, 1, dims);
    hid_t dset = H5Dcreate2(grp, name, H5T_STD_U32LE, space, H5P_DEFAULT, dcpl,
                            H5P_DEFAULT);
    if (dset >= 0)
    {
        if (len > 0)
            H5Dwrite(dset, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     data);
        H5Dclose(dset);
    }
    H5Sclose(space);
    if (dcpl >= 0)
        H5Pclose(dcpl);
    return dset >= 0 ? 0 : -1;
}

static int write_u16_dataset(hid_t grp, const char *name, const uint16_t *data,
                               hsize_t len, int compression)
{
    hsize_t dims[1] = {len > 0 ? len : 1};
    hid_t space = H5Screate_simple(1, dims, NULL);
    if (space < 0)
        return -1;
    hid_t dcpl = make_chunked_dcpl(compression, 1, dims);
    hid_t dset = H5Dcreate2(grp, name, H5T_STD_U16LE, space, H5P_DEFAULT, dcpl,
                            H5P_DEFAULT);
    if (dset >= 0)
    {
        if (len > 0)
            H5Dwrite(dset, H5T_NATIVE_UINT16, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     data);
        H5Dclose(dset);
    }
    H5Sclose(space);
    if (dcpl >= 0)
        H5Pclose(dcpl);
    return dset >= 0 ? 0 : -1;
}

int write_read_level_patterns_h5(const char *filename, const char *context,
                                 int tile_size, int min_tile_reads,
                                 int compression, const ReadLevelTiles *tiles,
                                 const uint32_t *hist)
{
    if (!filename || !context || !tiles || !hist || tiles->n_tiles == 0)
        return 0;

    int k = tile_size;
    int n_patterns = 1 << k;
    size_t n_emit = 0;
    for (size_t t = 0; t < tiles->n_tiles; ++t)
    {
        uint32_t total = 0;
        for (int p = 0; p < n_patterns; ++p)
            total += hist[t * (size_t)n_patterns + (size_t)p];
        if (total >= (uint32_t)min_tile_reads)
            n_emit++;
    }

    if (n_emit == 0)
        return 0;

    size_t nnz_cap = n_emit * (size_t)n_patterns;
    uint32_t *pattern_tile_id = malloc(nnz_cap * sizeof(uint32_t));
    uint16_t *pattern_id = malloc(nnz_cap * sizeof(uint16_t));
    uint32_t *pattern_count = malloc(nnz_cap * sizeof(uint32_t));
    uint32_t *tile_start_pos = malloc(n_emit * sizeof(uint32_t));
    uint32_t *tile_cpg_positions =
        malloc(n_emit * (size_t)k * sizeof(uint32_t));
    uint32_t *tile_n_reads = malloc(n_emit * sizeof(uint32_t));
    if (!pattern_tile_id || !pattern_id || !pattern_count || !tile_start_pos ||
        !tile_cpg_positions || !tile_n_reads)
    {
        free(pattern_tile_id);
        free(pattern_id);
        free(pattern_count);
        free(tile_start_pos);
        free(tile_cpg_positions);
        free(tile_n_reads);
        return -1;
    }

    size_t emit_idx = 0;
    size_t nnz = 0;
    for (size_t t = 0; t < tiles->n_tiles; ++t)
    {
        uint32_t total = 0;
        for (int p = 0; p < n_patterns; ++p)
            total += hist[t * (size_t)n_patterns + (size_t)p];
        if (total < (uint32_t)min_tile_reads)
            continue;

        tile_start_pos[emit_idx] = tiles->tile_start_pos[t];
        tile_n_reads[emit_idx] = total;
        for (int j = 0; j < k; ++j)
            tile_cpg_positions[emit_idx * (size_t)k + (size_t)j] =
                tiles->tile_cpg_positions[t * (size_t)k + (size_t)j];

        for (int p = 0; p < n_patterns; ++p)
        {
            uint32_t cnt = hist[t * (size_t)n_patterns + (size_t)p];
            if (cnt == 0)
                continue;
            pattern_tile_id[nnz] = (uint32_t)emit_idx;
            pattern_id[nnz] = (uint16_t)p;
            pattern_count[nnz] = cnt;
            nnz++;
        }
        emit_idx++;
    }

    export_lock();
    hid_t file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0)
        goto fail;

    hid_t grp = H5Gcreate2(file, READ_LEVEL_GROUP, H5P_DEFAULT, H5P_DEFAULT,
                           H5P_DEFAULT);
    if (grp < 0)
        goto fail_close;

    {
        hid_t space = H5Screate(H5S_SCALAR);
        write_str_attr(grp, "context", context);

        int32_t val = tile_size;
        hid_t attr = H5Acreate2(grp, "tile_size", H5T_STD_I32LE, space,
                                H5P_DEFAULT, H5P_DEFAULT);
        if (attr >= 0)
        {
            H5Awrite(attr, H5T_NATIVE_INT, &val);
            H5Aclose(attr);
        }

        val = min_tile_reads;
        attr = H5Acreate2(grp, "min_tile_reads", H5T_STD_I32LE, space,
                          H5P_DEFAULT, H5P_DEFAULT);
        if (attr >= 0)
        {
            H5Awrite(attr, H5T_NATIVE_INT, &val);
            H5Aclose(attr);
        }
        H5Sclose(space);
    }

    write_str_attr(grp, "pattern_encoding", READ_LEVEL_PATTERN_ENCODING);
    write_str_attr(grp, "schema_version", READ_LEVEL_SCHEMA_VERSION);

    write_u32_dataset(grp, "tile_start_pos", tile_start_pos, n_emit, compression);
    write_u32_dataset(grp, "tile_n_reads", tile_n_reads, n_emit, compression);

    {
        hsize_t dims2[2] = {n_emit, (hsize_t)k};
        hid_t space2 = H5Screate_simple(2, dims2, NULL);
        hid_t dcpl2 = make_chunked_dcpl(compression, 2, dims2);
        hid_t dset = H5Dcreate2(grp, "tile_cpg_positions", H5T_STD_U32LE, space2,
                                H5P_DEFAULT, dcpl2, H5P_DEFAULT);
        if (dset >= 0)
        {
            H5Dwrite(dset, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     tile_cpg_positions);
            H5Dclose(dset);
        }
        H5Sclose(space2);
        if (dcpl2 >= 0)
            H5Pclose(dcpl2);
    }

    if (nnz > 0)
    {
        write_u32_dataset(grp, "pattern_tile_id", pattern_tile_id, nnz, compression);
        write_u16_dataset(grp, "pattern_id", pattern_id, nnz, compression);
        write_u32_dataset(grp, "pattern_count", pattern_count, nnz, compression);
    }

    H5Gclose(grp);
    H5Fclose(file);
    export_unlock();

    free(pattern_tile_id);
    free(pattern_id);
    free(pattern_count);
    free(tile_start_pos);
    free(tile_cpg_positions);
    free(tile_n_reads);
    log_time("Wrote read-level patterns: %s (%zu tiles, %zu nnz)\n", filename,
             n_emit, nnz);
    return 0;

fail_close:
    H5Fclose(file);
fail:
    export_unlock();
    free(pattern_tile_id);
    free(pattern_id);
    free(pattern_count);
    free(tile_start_pos);
    free(tile_cpg_positions);
    free(tile_n_reads);
    return -1;
}

void free_read_level_chrom_data(ReadLevelChromData *data)
{
    if (!data)
        return;
    for (int i = 0; i < data->n_active; ++i)
    {
        free_read_level_tiles(&data->ctx[i].tiles);
        free(data->ctx[i].hist);
        data->ctx[i].hist = NULL;
        data->ctx[i].active = 0;
    }
    data->n_active = 0;
}
