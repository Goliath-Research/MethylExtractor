#include "methyl_extractor.h"

#define MHAP_READ_CAP0 1024
#define MHAP_OBS_CAP0 4096

void mhap_store_init(MhapStore *s)
{
    memset(s, 0, sizeof(*s));
}

void mhap_store_free(MhapStore *s)
{
    if (!s)
        return;
    free(s->read_start);
    free(s->read_strand);
    free(s->n_cpg);
    free(s->cpg_pos);
    free(s->meth);
    memset(s, 0, sizeof(*s));
}

static int grow_reads(MhapStore *s, size_t need)
{
    if (s->n_reads + need <= s->cap_reads)
        return 0;
    size_t cap = s->cap_reads ? s->cap_reads : MHAP_READ_CAP0;
    while (cap < s->n_reads + need)
        cap *= 2;
    int32_t *rs = realloc(s->read_start, cap * sizeof(int32_t));
    int8_t *st = realloc(s->read_strand, cap * sizeof(int8_t));
    int32_t *nc = realloc(s->n_cpg, cap * sizeof(int32_t));
    if (!rs || !st || !nc)
    {
        free(rs);
        free(st);
        free(nc);
        return -1;
    }
    s->read_start = rs;
    s->read_strand = st;
    s->n_cpg = nc;
    s->cap_reads = cap;
    return 0;
}

static int grow_obs(MhapStore *s, size_t need)
{
    if (s->n_obs + need <= s->cap_obs)
        return 0;
    size_t cap = s->cap_obs ? s->cap_obs : MHAP_OBS_CAP0;
    while (cap < s->n_obs + need)
        cap *= 2;
    int32_t *pos = realloc(s->cpg_pos, cap * sizeof(int32_t));
    uint8_t *meth = realloc(s->meth, cap * sizeof(uint8_t));
    if (!pos || !meth)
    {
        free(pos);
        free(meth);
        return -1;
    }
    s->cpg_pos = pos;
    s->meth = meth;
    s->cap_obs = cap;
    return 0;
}

int mhap_store_add(MhapStore *s, int32_t start, int8_t strand,
                   const int32_t *pos, const uint8_t *meth, int n)
{
    if (!s || n <= 0 || !pos || !meth)
        return 0;
    if (grow_reads(s, 1) != 0 || grow_obs(s, (size_t)n) != 0)
        return -1;
    s->read_start[s->n_reads] = start;
    s->read_strand[s->n_reads] = strand;
    s->n_cpg[s->n_reads] = n;
    memcpy(s->cpg_pos + s->n_obs, pos, (size_t)n * sizeof(int32_t));
    memcpy(s->meth + s->n_obs, meth, (size_t)n * sizeof(uint8_t));
    s->n_obs += (size_t)n;
    s->n_reads++;
    return 0;
}

int mhap_store_merge(MhapStore *dst, const MhapStore *src)
{
    if (!dst || !src || src->n_reads == 0)
        return 0;
    if (grow_reads(dst, src->n_reads) != 0 || grow_obs(dst, src->n_obs) != 0)
        return -1;
    memcpy(dst->read_start + dst->n_reads, src->read_start,
           src->n_reads * sizeof(int32_t));
    memcpy(dst->read_strand + dst->n_reads, src->read_strand,
           src->n_reads * sizeof(int8_t));
    memcpy(dst->n_cpg + dst->n_reads, src->n_cpg, src->n_reads * sizeof(int32_t));
    memcpy(dst->cpg_pos + dst->n_obs, src->cpg_pos, src->n_obs * sizeof(int32_t));
    memcpy(dst->meth + dst->n_obs, src->meth, src->n_obs * sizeof(uint8_t));
    dst->n_reads += src->n_reads;
    dst->n_obs += src->n_obs;
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
        hsize_t chunk[1] = {dims[0] > 0 ? dims[0] : 1};
        if (H5Pset_chunk(dcpl, rank, chunk) >= 0)
        {
            htri_t zstd_avail = H5Zfilter_avail(ZSTD_FILTER);
            if (zstd_avail > 0)
            {
                unsigned int cd_values[1] = {(unsigned int)compression};
                if (H5Pset_filter(dcpl, ZSTD_FILTER, H5Z_FLAG_MANDATORY, 1,
                                  cd_values) < 0)
                {
                    int gzip_level = compression > 9 ? 9 : compression;
                    H5Pset_deflate(dcpl, gzip_level);
                }
            }
            else
            {
                int gzip_level = compression > 9 ? 9 : compression;
                H5Pset_deflate(dcpl, gzip_level);
            }
        }
    }
    return dcpl;
}

static int write_i32_dataset(hid_t grp, const char *name, const int32_t *data,
                             hsize_t len, int compression)
{
    hsize_t dims[1] = {len > 0 ? len : 1};
    hid_t space = H5Screate_simple(1, dims, NULL);
    if (space < 0)
        return -1;
    hid_t dcpl = make_chunked_dcpl(compression, 1, dims);
    hid_t dset = H5Dcreate2(grp, name, H5T_STD_I32LE, space, H5P_DEFAULT, dcpl,
                            H5P_DEFAULT);
    if (dset >= 0)
    {
        if (len > 0)
            H5Dwrite(dset, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        H5Dclose(dset);
    }
    H5Sclose(space);
    if (dcpl >= 0)
        H5Pclose(dcpl);
    return dset >= 0 ? 0 : -1;
}

static int write_i8_dataset(hid_t grp, const char *name, const int8_t *data,
                            hsize_t len, int compression)
{
    hsize_t dims[1] = {len > 0 ? len : 1};
    hid_t space = H5Screate_simple(1, dims, NULL);
    if (space < 0)
        return -1;
    hid_t dcpl = make_chunked_dcpl(compression, 1, dims);
    hid_t dset = H5Dcreate2(grp, name, H5T_STD_I8LE, space, H5P_DEFAULT, dcpl,
                            H5P_DEFAULT);
    if (dset >= 0)
    {
        if (len > 0)
            H5Dwrite(dset, H5T_NATIVE_INT8, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        H5Dclose(dset);
    }
    H5Sclose(space);
    if (dcpl >= 0)
        H5Pclose(dcpl);
    return dset >= 0 ? 0 : -1;
}

static int write_u8_dataset(hid_t grp, const char *name, const uint8_t *data,
                            hsize_t len, int compression)
{
    hsize_t dims[1] = {len > 0 ? len : 1};
    hid_t space = H5Screate_simple(1, dims, NULL);
    if (space < 0)
        return -1;
    hid_t dcpl = make_chunked_dcpl(compression, 1, dims);
    hid_t dset = H5Dcreate2(grp, name, H5T_STD_U8LE, space, H5P_DEFAULT, dcpl,
                            H5P_DEFAULT);
    if (dset >= 0)
    {
        if (len > 0)
            H5Dwrite(dset, H5T_NATIVE_UINT8, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        H5Dclose(dset);
    }
    H5Sclose(space);
    if (dcpl >= 0)
        H5Pclose(dcpl);
    return dset >= 0 ? 0 : -1;
}

int write_mhap_h5(const char *filename, const char *chrom, const char *context,
                  int compression, const MhapStore *store)
{
    if (!filename || !chrom || !context || !store || store->n_reads == 0)
        return 0;

    export_lock();
    hid_t file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0)
    {
        export_unlock();
        return -1;
    }
    hid_t grp = H5Gcreate2(file, MHAP_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (grp < 0)
    {
        H5Fclose(file);
        export_unlock();
        return -1;
    }

    write_str_attr(grp, "schema_version", MHAP_SCHEMA_VERSION);
    write_str_attr(grp, "context", context);
    write_str_attr(grp, "chrom", chrom);

    write_i32_dataset(grp, "read_start", store->read_start, store->n_reads,
                      compression);
    write_i8_dataset(grp, "read_strand", store->read_strand, store->n_reads,
                     compression);
    write_i32_dataset(grp, "n_cpg", store->n_cpg, store->n_reads, compression);
    write_i32_dataset(grp, "cpg_pos", store->cpg_pos, store->n_obs, compression);
    write_u8_dataset(grp, "meth", store->meth, store->n_obs, compression);

    H5Gclose(grp);
    H5Fclose(file);
    export_unlock();
    return 0;
}
