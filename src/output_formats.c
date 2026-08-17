#include "methyl_extractor.h"

void merge_filter_stats(FilterStats *dst, const FilterStats *src)
{
    if (!dst || !src)
        return;
    dst->reads_seen += src->reads_seen;
    dst->reads_used += src->reads_used;
    dst->reads_dropped_unmapped += src->reads_dropped_unmapped;
    dst->reads_dropped_secondary += src->reads_dropped_secondary;
    dst->reads_dropped_qc_fail += src->reads_dropped_qc_fail;
    dst->reads_dropped_duplicate += src->reads_dropped_duplicate;
    dst->reads_dropped_supplementary += src->reads_dropped_supplementary;
    dst->reads_dropped_low_mapq += src->reads_dropped_low_mapq;
    dst->reads_dropped_multimap += src->reads_dropped_multimap;
    dst->reads_dropped_no_strand += src->reads_dropped_no_strand;
    dst->bases_skipped_overlap_clip += src->bases_skipped_overlap_clip;
    dst->bases_skipped_low_phred += src->bases_skipped_low_phred;
    dst->bases_skipped_non_cytosine += src->bases_skipped_non_cytosine;
    dst->bases_counted_methylated += src->bases_counted_methylated;
    dst->bases_counted_unmethylated += src->bases_counted_unmethylated;
}

void coverage_histogram_init(CoverageHistogram *hist)
{
    if (!hist)
        return;
    memset(hist, 0, sizeof(*hist));
}

void coverage_histogram_add(CoverageHistogram *hist, uint32_t coverage)
{
    if (!hist)
        return;
    if (coverage >= COV_HIST_BINS)
        hist->overflow++;
    else
        hist->bins[coverage]++;
    hist->total_sites++;
}

double coverage_histogram_percentile(const CoverageHistogram *hist, double percentile)
{
    if (!hist || hist->total_sites == 0)
        return 0.0;
    if (percentile < 0.0)
        percentile = 0.0;
    if (percentile > 1.0)
        percentile = 1.0;

    uint64_t target = (uint64_t)ceil(percentile * (double)hist->total_sites);
    if (target == 0)
        target = 1;

    uint64_t cumulative = 0;
    for (int i = 0; i < COV_HIST_BINS; i++)
    {
        cumulative += hist->bins[i];
        if (cumulative >= target)
            return (double)i;
    }
    return (double)COV_HIST_BINS;
}

MethylStats calculate_statistics(MethylRecord *filtered_buffer, size_t n_records)
{
    MethylStats stats;
    memset(&stats, 0, sizeof(stats));

    if (n_records == 0)
        return stats;

    stats.num_positions = n_records;
    stats.sites_in_reference = n_records;

    uint64_t total_mC = 0;
    uint64_t total_uC = 0;

    for (size_t i = 0; i < n_records; i++)
    {
        total_mC += filtered_buffer[i].mC;
        total_uC += filtered_buffer[i].uC;
    }

    stats.total_methylated = total_mC;
    stats.total_unmethylated = total_uC;

    uint64_t total_coverage = total_mC + total_uC;
    if (total_coverage > 0)
    {
        stats.avg_methylation_level = (double)total_mC / (double)total_coverage;
        stats.avg_coverage = (double)total_coverage / (double)n_records;
    }

    return stats;
}

MethylStats analyze_buffer(MethylRecord *buffer, size_t n_records, int min_cov)
{
    MethylStats stats;
    memset(&stats, 0, sizeof(stats));

    if (n_records == 0)
        return stats;

    stats.sites_in_reference = n_records;

    CoverageHistogram hist;
    coverage_histogram_init(&hist);

    uint64_t total_coverage_all = 0;
    uint64_t plus_mC = 0, plus_uC = 0, minus_mC = 0, minus_uC = 0;

    for (size_t i = 0; i < n_records; i++)
    {
        uint32_t cov = (uint32_t)buffer[i].mC + (uint32_t)buffer[i].uC;
        coverage_histogram_add(&hist, cov);
        total_coverage_all += cov;

        if (cov > 0)
            stats.sites_with_any_coverage++;
        if (cov > 0 && cov < (uint32_t)min_cov)
            stats.sites_below_min_cov++;

        if (buffer[i].tnc.strand)
        {
            minus_mC += buffer[i].mC;
            minus_uC += buffer[i].uC;
        }
        else
        {
            plus_mC += buffer[i].mC;
            plus_uC += buffer[i].uC;
        }
    }

    stats.avg_coverage_all_sites =
        (double)total_coverage_all / (double)n_records;
    stats.coverage_p10 = coverage_histogram_percentile(&hist, 0.10);
    stats.coverage_median = coverage_histogram_percentile(&hist, 0.50);
    stats.coverage_p90 = coverage_histogram_percentile(&hist, 0.90);

    stats.methylated_plus = plus_mC;
    stats.unmethylated_plus = plus_uC;
    stats.methylated_minus = minus_mC;
    stats.unmethylated_minus = minus_uC;

    uint64_t plus_total = plus_mC + plus_uC;
    uint64_t minus_total = minus_mC + minus_uC;
    if (plus_total > 0)
        stats.methylation_plus = (double)plus_mC / (double)plus_total;
    if (minus_total > 0)
        stats.methylation_minus = (double)minus_mC / (double)minus_total;

    return stats;
}

size_t flush_buffer(const char *filename, MethylRecord *buffer,
                    size_t n_records, int compression, int chunk_size,
                    int append_mode, int min_cov, int cap_cov,
                    OutputFormat output_format, const ExtractionMeta *meta,
                    const FilterStats *filter_stats, MethylStats *out_stats)
{
    size_t records_written = 0;
    hsize_t dims[1] = {0};
    FILE *txt_fp = NULL;
    hid_t file = -1, dataset = -1, space = -1, type = -1, mem_type = -1,
          dcpl = -1;
    herr_t status = -1;
    int hdf5_locked = 0;

    MethylStats pre_stats = analyze_buffer(buffer, n_records, min_cov);

    double avg_cov = 0.0;
    for (size_t i = 0; i < n_records; i++)
    {
        int total = buffer[i].mC + buffer[i].uC;
        if (total >= min_cov)
        {
            dims[0]++;
            avg_cov += ((double)total - avg_cov) / dims[0];
        }
    }

    if (dims[0] == 0)
    {
        log_time("No records passed the coverage filter\n");
        if (out_stats)
        {
            *out_stats = pre_stats;
            out_stats->num_positions = 0;
        }
        return 0;
    }

    MethylRecord *filtered_buffer = malloc(dims[0] * sizeof(MethylRecord));
    if (!filtered_buffer)
    {
        fprintf(stderr, "Failed to allocate memory for filtered_buffer\n");
        goto cleanup;
    }

    size_t sites_capped = 0;
    size_t j = 0;
    for (size_t i = 0; i < n_records; i++)
    {
        int coverage = buffer[i].mC + buffer[i].uC;
        if (coverage >= min_cov)
        {
            if (cap_cov && coverage > avg_cov)
            {
                double prop = (double)buffer[i].mC / coverage;
                buffer[i].mC = (uint16_t)round(avg_cov * prop);
                buffer[i].uC = (uint16_t)(avg_cov - buffer[i].mC);
                sites_capped++;
            }
            filtered_buffer[j++] = buffer[i];
        }
    }

    if (output_format == OUTPUT_TXT || output_format == OUTPUT_BOTH)
    {
        char txt_filename[1024];
        const char *chr_num = strrchr(filename, '/');
        const char *dir_end = chr_num;
        if (chr_num)
            chr_num++;
        else
            chr_num = filename;

        if (dir_end)
        {
            snprintf(txt_filename, sizeof(txt_filename), "%.*s%.*s.txt",
                     (int)(dir_end - filename + 1), filename,
                     (int)(strrchr(chr_num, '.') - chr_num), chr_num);
        }
        else
        {
            snprintf(txt_filename, sizeof(txt_filename), "%.*s.txt",
                     (int)(strrchr(chr_num, '.') - chr_num), chr_num);
        }

        log_time("Starting to write text file: %s\n", txt_filename);
        txt_fp = fopen(txt_filename, "w");
        if (!txt_fp)
        {
            fprintf(stderr, "Failed to open text file for writing: %s\n",
                    txt_filename);
            goto cleanup;
        }

        for (size_t i = 0; i < dims[0]; i++)
        {
            char tnc_str[4];
            decode_trinucleotide(filtered_buffer[i].tnc.tnc, tnc_str);
            char strand = filtered_buffer[i].tnc.strand ? '-' : '+';
            int context = filtered_buffer[i].tnc.context;
            fprintf(txt_fp, "%u\t%c\t%u\t%u\t%s\t%s\n", filtered_buffer[i].pos,
                    strand, filtered_buffer[i].mC, filtered_buffer[i].uC,
                    get_context_string(context), tnc_str);
        }

        log_time("Finished writing text file: %s\n", txt_filename);
        fclose(txt_fp);
        txt_fp = NULL;
    }

    if (output_format == OUTPUT_HDF5 || output_format == OUTPUT_BOTH)
    {
        export_lock();
        hdf5_locked = 1;
        log_time("Starting to write HDF5 file: %s\n", filename);
        type = H5Tcreate(H5T_COMPOUND, sizeof(MethylRecord));
        H5Tinsert(type, "pos", HOFFSET(MethylRecord, pos), H5T_NATIVE_UINT32);
        H5Tinsert(type, "mC", HOFFSET(MethylRecord, mC), H5T_NATIVE_UINT16);
        H5Tinsert(type, "uC", HOFFSET(MethylRecord, uC), H5T_NATIVE_UINT16);
        H5Tinsert(type, "tnc", HOFFSET(MethylRecord, tnc), H5T_NATIVE_UINT8);

        file = append_mode
                   ? H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT)
                   : H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (file < 0)
        {
            fprintf(stderr, "Failed to %s HDF5 file: %s\n",
                    append_mode ? "open" : "create", filename);
            goto cleanup;
        }

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
            htri_t zstd_avail = H5Zfilter_avail(ZSTD_FILTER);
            int use_gzip = 1;
            if (zstd_avail > 0)
            {
                unsigned int cd_values[1] = {(unsigned int)compression};
                status = H5Pset_filter(dcpl, ZSTD_FILTER, H5Z_FLAG_MANDATORY, 1,
                                       cd_values);
                if (status < 0)
                    fprintf(stderr, "Warning: failed to set Zstd filter, falling "
                                    "back to gzip\n");
                else
                    use_gzip = 0;
            }
            else
                fprintf(stderr, "Note: Zstd HDF5 filter not available "
                                "(HDF5_PLUGIN_PATH not set?); using gzip\n");

            if (use_gzip)
            {
                int gzip_level = compression > 9 ? 9 : compression;
                if (H5Pset_deflate(dcpl, gzip_level) < 0)
                    fprintf(stderr, "Warning: gzip compression failed, writing "
                                    "uncompressed HDF5\n");
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
                H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count,
                                    NULL);
                status = H5Dwrite(dataset, mem_type, space, file_space, H5P_DEFAULT,
                                  filtered_buffer);
                H5Sclose(file_space);
            }
            else
            {
                dataset = H5Dcreate2(file, "methylation_data", type, space, H5P_DEFAULT,
                                     dcpl, H5P_DEFAULT);
                if (dataset >= 0)
                    status = H5Dwrite(dataset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                      filtered_buffer);
            }
        }
        else
        {
            dataset = H5Dcreate2(file, "methylation_data", type, space, H5P_DEFAULT,
                                 dcpl, H5P_DEFAULT);
            if (dataset >= 0)
                status = H5Dwrite(dataset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                  filtered_buffer);
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

        log_time("Finished writing HDF5 file: %s\n", filename);
    }

    MethylStats stats = calculate_statistics(filtered_buffer, dims[0]);
    stats.sites_in_reference = pre_stats.sites_in_reference;
    stats.sites_with_any_coverage = pre_stats.sites_with_any_coverage;
    stats.sites_below_min_cov = pre_stats.sites_below_min_cov;
    stats.sites_capped = sites_capped;
    stats.avg_coverage_all_sites = pre_stats.avg_coverage_all_sites;
    stats.coverage_p10 = pre_stats.coverage_p10;
    stats.coverage_median = pre_stats.coverage_median;
    stats.coverage_p90 = pre_stats.coverage_p90;
    stats.methylation_plus = pre_stats.methylation_plus;
    stats.methylation_minus = pre_stats.methylation_minus;
    stats.methylated_plus = pre_stats.methylated_plus;
    stats.unmethylated_plus = pre_stats.unmethylated_plus;
    stats.methylated_minus = pre_stats.methylated_minus;
    stats.unmethylated_minus = pre_stats.unmethylated_minus;

    if (meta)
        write_context_qc_json(filename, meta, filter_stats, stats);
    else
        write_context_qc_json(filename, NULL, filter_stats, stats);

    if (out_stats)
        *out_stats = stats;

    records_written = dims[0];
    log_time("Finished writing output files, wrote %llu filtered records\n",
             (unsigned long long)dims[0]);

cleanup:
    if (filtered_buffer)
        free(filtered_buffer);
    if (txt_fp)
        fclose(txt_fp);
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
    if (hdf5_locked)
        export_unlock();
    return records_written;
}

void cleanup_hdf5(void) { H5close(); }
