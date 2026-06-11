#include "methyl_extractor.h"

MethylStats calculate_statistics(MethylRecord *filtered_buffer,
                                 size_t n_records)
{
    MethylStats stats = {0, 0, 0, 0.0, 0.0};

    if (n_records == 0)
        return stats;

    stats.num_positions = n_records;

    // Use uint64_t to avoid overflow for large datasets
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

int write_statistics_json(const char *base_filename, MethylStats stats)
{
    char json_filename[1024];

    // Create JSON filename by replacing the extension with .json
    const char *dot = strrchr(base_filename, '.');
    if (dot && (strcmp(dot, ".txt") == 0 || strcmp(dot, ".h5") == 0))
    {
        size_t base_len = dot - base_filename;
        strncpy(json_filename, base_filename, base_len);
        json_filename[base_len] = '\0';
        strcat(json_filename, ".json");
    }
    else
        snprintf(json_filename, sizeof(json_filename), "%s.json", base_filename);

    // Create JSON object
    cJSON *root = cJSON_CreateObject();
    if (!root)
    {
        fprintf(stderr, "Failed to create JSON root object\n");
        return -1;
    }

    // Add statistics to JSON
    cJSON_AddNumberToObject(root, "num_positions", (double)stats.num_positions);
    cJSON_AddNumberToObject(root, "total_methylated",
                            (double)stats.total_methylated);
    cJSON_AddNumberToObject(root, "total_unmethylated",
                            (double)stats.total_unmethylated);
    cJSON_AddNumberToObject(root, "avg_methylation_level",
                            stats.avg_methylation_level);
    cJSON_AddNumberToObject(root, "avg_coverage", stats.avg_coverage);

    // Write JSON to file
    char *json_string = cJSON_Print(root);
    if (!json_string)
    {
        fprintf(stderr, "Failed to print JSON\n");
        cJSON_Delete(root);
        return -1;
    }

    FILE *json_fp = fopen(json_filename, "w");
    if (!json_fp)
    {
        fprintf(stderr, "Failed to open JSON file for writing: %s\n",
                json_filename);
        free(json_string);
        cJSON_Delete(root);
        return -1;
    }

    fprintf(json_fp, "%s\n", json_string);
    fclose(json_fp);

    log_time("Wrote statistics to: %s\n", json_filename);

    // Cleanup
    free(json_string);
    cJSON_Delete(root);

    return 0;
}

size_t flush_buffer(const char *filename, MethylRecord *buffer,
                    size_t n_records, int compression, int chunk_size,
                    int append_mode, int min_cov, int cap_cov,
                    OutputFormat output_format)
{
    size_t records_written = 0;
    hsize_t dims[1] = {0};
    FILE *txt_fp = NULL;
    hid_t file = -1, dataset = -1, space = -1, type = -1, mem_type = -1,
          dcpl = -1;
    herr_t status = -1;

    // Calculate the number of records to write (filtered by coverage)
    // and the average coverage
    double avg_cov = 0.0;
    for (size_t i = 0; i < n_records; i++)
    {
        int total = buffer[i].mC + buffer[i].uC;
        if (total >= min_cov)
        {
            dims[0]++;
            // Online mean update using dims[0] as valid_count
            avg_cov += ((double)total - avg_cov) / dims[0];
        }
    }

    // If no records pass the filter, return early
    if (dims[0] == 0)
    {
        log_time("No records passed the coverage filter\n");
        return 0;
    }

    // Prepare filtered buffer
    MethylRecord *filtered_buffer = malloc(dims[0] * sizeof(MethylRecord));
    if (!filtered_buffer)
    {
        fprintf(stderr, "Failed to allocate memory for filtered_buffer\n");
        goto cleanup;
    }

    // fill filtered buffer
    size_t j = 0;
    for (size_t i = 0; i < n_records; i++)
    {
        int coverage = buffer[i].mC + buffer[i].uC;
        if (coverage >= min_cov)
        {
            // Only cap coverage if enabled
            if (cap_cov && coverage > avg_cov)
            {
                double prop = (double)buffer[i].mC / coverage;
                buffer[i].mC = (uint16_t)round((avg_cov * prop));
                buffer[i].uC = (uint16_t)(avg_cov - buffer[i].mC);
            }
            filtered_buffer[j++] = buffer[i];
        }
    }

    // Handle TXT output
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

    // Note: Parquet output removed to maintain pure C implementation
    // Use HDF5 format with proper compression instead

    // Handle HDF5 output
    if (output_format == OUTPUT_HDF5 || output_format == OUTPUT_BOTH)
    {
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

        // Set compression. Prefer Zstd, but only if the filter is actually
        // available at runtime; otherwise fall back to gzip deterministically.
        // (H5Pset_filter alone does not detect a missing plugin, and
        // H5Z_FLAG_OPTIONAL would silently write uncompressed data instead.)
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
                {
                    log_time("Using Zstd compression level %d\n", compression);
                    use_gzip = 0;
                }
            }
            else
                fprintf(stderr, "Note: Zstd HDF5 filter not available "
                                "(HDF5_PLUGIN_PATH not set?); using gzip\n");

            if (use_gzip)
            {
                // gzip/deflate supports levels 0-9
                int gzip_level = compression > 9 ? 9 : compression;
                if (H5Pset_deflate(dcpl, gzip_level) < 0)
                    fprintf(stderr, "Warning: gzip compression failed, writing "
                                    "uncompressed HDF5\n");
                else
                    log_time("Using gzip compression level %d\n", gzip_level);
            }
        }
        else
        {
            // compression == 0 → truly raw
            log_time("Writing raw uncompressed HDF5\n");
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

    // Calculate and write statistics (same for both formats since content is
    // identical)
    MethylStats stats = calculate_statistics(filtered_buffer, dims[0]);
    write_statistics_json(filename, stats);

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
    return records_written;
}

void cleanup_hdf5(void) { H5close(); }
