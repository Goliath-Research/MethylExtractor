#include "methyl_extractor.h"

static void utc_timestamp(char *buf, size_t buflen)
{
    time_t now = time(NULL);
    struct tm tm_utc;
    gmtime_r(&now, &tm_utc);
    strftime(buf, buflen, "%Y-%m-%dT%H:%M:%SZ", &tm_utc);
}

static void json_path_from_base(const char *base_filename, char *json_filename,
                                size_t json_filename_len)
{
    const char *dot = strrchr(base_filename, '.');
    if (dot && (strcmp(dot, ".txt") == 0 || strcmp(dot, ".h5") == 0))
    {
        size_t base_len = (size_t)(dot - base_filename);
        if (base_len >= json_filename_len)
            base_len = json_filename_len - 1;
        memcpy(json_filename, base_filename, base_len);
        json_filename[base_len] = '\0';
        strncat(json_filename, ".json", json_filename_len - base_len - 1);
    }
    else
        snprintf(json_filename, json_filename_len, "%s.json", base_filename);
}

void path_basename(const char *path, char *out, size_t out_len)
{
    if (!path || !out || out_len == 0)
        return;
    const char *base = strrchr(path, '/');
    base = base ? base + 1 : path;
    strncpy(out, base, out_len - 1);
    out[out_len - 1] = '\0';
    char *dot = strrchr(out, '.');
    if (dot)
        *dot = '\0';
}

static cJSON *filter_stats_to_json(const FilterStats *fs)
{
    cJSON *obj = cJSON_CreateObject();
    if (!obj || !fs)
        return obj;

    cJSON_AddNumberToObject(obj, "reads_seen", (double)fs->reads_seen);
    cJSON_AddNumberToObject(obj, "reads_used", (double)fs->reads_used);
    cJSON_AddNumberToObject(obj, "reads_dropped_unmapped",
                            (double)fs->reads_dropped_unmapped);
    cJSON_AddNumberToObject(obj, "reads_dropped_secondary",
                            (double)fs->reads_dropped_secondary);
    cJSON_AddNumberToObject(obj, "reads_dropped_qc_fail",
                            (double)fs->reads_dropped_qc_fail);
    cJSON_AddNumberToObject(obj, "reads_dropped_duplicate",
                            (double)fs->reads_dropped_duplicate);
    cJSON_AddNumberToObject(obj, "reads_dropped_supplementary",
                            (double)fs->reads_dropped_supplementary);
    cJSON_AddNumberToObject(obj, "reads_dropped_low_mapq",
                            (double)fs->reads_dropped_low_mapq);
    cJSON_AddNumberToObject(obj, "reads_dropped_multimap",
                            (double)fs->reads_dropped_multimap);
    cJSON_AddNumberToObject(obj, "reads_dropped_no_strand",
                            (double)fs->reads_dropped_no_strand);
    cJSON_AddNumberToObject(obj, "bases_skipped_overlap_clip",
                            (double)fs->bases_skipped_overlap_clip);
    cJSON_AddNumberToObject(obj, "bases_skipped_low_phred",
                            (double)fs->bases_skipped_low_phred);
    cJSON_AddNumberToObject(obj, "bases_skipped_non_cytosine",
                            (double)fs->bases_skipped_non_cytosine);
    cJSON_AddNumberToObject(obj, "bases_counted_methylated",
                            (double)fs->bases_counted_methylated);
    cJSON_AddNumberToObject(obj, "bases_counted_unmethylated",
                            (double)fs->bases_counted_unmethylated);

    if (fs->reads_seen > 0)
    {
        cJSON_AddNumberToObject(obj, "read_retention_rate",
                                (double)fs->reads_used / (double)fs->reads_seen);
    }

    return obj;
}

static cJSON *methyl_stats_to_json(const MethylStats *stats)
{
    cJSON *obj = cJSON_CreateObject();
    if (!obj || !stats)
        return obj;

    cJSON_AddNumberToObject(obj, "num_positions", (double)stats->num_positions);
    cJSON_AddNumberToObject(obj, "methylation_level", stats->avg_methylation_level);
    cJSON_AddNumberToObject(obj, "mean_coverage", stats->avg_coverage);
    cJSON_AddNumberToObject(obj, "total_methylated", (double)stats->total_methylated);
    cJSON_AddNumberToObject(obj, "total_unmethylated",
                            (double)stats->total_unmethylated);
    return obj;
}

int write_context_qc_json(const char *base_filename, const ExtractionMeta *meta,
                          const FilterStats *filter_stats, MethylStats stats)
{
    uint64_t t0 = now_ms();
    export_lock();
    char json_filename[1024];
    json_path_from_base(base_filename, json_filename, sizeof(json_filename));

    cJSON *root = cJSON_CreateObject();
    if (!root)
    {
        fprintf(stderr, "Failed to create JSON root object\n");
        export_unlock();
        return -1;
    }

    // Legacy top-level fields for backward compatibility.
    cJSON_AddNumberToObject(root, "num_positions", (double)stats.num_positions);
    cJSON_AddNumberToObject(root, "total_methylated",
                            (double)stats.total_methylated);
    cJSON_AddNumberToObject(root, "total_unmethylated",
                            (double)stats.total_unmethylated);
    cJSON_AddNumberToObject(root, "avg_methylation_level",
                            stats.avg_methylation_level);
    cJSON_AddNumberToObject(root, "avg_coverage", stats.avg_coverage);

    cJSON *metadata = cJSON_CreateObject();
    char ts[32];
    utc_timestamp(ts, sizeof(ts));
    cJSON_AddStringToObject(metadata, "schema_name", EXTRACTION_CONTEXT_QC_SCHEMA);
    cJSON_AddStringToObject(metadata, "schema_version", EXTRACTION_CONTEXT_QC_VERSION);
    cJSON_AddStringToObject(metadata, "exported_at_utc", ts);
    if (meta && meta->chromosome)
        cJSON_AddStringToObject(metadata, "chromosome", meta->chromosome);
    if (meta && meta->context)
        cJSON_AddStringToObject(metadata, "context", meta->context);
    if (meta && meta->output_path)
        cJSON_AddStringToObject(metadata, "companion_output", meta->output_path);

    cJSON *filters = cJSON_CreateObject();
    if (meta)
    {
        cJSON_AddNumberToObject(filters, "min_mapq", meta->min_mapq);
        cJSON_AddNumberToObject(filters, "min_phred", meta->min_phred);
        cJSON_AddNumberToObject(filters, "min_cov", meta->min_cov);
        cJSON_AddNumberToObject(filters, "cap_cov", meta->cap_cov);
    }
    cJSON_AddItemToObject(metadata, "filters", filters);
    cJSON_AddItemToObject(root, "metadata", metadata);

    cJSON *sites = cJSON_CreateObject();
    cJSON_AddNumberToObject(sites, "sites_in_reference",
                            (double)stats.sites_in_reference);
    cJSON_AddNumberToObject(sites, "sites_passing_min_cov",
                            (double)stats.num_positions);
    cJSON_AddNumberToObject(sites, "sites_with_any_coverage",
                            (double)stats.sites_with_any_coverage);
    cJSON_AddNumberToObject(sites, "sites_below_min_cov",
                            (double)stats.sites_below_min_cov);
    cJSON_AddNumberToObject(sites, "sites_capped", (double)stats.sites_capped);
    if (stats.sites_in_reference > 0)
    {
        cJSON_AddNumberToObject(sites, "fraction_sites_covered",
                                (double)stats.num_positions /
                                    (double)stats.sites_in_reference);
    }
    cJSON_AddItemToObject(root, "sites", sites);

    cJSON *coverage = cJSON_CreateObject();
    cJSON_AddNumberToObject(coverage, "mean_on_passing_sites", stats.avg_coverage);
    cJSON_AddNumberToObject(coverage, "mean_all_sites", stats.avg_coverage_all_sites);
    cJSON_AddNumberToObject(coverage, "median_all_sites", stats.coverage_median);
    cJSON_AddNumberToObject(coverage, "p10_all_sites", stats.coverage_p10);
    cJSON_AddNumberToObject(coverage, "p90_all_sites", stats.coverage_p90);
    cJSON_AddItemToObject(root, "coverage", coverage);

    cJSON *methylation = cJSON_CreateObject();
    cJSON_AddNumberToObject(methylation, "level", stats.avg_methylation_level);
    cJSON_AddNumberToObject(methylation, "total_methylated",
                            (double)stats.total_methylated);
    cJSON_AddNumberToObject(methylation, "total_unmethylated",
                            (double)stats.total_unmethylated);
    cJSON *by_strand = cJSON_CreateObject();
    cJSON_AddNumberToObject(by_strand, "plus", stats.methylation_plus);
    cJSON_AddNumberToObject(by_strand, "minus", stats.methylation_minus);
    cJSON_AddItemToObject(methylation, "by_strand", by_strand);
    cJSON_AddItemToObject(root, "methylation", methylation);

    if (filter_stats)
    {
        cJSON *rf = filter_stats_to_json(filter_stats);
        cJSON_AddItemToObject(root, "read_filtering", rf);
    }

    char *json_string = cJSON_Print(root);
    if (!json_string)
    {
        fprintf(stderr, "Failed to print JSON\n");
        cJSON_Delete(root);
        export_unlock();
        return -1;
    }

    FILE *json_fp = fopen(json_filename, "w");
    if (!json_fp)
    {
        fprintf(stderr, "Failed to open JSON file for writing: %s\n",
                json_filename);
        free(json_string);
        cJSON_Delete(root);
        export_unlock();
        return -1;
    }

    fprintf(json_fp, "%s\n", json_string);
    fclose(json_fp);

    log_time("Wrote context QC JSON to: %s\n", json_filename);

    free(json_string);
    cJSON_Delete(root);
    if (meta && meta->timing)
        meta->timing->qc_json_ms += elapsed_ms(t0);
    export_unlock();
    return 0;
}

static int context_index(const char *context)
{
    if (!context)
        return -1;
    if (strcmp(context, "CG") == 0)
        return 0;
    if (strcmp(context, "CHG") == 0)
        return 1;
    if (strcmp(context, "CHH") == 0)
        return 2;
    return -1;
}

static void aggregate_context(uint64_t *total_m, uint64_t *total_u,
                              uint64_t *sites_passing, uint64_t *sites_reference,
                              const ContextExport *ctx)
{
    int idx = context_index(ctx->context);
    if (idx < 0)
        return;
    total_m[idx] += ctx->stats.total_methylated;
    total_u[idx] += ctx->stats.total_unmethylated;
    sites_passing[idx] += ctx->stats.num_positions;
    sites_reference[idx] += ctx->stats.sites_in_reference;
}

static cJSON *build_summary(const ChromosomeExport *chromosomes, int n_chromosomes,
                            int keep_chg, int keep_chh)
{
    uint64_t total_m[3] = {0};
    uint64_t total_u[3] = {0};
    uint64_t sites_passing[3] = {0};
    uint64_t sites_reference[3] = {0};

    for (int i = 0; i < n_chromosomes; i++)
    {
        for (int j = 0; j < chromosomes[i].n_contexts; j++)
            aggregate_context(total_m, total_u, sites_passing, sites_reference,
                              &chromosomes[i].contexts[j]);
    }

    cJSON *summary = cJSON_CreateObject();
    cJSON_AddNumberToObject(summary, "chromosomes_processed", n_chromosomes);

    if (sites_passing[0] > 0)
    {
        cJSON_AddNumberToObject(summary, "cpg_sites_passing_min_cov",
                                (double)sites_passing[0]);
        cJSON_AddNumberToObject(summary, "cpg_sites_in_reference",
                                (double)sites_reference[0]);
        uint64_t cov_total = total_m[0] + total_u[0];
        cJSON_AddNumberToObject(summary, "cpg_weighted_mean_coverage",
                                (double)cov_total / (double)sites_passing[0]);
        if (total_m[0] + total_u[0] > 0)
            cJSON_AddNumberToObject(summary, "cpg_methylation_level",
                                    (double)total_m[0] / (double)(total_m[0] + total_u[0]));
        if (sites_reference[0] > 0)
            cJSON_AddNumberToObject(summary, "cpg_fraction_sites_covered",
                                    (double)sites_passing[0] /
                                        (double)sites_reference[0]);
    }

    if (keep_chg && total_m[1] + total_u[1] > 0)
        cJSON_AddNumberToObject(summary, "chg_methylation_level",
                                (double)total_m[1] / (double)(total_m[1] + total_u[1]));
    if (keep_chh && total_m[2] + total_u[2] > 0)
        cJSON_AddNumberToObject(summary, "chh_methylation_level",
                                (double)total_m[2] / (double)(total_m[2] + total_u[2]));

    return summary;
}

int write_extraction_manifest(const SampleRunInfo *run,
                              const ChromosomeExport *chromosomes, int n_chromosomes)
{
    if (!run || !run->out_dir)
        return -1;

    char sample_id[256];
    path_basename(run->out_dir, sample_id, sizeof(sample_id));

    char manifest_path[1200];
    snprintf(manifest_path, sizeof(manifest_path), "%s/%s.extraction_manifest.json",
             run->out_dir, sample_id);

    cJSON *root = cJSON_CreateObject();
    if (!root)
        return -1;

    cJSON *metadata = cJSON_CreateObject();
    char ts[32];
    utc_timestamp(ts, sizeof(ts));
    cJSON_AddStringToObject(metadata, "schema_name", EXTRACTION_MANIFEST_SCHEMA);
    cJSON_AddStringToObject(metadata, "schema_version", EXTRACTION_MANIFEST_VERSION);
    cJSON_AddStringToObject(metadata, "exported_at_utc", ts);
    cJSON_AddStringToObject(metadata, "sample_id", sample_id);
    cJSON_AddStringToObject(metadata, "sample_dir", run->out_dir);
    if (run->bam_file)
        cJSON_AddStringToObject(metadata, "bam_file", run->bam_file);
    if (run->reference)
        cJSON_AddStringToObject(metadata, "reference", run->reference);

    cJSON *contexts = cJSON_CreateArray();
    cJSON_AddItemToArray(contexts, cJSON_CreateString("CG"));
    if (run->keep_chg)
        cJSON_AddItemToArray(contexts, cJSON_CreateString("CHG"));
    if (run->keep_chh)
        cJSON_AddItemToArray(contexts, cJSON_CreateString("CHH"));
    cJSON_AddItemToObject(metadata, "contexts_extracted", contexts);

    cJSON *filters = cJSON_CreateObject();
    cJSON_AddNumberToObject(filters, "min_mapq", run->min_mapq);
    cJSON_AddNumberToObject(filters, "min_phred", run->min_phred);
    cJSON_AddNumberToObject(filters, "min_cov", run->min_cov);
    cJSON_AddNumberToObject(filters, "cap_cov", run->cap_cov);
    cJSON_AddNumberToObject(filters, "split_context_files", run->split_context_files);
    cJSON_AddItemToObject(metadata, "filters", filters);
    cJSON_AddItemToObject(root, "metadata", metadata);

    cJSON_AddItemToObject(root, "summary",
                          build_summary(chromosomes, n_chromosomes, run->keep_chg,
                                        run->keep_chh));

    FilterStats total_filters;
    memset(&total_filters, 0, sizeof(total_filters));
    for (int i = 0; i < n_chromosomes; i++)
        merge_filter_stats(&total_filters, &chromosomes[i].filter_stats);
    cJSON_AddItemToObject(root, "read_filtering",
                          filter_stats_to_json(&total_filters));

    cJSON *per_chromosome = cJSON_CreateObject();
    for (int i = 0; i < n_chromosomes; i++)
    {
        const ChromosomeExport *chr = &chromosomes[i];
        cJSON *chr_obj = cJSON_CreateObject();
        for (int j = 0; j < chr->n_contexts; j++)
        {
            const ContextExport *ctx = &chr->contexts[j];
            cJSON *ctx_obj = methyl_stats_to_json(&ctx->stats);
            cJSON_AddStringToObject(ctx_obj, "companion_output", ctx->output_path);
            cJSON_AddItemToObject(chr_obj, ctx->context, ctx_obj);
        }
        cJSON_AddItemToObject(chr_obj, "read_filtering",
                              filter_stats_to_json(&chr->filter_stats));
        cJSON_AddItemToObject(per_chromosome, chr->chromosome, chr_obj);
    }
    cJSON_AddItemToObject(root, "per_chromosome", per_chromosome);

    char *json_string = cJSON_Print(root);
    if (!json_string)
    {
        cJSON_Delete(root);
        return -1;
    }

    FILE *fp = fopen(manifest_path, "w");
    if (!fp)
    {
        fprintf(stderr, "Failed to open extraction manifest: %s\n", manifest_path);
        free(json_string);
        cJSON_Delete(root);
        return -1;
    }

    fprintf(fp, "%s\n", json_string);
    fclose(fp);
    log_time("Wrote extraction manifest to: %s\n", manifest_path);

    free(json_string);
    cJSON_Delete(root);
    return 0;
}

int write_extraction_timing(const char *out_dir, const SampleTiming *sample,
                            const ChromosomeTiming *chromosomes, int n_chromosomes)
{
    if (!out_dir || !sample)
        return -1;

    char sample_id[256];
    path_basename(out_dir, sample_id, sizeof(sample_id));

    char path[1200];
    snprintf(path, sizeof(path), "%s/%s.timing.json", out_dir, sample_id);

    cJSON *root = cJSON_CreateObject();
    if (!root)
        return -1;

    cJSON_AddStringToObject(root, "schema_name", EXTRACTION_TIMING_SCHEMA);
    cJSON_AddStringToObject(root, "schema_version", EXTRACTION_TIMING_VERSION);
    cJSON_AddNumberToObject(root, "total_ms", (double)sample->total_ms);
    cJSON_AddNumberToObject(root, "bam_ms", (double)sample->bam_ms);
    cJSON_AddNumberToObject(root, "write_ms", (double)sample->write_ms);
    double tot = sample->total_ms > 0 ? (double)sample->total_ms : 1.0;
    cJSON_AddNumberToObject(root, "bam_fraction", (double)sample->bam_ms / tot);
    cJSON_AddNumberToObject(root, "write_fraction", (double)sample->write_ms / tot);
    cJSON_AddNumberToObject(root, "chrom_parallel", sample->chrom_parallel);
    cJSON_AddNumberToObject(root, "region_threads", sample->region_threads);
    cJSON_AddNumberToObject(root, "max_rss_gb", sample->max_rss_gb);
    cJSON_AddNumberToObject(root, "n_chromosomes", sample->n_chromosomes);

    cJSON *arr = cJSON_CreateArray();
    for (int i = 0; i < n_chromosomes; i++)
    {
        const ChromosomeTiming *ct = &chromosomes[i];
        cJSON *obj = cJSON_CreateObject();
        cJSON_AddStringToObject(obj, "name", ct->chromosome);
        cJSON_AddNumberToObject(obj, "site_enum_ms", (double)ct->site_enum_ms);
        cJSON_AddNumberToObject(obj, "bam_scan_ms", (double)ct->bam_scan_ms);
        cJSON_AddNumberToObject(obj, "merge_ms", (double)ct->merge_ms);
        cJSON_AddNumberToObject(obj, "hdf5_write_ms", (double)ct->hdf5_write_ms);
        cJSON_AddNumberToObject(obj, "qc_json_ms", (double)ct->qc_json_ms);
        cJSON_AddNumberToObject(obj, "total_ms", (double)ct->total_ms);
        cJSON_AddNumberToObject(obj, "rss_est_bytes", (double)ct->rss_est_bytes);
        cJSON_AddItemToArray(arr, obj);
    }
    cJSON_AddItemToObject(root, "chromosomes", arr);

    char *json_string = cJSON_Print(root);
    if (!json_string)
    {
        cJSON_Delete(root);
        return -1;
    }

    FILE *fp = fopen(path, "w");
    if (!fp)
    {
        fprintf(stderr, "Failed to open timing JSON: %s\n", path);
        free(json_string);
        cJSON_Delete(root);
        return -1;
    }
    fprintf(fp, "%s\n", json_string);
    fclose(fp);
    log_time("Wrote extraction timing to: %s\n", path);
    free(json_string);
    cJSON_Delete(root);
    return 0;
}
