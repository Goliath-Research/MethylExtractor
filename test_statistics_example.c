#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cjson/cJSON.h"

// Mock structures for testing
typedef struct {
    unsigned tnc     : 5;
    unsigned context : 2;
    unsigned strand  : 1;
} tnc_bitfield_t;

typedef struct {
    uint32_t pos;
    uint16_t mC;
    uint16_t uC;
    tnc_bitfield_t tnc;
    uint8_t _pad[1];
} MethylRecord;

typedef struct {
    size_t num_positions;
    uint64_t total_methylated;
    uint64_t total_unmethylated;
    double avg_methylation_level;
    double avg_coverage;
} MethylStats;

// Function to calculate statistics from filtered buffer
static MethylStats calculate_statistics(MethylRecord *filtered_buffer, size_t n_records)
{
    MethylStats stats = {0, 0, 0, 0.0, 0.0};
    
    if (n_records == 0) {
        return stats;
    }
    
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

// Function to write statistics JSON file
static int write_statistics_json(const char *base_filename, MethylStats stats)
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
    {
        // Fallback: append .json to the base filename
        snprintf(json_filename, sizeof(json_filename), "%s.json", base_filename);
    }
    
    // Create JSON object
    cJSON *root = cJSON_CreateObject();
    if (!root) 
    {
        fprintf(stderr, "Failed to create JSON root object\n");
        return -1;
    }
    
    // Add statistics to JSON
    cJSON_AddNumberToObject(root, "num_positions", (double)stats.num_positions);
    cJSON_AddNumberToObject(root, "total_methylated", (double)stats.total_methylated);
    cJSON_AddNumberToObject(root, "total_unmethylated", (double)stats.total_unmethylated);
    cJSON_AddNumberToObject(root, "avg_methylation_level", stats.avg_methylation_level);
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
        fprintf(stderr, "Failed to open JSON file for writing: %s\n", json_filename);
        free(json_string);
        cJSON_Delete(root);
        return -1;
    }
    
    fprintf(json_fp, "%s\n", json_string);
    fclose(json_fp);
    
    printf("Wrote statistics to: %s\n", json_filename);
    printf("JSON content:\n%s\n", json_string);
    
    // Cleanup
    free(json_string);
    cJSON_Delete(root);
    
    return 0;
}

int main() 
{
    // Create sample methylation data
    MethylRecord sample_data[] = {
        {1000, 15, 5, {0}, {0}},   // 75% methylation, coverage 20
        {1001, 8, 12, {0}, {0}},   // 40% methylation, coverage 20
        {1002, 20, 0, {0}, {0}},   // 100% methylation, coverage 20
        {1003, 0, 25, {0}, {0}},   // 0% methylation, coverage 25
        {1004, 10, 10, {0}, {0}},  // 50% methylation, coverage 20
    };
    
    size_t n_records = sizeof(sample_data) / sizeof(sample_data[0]);
    
    printf("Sample methylation data:\n");
    for (size_t i = 0; i < n_records; i++) 
    {
        printf("Position %u: mC=%u, uC=%u, Coverage=%u, Methylation=%.1f%%\n",
               sample_data[i].pos,
               sample_data[i].mC,
               sample_data[i].uC,
               sample_data[i].mC + sample_data[i].uC,
               (double)sample_data[i].mC / (sample_data[i].mC + sample_data[i].uC) * 100.0);
    }
    
    // Calculate statistics
    MethylStats stats = calculate_statistics(sample_data, n_records);
    
    printf("\nCalculated statistics:\n");
    printf("Number of positions: %zu\n", stats.num_positions);
    printf("Total methylated: %llu\n", (unsigned long long)stats.total_methylated);
    printf("Total unmethylated: %llu\n", (unsigned long long)stats.total_unmethylated);
    printf("Average methylation level: %.4f (%.2f%%)\n", stats.avg_methylation_level, stats.avg_methylation_level * 100.0);
    printf("Average coverage: %.2f\n", stats.avg_coverage);
    
    // Write statistics JSON files for different output formats
    printf("\nGenerating JSON statistics files:\n");
    
    write_statistics_json("output.txt", stats);
    write_statistics_json("output.h5", stats);
    write_statistics_json("chromosome1-CG.h5", stats);
    
    return 0;
} 