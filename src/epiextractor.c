#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#include <time.h> // For clock_gettime, CLOCK_REALTIME, and CLOCK_MONOTONIC
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <errno.h>
#include <limits.h>
#include <pthread.h>
#include <hdf5/serial/hdf5.h>
#include <hdf5/serial/hdf5_hl.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/stat.h>
#include <ctype.h>
#include <sys/sysinfo.h> // For getting system information like available memory and CPU count

#define BUFFER_SIZE 100000
#define MAX(a, b) ((a) > (b) ? (a) : (b))

// Forward declarations of inline functions
static inline int isCpG(char *seq, hts_pos_t pos, int seqlen);
static inline int isCHG(char *seq, hts_pos_t pos, int seqlen);
static inline int isCHH(char *seq, hts_pos_t pos, int seqlen);

// Data structures
typedef struct 
{
    char *FastaName;
    char *BAMName;
    char *hdf5_output_dir;
    int hdf5_compression_level;
    int hdf5_chunk_size;
    int nThreads;           // Now will be calculated based on system resources
    size_t chunkSize;       // Adjusted based on memory availability
    size_t bufferSizePerThread; // Buffer size per thread for in-memory processing
    int minMapq;
    int minPhred;
    uint32_t minDepth;
    uint32_t maxDepth;
    int keepCpG;
    int keepCHG;
    int keepCHH;
    hts_idx_t *bai;    // BAM index
    faidx_t *fai;      // FASTA index
    bam_hdr_t *hdr;    // BAM header
} Config;

// Add a new structure for accessing bit fields of flags and tnc
typedef struct {
    unsigned int context : 2;  // Bits 0-1 for methylation context
    unsigned int strand : 1;   // Bit 2 for strand
    unsigned int chrom : 5;    // Bits 3-7 for chromosome number
} FlagsBitField;

typedef struct {
    unsigned int base3 : 2;    // Bits 0-1 for third base
    unsigned int base2 : 2;    // Bits 2-3 for second base
    unsigned int base1 : 2;    // Bits 4-5 for first base
    unsigned int reserved : 2; // Bits 6-7 reserved
} TncBitField;

// MethylData will keep storage as uint8_t for HDF5 compatibility
typedef struct 
{
    uint8_t  flags;    // Stored as uint8_t for HDF5
    uint32_t pos;
    uint32_t nmethyl;
    uint32_t nunmethyl;
    uint8_t  tnc;      // Stored as uint8_t for HDF5
} MethylData;

// Add inline functions to set and get values using bit fields
static inline void setMethylFlags(MethylData* mdata, uint8_t chrom, uint8_t strand, uint8_t context) 
{
    FlagsBitField* fbf = (FlagsBitField*)&mdata->flags;
    fbf->chrom = chrom;
    fbf->strand = strand;
    fbf->context = context;
}

static inline uint8_t getChrom(MethylData* mdata) 
{
    FlagsBitField* fbf = (FlagsBitField*)&mdata->flags;
    return fbf->chrom;
}

/* static inline uint8_t _getStrand(MethylData* mdata) 
{
    FlagsBitField* fbf = (FlagsBitField*)&mdata->flags;
    return fbf->strand;
} */

static inline uint8_t getContext(MethylData* mdata) 
{
    FlagsBitField* fbf = (FlagsBitField*)&mdata->flags;
    return fbf->context;
}

static inline void setTnc(MethylData* mdata, uint8_t base1, uint8_t base2, uint8_t base3) 
{
    TncBitField* tbf = (TncBitField*)&mdata->tnc;
    tbf->base1 = base1;
    tbf->base2 = base2;
    tbf->base3 = base3;
}

static inline uint8_t getTncBase1(MethylData* mdata) 
{
    TncBitField* tbf = (TncBitField*)&mdata->tnc;
    return tbf->base1;
}

static inline uint8_t getTncBase2(MethylData* mdata) 
{
    TncBitField* tbf = (TncBitField*)&mdata->tnc;
    return tbf->base2;
}

static inline uint8_t getTncBase3(MethylData* mdata) 
{
    TncBitField* tbf = (TncBitField*)&mdata->tnc;
    return tbf->base3;
}

static int getStrand(const bam1_t *b) 
{
    char *XG = (char *) bam_aux_get(b, "XG");
    //Only bismark uses the XG tag like this. Some other aligners use it for other purposes...
    if(XG != NULL && *(XG+1) != 'C' && *(XG+1) != 'G') XG = NULL;
    if(XG == NULL) { //Can't handle non-directional libraries!
        if(b->core.flag & BAM_FPAIRED) {
            if((b->core.flag & 0x50) == 0x50) return 2; //Read1, reverse comp. == OB
            else if(b->core.flag & 0x40) return 1; //Read1, forward == OT
            else if((b->core.flag & 0x90) == 0x90) return 1; //Read2, reverse comp. == OT
            else if(b->core.flag & 0x80) return 2; //Read2, forward == OB
            return 0; //One of the above should be set!
        } else {
            if(b->core.flag & 0x10) return 2; //Reverse comp. == OB
            return 1; //OT
        }
    } else {
        if(*(XG+1) == 'C') { //OT or CTOT, due to C->T converted genome
            if((b->core.flag & 0x51) == 0x41) return 1; //Read#1 forward == OT
            else if((b->core.flag & 0x51) == 0x51) return 3; //Read #1 reverse == CTOT
            else if((b->core.flag & 0x91) == 0x81) return 3; //Read #2 forward == CTOT
            else if((b->core.flag & 0x91) == 0x91) return 1; //Read #2 reverse == OT
            else if(b->core.flag & 0x10) return 3; //Single-end reverse == CTOT
            else return 1; //Single-end forward == OT
        } else {
            if((b->core.flag & 0x51) == 0x41) return 4; //Read#1 forward == CTOB
            else if((b->core.flag & 0x51) == 0x51) return 2; //Read #1 reverse == OB
            else if((b->core.flag & 0x91) == 0x81) return 2; //Read #2 forward == OB
            else if((b->core.flag & 0x91) == 0x91) return 4; //Read #2 reverse == CTOB
            else if(b->core.flag & 0x10) return 2; //Single-end reverse == OB
            else return 4; //Single-end forward == CTOB
        }
    }
}

static int updateMetrics(Config *config, const bam_pileup1_t *plp) 
{
    uint8_t base = bam_seqi(bam_get_seq(plp->b), plp->qpos);
    int strand = getStrand(plp->b); //1=OT, 2=OB, 3=CTOT, 4=CTOB

    if(strand==0) return 0;
    //Is the phred score even high enough?
    if(bam_get_qual(plp->b)[plp->qpos] < config->minPhred) return 0;

    if(base == 2 && (strand==1 || strand==3)) return 1; //C on an OT/CTOT alignment
    else if(base == 8 && (strand==1 || strand==3)) return -1; //T on an OT/CTOT alignment
    else if(base == 4 && (strand==2 || strand==4)) return 1; //G on an OB/CTOB alignment
    else if(base == 1 && (strand==2 || strand==4)) return -1; //A on an OB/CTOB alignment
    return 0;
}

typedef struct {
    Config *config;
    bam_hdr_t *hdr;
    htsFile *fp;
    hts_itr_t *iter;
    char *seq;
    int32_t offset;
    int lseq;
    void *ohash;
} mplp_data;

// Forward declarations of static functions
static hid_t create_hdf5_file(const char* chrom, const char* output_dir, Config* config);
static int write_to_hdf5(hid_t handles, MethylData* buffer, size_t count);
static void close_hdf5_file(hid_t handles);
void *extractCalls(void *foo);

// Add these function declarations at the top with other declarations
static int plp_get_read(void *data, bam1_t *b);
static void plp_set_region(void *data, int tid, hts_pos_t start, hts_pos_t end);

// Add this structure definition with other structs
typedef struct {
    samFile *fp;
    hts_itr_t *iter;
    int tid;
    hts_pos_t start, end;
} aux_t;

// Add reference counting to HDF5 handles
typedef struct {
    hid_t file_id;
    hid_t dataset;
    hid_t compound_type;
    int ref_count;  // Add reference counting
    pthread_mutex_t lock;  // Per-handle lock
} HDF5Handles;

// Add global mutex for handle management
static pthread_mutex_t handle_mutex = PTHREAD_MUTEX_INITIALIZER;

// Modify ChromTransition to include error tracking
typedef struct {
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    int active_threads;         // Total active threads in the program
    char current_chrom[32];     // Current chromosome being processed
    int transition_pending;     // Flag for pending transition
    int threads_in_transition;  // Threads waiting to transition
    int working_threads;        // Threads actively working on current chromosome
    int error_count;            // Track errors across threads
} ChromTransition;

static ChromTransition chrom_transition = {
    .mutex = PTHREAD_MUTEX_INITIALIZER,
    .cond = PTHREAD_COND_INITIALIZER,
    .active_threads = 0,
    .current_chrom = "",
    .transition_pending = 0,
    .threads_in_transition = 0,
    .working_threads = 0,    // Initialize to 0
    .error_count = 0
};

// Implementation of inline functions
static inline int isCpG(char *seq, hts_pos_t pos, int seqlen) 
{
    if (pos >= seqlen) return 0;
    if (*(seq+pos) == 'C' || *(seq+pos) == 'c') 
    {
        if (pos+1 >= seqlen) return 0;
        if (*(seq+pos+1) == 'G' || *(seq+pos+1) == 'g') return 1;
        return 0;
    } 
    else if (*(seq+pos) == 'G' || *(seq+pos) == 'g') 
    {
        if (pos == 0) return 0;
        if (*(seq+pos-1) == 'C' || *(seq+pos-1) == 'c') return -1;
        return 0;
    }
    return 0;
}

static inline int isCHG(char *seq, hts_pos_t pos, int seqlen) 
{
    if (pos >= seqlen) return 0;
    if (*(seq+pos) == 'C' || *(seq+pos) == 'c') 
    {
        if (pos+2 >= seqlen) return 0;
        if (*(seq+pos+2) == 'G' || *(seq+pos+2) == 'g') return 1;
        return 0;
    } 
    else if (*(seq+pos) == 'G' || *(seq+pos) == 'g') 
    {
        if (pos <= 1) return 0;
        if (*(seq+pos-2) == 'C' || *(seq+pos-2) == 'c') return -1;
        return 0;
    }
    return 0;
}

static inline int isCHH(char *seq, hts_pos_t pos, int seqlen) 
{
    if (pos >= seqlen) return 0;
    if (*(seq+pos) == 'C' || *(seq+pos) == 'c') return 1;
    else if (*(seq+pos) == 'G' || *(seq+pos) == 'g') return -1;
    return 0;
}

static inline uint8_t encodeStrand(char strand) 
{
    return (strand == '+' || strand == 'F') ? 1 : 0;
}

static inline uint8_t encodeContext(int type) 
{
    return (uint8_t)type;  // Directly return type as it maps to CpG=0, CHG=1, CHH=2
}

static inline uint8_t chromToNumber(const char* chrom) 
{
    if (!chrom) return 0;
    if (strcmp(chrom, "X") == 0 || strcmp(chrom, "chrX") == 0) return 23;
    if (strcmp(chrom, "Y") == 0 || strcmp(chrom, "chrY") == 0) return 24;
    if (strncmp(chrom, "chr", 3) == 0) chrom += 3;
    char* end;
    long num = strtol(chrom, &end, 10);
    if (*end == '\0' && num > 0 && num <= 24) return num;
    return 0;
}

static inline void chromName(const char* chrom, char* chrom_name, size_t chrom_name_size) 
{
    if (!chrom || !chrom_name || chrom_name_size == 0) 
    {
        if (chrom_name && chrom_name_size > 0) chrom_name[0] = '\0';
        return;
    }
    if (strcmp(chrom, "X") == 0 || strcmp(chrom, "chrX") == 0) 
    {
        strncpy(chrom_name, "X", chrom_name_size - 1);
        chrom_name[chrom_name_size - 1] = '\0';
        return;
    }
    if (strcmp(chrom, "Y") == 0 || strcmp(chrom, "chrY") == 0) 
    {
        strncpy(chrom_name, "Y", chrom_name_size - 1);
        chrom_name[chrom_name_size - 1] = '\0';
        return;
    }
    if (strncmp(chrom, "chr", 3) == 0) chrom += 3;
    strncpy(chrom_name, chrom, chrom_name_size - 1);
    chrom_name[chrom_name_size - 1] = '\0';
    char* end;
    long num = strtol(chrom, &end, 10);
    if (*end != '\0' || num <= 0 || num > 24) 
    {
        // If not a valid number or out of range, keep the name as is
        return;
    }
    // If it's a valid number, we've already copied it
    return;
}

static char revcomp(char b) 
{
    switch(b) 
    {
        case 'A': case 'a': return 'T';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        case 'T': case 't': return 'A';
        default: return 'N';
    }
}

static inline uint8_t encodeTrinucleotide(const char* seq, hts_pos_t pos, int seqlen, int direction) 
{
    if (pos < 0 || pos >= seqlen) return 0;
    char base1, base2, base3;
    if (direction > 0) 
    {
        base1 = (pos > 0) ? toupper(seq[pos - 1]) : 'N';
        base2 = toupper(seq[pos]);
        base3 = (pos + 1 < seqlen) ? toupper(seq[pos + 1]) : 'N';
    } 
    else 
    {
        base1 = (pos > 0) ? toupper(seq[pos - 1]) : 'N';
        base2 = toupper(seq[pos]);
        base3 = (pos + 1 < seqlen) ? toupper(seq[pos + 1]) : 'N';
        base1 = (base1 == 'N') ? 'N' : revcomp(base1);
        base2 = (base2 == 'N') ? 'N' : revcomp(base2);
        base3 = (base3 == 'N') ? 'N' : revcomp(base3);
    }
    uint8_t b1 = (base1 == 'C') ? 1 : (base1 == 'G') ? 2 : (base1 == 'T') ? 3 : 0;
    uint8_t b2 = (base2 == 'C') ? 1 : (base2 == 'G') ? 2 : (base2 == 'T') ? 3 : 0;
    uint8_t b3 = (base3 == 'C') ? 1 : (base3 == 'G') ? 2 : (base3 == 'T') ? 3 : 0;
    MethylData temp;
    setTnc(&temp, b1, b2, b3);
    return temp.tnc;
}

static void init_hdf5_handles(HDF5Handles* h) 
{
    h->file_id = -1;
    h->dataset = -1;
    h->compound_type = -1;
    h->ref_count = 1;
    pthread_mutex_init(&h->lock, NULL);
}

static void unref_hdf5_handles(HDF5Handles* h) 
{
    if (!h) return;
    
    pthread_mutex_lock(&handle_mutex);  // Global lock for reference counting
    pthread_mutex_lock(&h->lock);       // Per-handle lock
    
    h->ref_count--;
    int should_free = (h->ref_count == 0);
    
    if (should_free) 
    {
        // Close all HDF5 resources under the handle's own lock
        if (h->dataset >= 0) 
        {
            H5Dclose(h->dataset);
            h->dataset = -1;
        }
        if (h->compound_type >= 0) 
        {
            H5Tclose(h->compound_type);
            h->compound_type = -1;
        }
        if (h->file_id >= 0) 
        {
            // Try to flush one last time
            H5Fflush(h->file_id, H5F_SCOPE_GLOBAL);
            H5Fclose(h->file_id);
            h->file_id = -1;
        }
    }
    
    pthread_mutex_unlock(&h->lock);
    
    if (should_free) 
    {
        pthread_mutex_destroy(&h->lock);
        free(h);
    }
    
    pthread_mutex_unlock(&handle_mutex);
}

static hid_t create_hdf5_file(const char* chrom, const char* output_dir, Config* config) 
{
    if (!chrom || chrom[0] == '\0') 
    {
        return -1;
    }
    char filepath[512];
    char chrom_name[32];
    // Convert chromosome name to number or letter without 'chr' prefix using chromName function
    chromName(chrom, chrom_name, sizeof(chrom_name));
    snprintf(filepath, sizeof(filepath), "%s/%s.h5", output_dir, chrom_name);
    // Check if file already exists and delete it to overwrite
    if (access(filepath, F_OK) == 0) 
    {
        if (unlink(filepath) != 0) 
        {
            return -1;
        }
    }
    char filename[2048];
    // Use the provided chrom as the full filename if it contains path separators, otherwise construct path
    if (strchr(chrom, '/') != NULL || strchr(chrom, '\\') != NULL) 
    {
        strncpy(filename, chrom, sizeof(filename) - 1);
        filename[sizeof(filename) - 1] = '\0';
    } 
    else 
    {
        int written = snprintf(filename, sizeof(filename), "%s/%s.h5", output_dir, chrom_name);
        if (written < 0 || (size_t)written >= sizeof(filename)) 
        {
            //fprintf(stderr, "Path too long for HDF5 file: %s/%s.h5\n", output_dir, chrom);
            return -1;
        }
    }
    
    // Create file access property list with thread-safe settings
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    if (fapl < 0) 
    {
        //fprintf(stderr, "Failed to create file access property list\n");
        return -1;
    }
    
    // Set thread-safe file access properties
    H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);
    H5Pset_cache(fapl, 0, 521, 8*1024*1024, 0.3);
    
    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    H5Pclose(fapl);
    
    if (file_id < 0) 
    {
        //fprintf(stderr, "Failed to create HDF5 file at %s\n", filename);
        return -1;
    }
    
    // Create compound datatype first
    hid_t compound_type = H5Tcreate(H5T_COMPOUND, sizeof(MethylData));
    if (compound_type < 0) 
    {
        //fprintf(stderr, "Failed to create compound datatype\n");
        H5Fclose(file_id);
        return -1;
    }
    
    H5Tinsert(compound_type, "flags", HOFFSET(MethylData, flags), H5T_NATIVE_UINT8);
    H5Tinsert(compound_type, "pos", HOFFSET(MethylData, pos), H5T_NATIVE_UINT32);
    H5Tinsert(compound_type, "nmethyl", HOFFSET(MethylData, nmethyl), H5T_NATIVE_UINT32);
    H5Tinsert(compound_type, "nunmethyl", HOFFSET(MethylData, nunmethyl), H5T_NATIVE_UINT32);
    H5Tinsert(compound_type, "tnc", HOFFSET(MethylData, tnc), H5T_NATIVE_UINT8);
    
    // Create dataset creation property list
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    if (dcpl < 0) 
    {
        //fprintf(stderr, "Failed to create dataset creation property list\n");
        H5Tclose(compound_type);
        H5Fclose(file_id);
        return -1;
    }
    
    // Set chunking with validated size
    hsize_t chunk_size = config->hdf5_chunk_size;
    if (chunk_size < 1024) chunk_size = 1024;  // Minimum chunk size
    if (chunk_size > 1048576) chunk_size = 1048576;  // Maximum chunk size (1MB)
    hsize_t chunk_dims[1] = {chunk_size};
    H5Pset_chunk(dcpl, 1, chunk_dims);
    
    if (config->hdf5_compression_level > 0) 
    {
        H5Pset_deflate(dcpl, config->hdf5_compression_level);
    }
    
    // Create initial dataspace
    hsize_t dims[1] = {0};
    hsize_t maxdims[1] = {H5S_UNLIMITED};
    hid_t space = H5Screate_simple(1, dims, maxdims);
    if (space < 0) 
    {
        //fprintf(stderr, "Failed to create dataspace\n");
        H5Pclose(dcpl);
        H5Tclose(compound_type);
        H5Fclose(file_id);
        return -1;
    }
    
    // Create dataset with validated properties
    hid_t dataset = H5Dcreate2(file_id, "methylation_data", compound_type, space,
                              H5P_DEFAULT, dcpl, H5P_DEFAULT);
    H5Pclose(dcpl);
    H5Sclose(space);
    
    if (dataset < 0) 
    {
        //fprintf(stderr, "Failed to create dataset\n");
        H5Tclose(compound_type);
        H5Fclose(file_id);
        return -1;
    }
    
    // Create metadata group and attributes
    hid_t meta_group = H5Gcreate2(file_id, "metadata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (meta_group >= 0) 
    {
        const char* version = "1.0";
        const char* format = "methylation_data";
        const char* date = __DATE__;
        H5LTset_attribute_string(file_id, "metadata", "version", version);
        H5LTset_attribute_string(file_id, "metadata", "format", format);
        H5LTset_attribute_string(file_id, "metadata", "date", date);
        H5LTset_attribute_string(file_id, "metadata", "chromosome", chrom);
        H5Gclose(meta_group);
    }
    
    // Allocate and initialize handles structure
    HDF5Handles* handles = calloc(1, sizeof(HDF5Handles));
    if (!handles) 
    {
        //fprintf(stderr, "Failed to allocate handles\n");
        H5Dclose(dataset);
        H5Tclose(compound_type);
        H5Fclose(file_id);
        return -1;
    }
    
    init_hdf5_handles(handles);  // Initialize with ref count 1
    handles->file_id = file_id;
    handles->dataset = dataset;
    handles->compound_type = compound_type;
    
    return (hid_t)handles;
}

static int write_to_hdf5(hid_t handles, MethylData* buffer, size_t count) 
{
    if (handles < 0 || !buffer || count == 0) return -1;
    
    HDF5Handles* h = (HDF5Handles*)handles;
    pthread_mutex_lock(&h->lock);  // Lock per-handle operations
    
    herr_t status = -1;
    
    if (h->dataset < 0 || h->compound_type < 0 || h->file_id < 0) 
    {
        pthread_mutex_unlock(&h->lock);
        return -1;
    }
    
    // Get current dimensions
    hid_t space = H5Dget_space(h->dataset);
    if (space < 0) 
    {
        pthread_mutex_unlock(&h->lock);
        return -1;
    }
    
    hsize_t dims[1], maxdims[1];
    if (H5Sget_simple_extent_dims(space, dims, maxdims) < 0) 
    {
        H5Sclose(space);
        pthread_mutex_unlock(&h->lock);
        return -1;
    }
    H5Sclose(space);
    
    // Extend the dataset
    hsize_t new_size = dims[0] + count;
    if (H5Dset_extent(h->dataset, &new_size) < 0) 
    {
        pthread_mutex_unlock(&h->lock);
        return -1;
    }
    
    // Get the space of extended dataset
    space = H5Dget_space(h->dataset);
    if (space < 0) 
    {
        pthread_mutex_unlock(&h->lock);
        return -1;
    }
    
    // Define hyperslab for writing
    hsize_t start[1] = {dims[0]};
    hsize_t block[1] = {count};
    if (H5Sselect_hyperslab(space, H5S_SELECT_SET, start, NULL, block, NULL) < 0) 
    {
        H5Sclose(space);
        pthread_mutex_unlock(&h->lock);
        return -1;
    }
    
    // Create memory space
    hid_t mem_space = H5Screate_simple(1, block, NULL);
    if (mem_space < 0) 
    {
        H5Sclose(space);
        pthread_mutex_unlock(&h->lock);
        return -1;
    }
    
    // Write the data with error checking
    status = H5Dwrite(h->dataset, h->compound_type, mem_space, space, H5P_DEFAULT, buffer);
    if (status < 0) 
    {
        //fprintf(stderr, "Failed to write data to HDF5 dataset (count: %zu)\n", count);
    }
    
    H5Sclose(mem_space);
    H5Sclose(space);
    
    // Ensure data is flushed to disk
    if (status >= 0) 
    {
        status = H5Fflush(h->file_id, H5F_SCOPE_LOCAL);
        if (status < 0) 
        {
            //fprintf(stderr, "Failed to flush HDF5 file to disk\n");
        }
    }
    
    pthread_mutex_unlock(&h->lock);
    return status;
}

static void close_hdf5_file(hid_t handles) 
{
    if (handles < 0) return;
    
    HDF5Handles* h = (HDF5Handles*)handles;
    pthread_mutex_lock(&handle_mutex);  // Global lock for handle management
    pthread_mutex_lock(&h->lock);       // Per-handle lock
    
    // Only flush if handles are still valid
    if (h->file_id >= 0) 
    {
        herr_t status = H5Fflush(h->file_id, H5F_SCOPE_GLOBAL);
        if (status < 0) 
        {
            //fprintf(stderr, "Warning: Final H5Fflush failed\n");
        }
    }
    
    pthread_mutex_unlock(&h->lock);
    pthread_mutex_unlock(&handle_mutex);
    
    unref_hdf5_handles(h);  // This will handle the actual closing
}

static int ensure_directory_exists(const char* path) 
{
    if (mkdir(path, 0755) != 0 && errno != EEXIST) 
    {
        //fprintf(stderr, "Failed to create directory %s: %s\n", path, strerror(errno));
        return -1;
    }
    return 0;
}

// Modify the check_chromosome_in_fai function to handle chromosome name mismatches
static int check_chromosome_in_fai(faidx_t *fai, const char *chrom) 
{
    if (fai == NULL) 
    {
        //fprintf(stderr, "Error: FASTA index (fai) is NULL in check_chromosome_in_fai\n");
        return 0;
    }
    if (chrom == NULL) 
    {
        //fprintf(stderr, "Error: Chromosome name is NULL in check_chromosome_in_fai\n");
        return 0;
    }
    int len = faidx_seq_len(fai, chrom);
    if (len >= 0) 
    {
        return 1;
    }
    // If not found, try converting the name (e.g., add or remove 'chr' prefix)
    char alt_chrom[256];
    if (strncmp(chrom, "chr", 3) == 0) 
    {
        // Remove 'chr' prefix
        strncpy(alt_chrom, chrom + 3, sizeof(alt_chrom) - 1);
        alt_chrom[sizeof(alt_chrom) - 1] = '\0';
        len = faidx_seq_len(fai, alt_chrom);
        if (len >= 0) 
        {
            //fprintf(stderr, "Note: Chromosome '%s' not found, but '%s' was found in FASTA file\n", chrom, alt_chrom);
            return 1;
        }
    } 
    else 
    {
        // Add 'chr' prefix
        snprintf(alt_chrom, sizeof(alt_chrom), "chr%s", chrom);
        len = faidx_seq_len(fai, alt_chrom);
        if (len >= 0) 
        {
            //fprintf(stderr, "Note: Chromosome '%s' not found, but '%s' was found in FASTA file\n", chrom, alt_chrom);
            return 1;
        }
    }
    //fprintf(stderr, "Error: Chromosome '%s' not found in FASTA file (tried '%s' as well)\n", chrom, alt_chrom);
    return 0;
}

// Add a structure for work units
typedef struct 
{
    char chrom[32];          // Chromosome name
    hts_pos_t start_pos;     // Start position of the chunk
    hts_pos_t end_pos;       // End position of the chunk
    int tid;                 // Target ID from BAM header
} WorkUnit;

// Add a thread-safe queue for work units
typedef struct 
{
    WorkUnit* units;         // Array of work units
    size_t capacity;         // Total capacity of the queue
    size_t head;             // Head index for dequeuing
    size_t tail;             // Tail index for enqueuing
    size_t current_count;    // Current number of items in queue
    char current_chrom[32];  // Current chromosome being processed
    int is_current_chrom_done; // Flag to indicate if current chromosome is done
    pthread_mutex_t mutex;   // Mutex for thread safety
    pthread_cond_t not_empty; // Condition variable for not empty
    pthread_cond_t not_full;  // Condition variable for not full
    pthread_cond_t chrom_done;// Condition variable for chromosome completion
} WorkQueue;

static WorkQueue work_queue;

// Dequeue a work unit, ensuring it's for the current chromosome
static int dequeue_work(WorkUnit* unit) 
{
    pthread_mutex_lock(&work_queue.mutex);
    while (work_queue.current_count == 0 && !work_queue.is_current_chrom_done) 
    {
        pthread_cond_wait(&work_queue.not_empty, &work_queue.mutex);
    }
    if (work_queue.current_count == 0) 
    {
        pthread_mutex_unlock(&work_queue.mutex);
        //fprintf(stderr, "Done, no more chromosomes to process\n");
        return -1; // Queue is empty and done
    }
    // Check if the next unit is for the current chromosome
    if (work_queue.current_chrom[0] == '\0' || strcmp(work_queue.current_chrom, work_queue.units[work_queue.head].chrom) == 0) 
    {
        *unit = work_queue.units[work_queue.head];
        if (work_queue.current_chrom[0] == '\0') 
        {
            strncpy(work_queue.current_chrom, unit->chrom, sizeof(work_queue.current_chrom) - 1);
            work_queue.current_chrom[sizeof(work_queue.current_chrom) - 1] = '\0';
            //fprintf(stderr, "Starting processing of chromosome %s\n", work_queue.current_chrom);
        }
        work_queue.head = (work_queue.head + 1) % work_queue.capacity;
        work_queue.current_count--;
        pthread_cond_signal(&work_queue.not_full);
        pthread_mutex_unlock(&work_queue.mutex);
        // Validate chromosome name
        if (unit->chrom[0] == '\0') 
        {
            //fprintf(stderr, "Error: Empty chromosome name in work unit\n");
            return -1;
        }
        return 0;
    } 
    else 
    {
        // Wait for current chromosome to be fully processed
        while (!work_queue.is_current_chrom_done) 
        {
            pthread_cond_wait(&work_queue.chrom_done, &work_queue.mutex);
        }
        // Check if there are more work units after current chromosome is done
        if (work_queue.current_count == 0) 
        {
            pthread_mutex_unlock(&work_queue.mutex);
            //fprintf(stderr, "Done, no more chromosomes to process\n");
            return -1; // No more work units, exit
        }
        // Transition to the next chromosome
        strncpy(work_queue.current_chrom, work_queue.units[work_queue.head].chrom, sizeof(work_queue.current_chrom) - 1);
        work_queue.current_chrom[sizeof(work_queue.current_chrom) - 1] = '\0';
        work_queue.is_current_chrom_done = 0;
        //fprintf(stderr, "Transitioning to chromosome %s\n", work_queue.current_chrom);
        *unit = work_queue.units[work_queue.head];
        work_queue.head = (work_queue.head + 1) % work_queue.capacity;
        work_queue.current_count--;
        pthread_cond_signal(&work_queue.not_full);
        pthread_mutex_unlock(&work_queue.mutex);
        // Validate chromosome name
        if (unit->chrom[0] == '\0') 
        {
            //fprintf(stderr, "Error: Empty chromosome name in work unit\n");
            return -1;
        }
        return 0;
    }
}

// Signal that a chromosome is done
static void signal_chrom_done(const char* chrom) 
{
    pthread_mutex_lock(&work_queue.mutex);
    if (strcmp(work_queue.current_chrom, chrom) == 0) 
    {
        // Check if all work units for the current chromosome are dequeued
        int remaining_units_for_chrom = 0;
        size_t index = work_queue.head;
        for (size_t i = 0; i < work_queue.current_count; i++) 
        {
            if (strcmp(work_queue.units[index].chrom, chrom) == 0) 
            {
                remaining_units_for_chrom++;
            }
            index = (index + 1) % work_queue.capacity;
        }
        if (remaining_units_for_chrom == 0) 
        {
            work_queue.is_current_chrom_done = 1;
            //fprintf(stderr, "All work units for chromosome %s are processed, signaling completion\n", chrom);
            pthread_cond_broadcast(&work_queue.chrom_done);
            // If no more work units in queue, signal all threads to exit
            if (work_queue.current_count == 0) 
            {
                ///fprintf(stderr, "Work queue is empty, no more chromosomes to process, signaling all threads to exit\n");
                pthread_cond_broadcast(&work_queue.not_empty);
            }
        } 
        else 
        {
            //fprintf(stderr, "Chromosome %s processing not complete, %d work units remaining\n", chrom, remaining_units_for_chrom);
        }
    }
    pthread_mutex_unlock(&work_queue.mutex);
}

// Add a structure for workers (threads) with memory buffers
typedef struct 
{
    int thread_id;           // Unique identifier for the thread
    char current_chrom[32];  // Current chromosome being processed by this thread
    hts_pos_t start_pos;     // Start position of the current chunk
    hts_pos_t end_pos;       // End position of the current chunk
    hid_t hdf5_handles;      // HDF5 file handles for the current chromosome
    MethylData* buffer;      // Buffer to store methylation data
    size_t buffer_pos;       // Current position in the buffer
    int had_error;           // Flag to indicate if an error occurred
} Worker;

// Modify ChromBuffer structure to include position tracking for ordered writing (around line 700)
typedef struct 
{
    MethylData* data;        // Dynamic array to store methylation data
    size_t capacity;         // Total capacity of the buffer
    size_t size;             // Current number of records in buffer
    pthread_mutex_t mutex;   // Mutex for thread-safe access
    int is_processing_done;  // Flag to indicate if processing for this chromosome is complete
    int active_threads;      // Number of threads currently processing this chromosome
    hts_pos_t* positions;    // Array to store positions for sorting if needed
    size_t pos_capacity;     // Capacity of positions array
    size_t pos_size;         // Current size of positions array
} ChromBuffer;

// Add a global map of chromosome buffers (around line 210)
#define MAX_CHROMOSOMES 32
// Replace array with single global buffer
static ChromBuffer global_chrom_buffer;

// Add global variables for chromosome sequences
// Define ChromSequence struct before global variables
typedef struct 
{
    char* chrom_name;
    char* sequence;
    int seq_len;
} ChromSequence;

static ChromSequence* chrom_sequences = NULL;
static int num_chromosomes = 0;

// Initialize chromosome buffers for all possible chromosomes with large capacity
static void initialize_global_chrom_buffer(size_t initial_capacity) 
{
    ChromBuffer* buffer = &global_chrom_buffer;
    if (buffer->capacity == 0) 
    {
        buffer->data = malloc(initial_capacity * sizeof(MethylData));
        if (!buffer->data) 
        {
            //fprintf(stderr, "Failed to allocate memory for global chromosome buffer\n");
            return;
        }
        buffer->positions = malloc(initial_capacity * sizeof(hts_pos_t));
        if (!buffer->positions) 
        {
            //fprintf(stderr, "Failed to allocate memory for positions array in global chromosome buffer\n");
            free(buffer->data);
            buffer->data = NULL;
            buffer->capacity = 0;
            return;
        }
        buffer->capacity = initial_capacity;
        buffer->pos_capacity = initial_capacity;
        buffer->size = 0;
        buffer->pos_size = 0;
        buffer->is_processing_done = 0;
        buffer->active_threads = 0;
        pthread_mutex_init(&buffer->mutex, NULL);
        //fprintf(stderr, "Initialized global chromosome buffer with capacity %zu\n", initial_capacity);
    }
}

// Free global chromosome buffer at program end
static void free_global_chrom_buffer() 
{
    ChromBuffer* buffer = &global_chrom_buffer;
    if (buffer->data) 
    {
        free(buffer->data);
        buffer->data = NULL;
    }
    if (buffer->positions) 
    {
        free(buffer->positions);
        buffer->positions = NULL;
    }
    buffer->capacity = 0;
    buffer->pos_capacity = 0;
    buffer->size = 0;
    buffer->pos_size = 0;
    pthread_mutex_destroy(&buffer->mutex);
    //fprintf(stderr, "Freed global chromosome buffer\n");
}

// Update add_to_chrom_buffer to store position for potential sorting (around line 760)
static int add_to_chrom_buffer(ChromBuffer* buffer, MethylData* mdata) 
{
    pthread_mutex_lock(&buffer->mutex);
    if (buffer->size >= buffer->capacity) 
    {
        size_t new_capacity = buffer->capacity * 2;
        MethylData* new_data = realloc(buffer->data, new_capacity * sizeof(MethylData));
        if (!new_data) 
        {
            fprintf(stderr, "Failed to resize global chromosome buffer from %zu to %zu\n", buffer->capacity, new_capacity);
            pthread_mutex_unlock(&buffer->mutex);
            return -1;
        }
        hts_pos_t* new_positions = realloc(buffer->positions, new_capacity * sizeof(hts_pos_t));
        if (!new_positions) 
        {
            free(new_data);
            pthread_mutex_unlock(&buffer->mutex);
            return -1;
        }
        buffer->data = new_data;
        buffer->positions = new_positions;
        buffer->capacity = new_capacity;
        buffer->pos_capacity = new_capacity;
    }
    buffer->data[buffer->size] = *mdata;
    buffer->positions[buffer->pos_size] = mdata->pos;
    buffer->size++;
    buffer->pos_size++;
    pthread_mutex_unlock(&buffer->mutex);
    return 0;
}

static ChromBuffer* current_buffer_for_sorting = NULL;
// Add global variables for tracking processed chromosomes for progress reporting
static int processed_chromosomes = 0;
static int32_t total_chromosomes_global;
static pthread_mutex_t progress_mutex = PTHREAD_MUTEX_INITIALIZER;

static int compare_positions(const void *a, const void *b) 
{
    size_t idx_a = *(size_t *)a;
    size_t idx_b = *(size_t *)b;
    hts_pos_t pos_a = current_buffer_for_sorting->positions[idx_a];
    hts_pos_t pos_b = current_buffer_for_sorting->positions[idx_b];
    return (pos_a > pos_b) - (pos_a < pos_b);
}

// Modify signal_thread_completion to update progress bar when HDF5 file is written
static int signal_thread_completion(ChromBuffer* buffer, const char* chrom, Config* config) 
{
    pthread_mutex_lock(&buffer->mutex);
    buffer->active_threads--;
    if (buffer->active_threads == 0) 
    {
        buffer->is_processing_done = 1;
        if (buffer->size > 0) 
        {
            // Create an array of indices for sorting
            size_t* indices = malloc(buffer->size * sizeof(size_t));
            if (!indices) 
            {
                pthread_mutex_unlock(&buffer->mutex);
                return -1;
            }
            for (size_t i = 0; i < buffer->size; i++) 
                indices[i] = i;

            // Use qsort for efficient sorting
            current_buffer_for_sorting = buffer;
            qsort(indices, buffer->size, sizeof(size_t), compare_positions);
            current_buffer_for_sorting = NULL;
            // Create a temporary sorted array
            MethylData* sorted_data = malloc(buffer->size * sizeof(MethylData));
            if (!sorted_data) 
            {
                free(indices);
                pthread_mutex_unlock(&buffer->mutex);
                return -1;
            }
            for (size_t i = 0; i < buffer->size; i++) 
                sorted_data[i] = buffer->data[indices[i]];
            free(indices);
            // Write sorted data to HDF5
            hid_t hdf5_handles = create_hdf5_file(chrom, config->hdf5_output_dir, config);
            if (hdf5_handles < 0) 
            {
                free(sorted_data);
                pthread_mutex_unlock(&buffer->mutex);
                return -1;
            }
            if (write_to_hdf5(hdf5_handles, sorted_data, buffer->size) < 0) 
            {
                close_hdf5_file(hdf5_handles);
                free(sorted_data);
                pthread_mutex_unlock(&buffer->mutex);
                return -1;
            }
            close_hdf5_file(hdf5_handles);
            free(sorted_data);
        } 
        else 
        {
            // If no data, still create an empty HDF5 file for consistency
            hid_t hdf5_handles = create_hdf5_file(chrom, config->hdf5_output_dir, config);
            if (hdf5_handles < 0) 
            {
                pthread_mutex_unlock(&buffer->mutex);
                return -1;
            }
            close_hdf5_file(hdf5_handles);
        }
        // Reset buffer counters for the next chromosome
        buffer->size = 0;
        buffer->pos_size = 0;
        pthread_mutex_lock(&progress_mutex);
        processed_chromosomes++;
        fprintf(stderr, "=");
        // Check if all chromosomes are processed to close the progress bar
        if (processed_chromosomes >= total_chromosomes_global) 
        {
            fprintf(stderr, "] Done\n");
        }
        pthread_mutex_unlock(&progress_mutex);
    }
    pthread_mutex_unlock(&buffer->mutex);
    return 0;
}

// Add missing plp_get_read function (around line 600, after aux_t structure definition)
static int plp_get_read(void *data, bam1_t *b) 
{
    aux_t *aux = (aux_t*)data;
    int rv = 0;
    
    while (rv >= 0) 
    {
        if ((rv = aux->iter ? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, NULL, b)) >= 0) 
        {
            if (b->core.tid < 0) return -1;  // unmapped
            if (aux->tid >= 0 && (aux->tid != b->core.tid || b->core.pos >= aux->end)) return -1;
        }
    }
    return rv;
}

// Add missing plp_set_region function (around line 605, after plp_get_read function)
static void plp_set_region(void *data, int tid, hts_pos_t start, hts_pos_t end) 
{
    aux_t *aux = (aux_t*)data;
    aux->tid = tid;
    aux->start = start;
    aux->end = end;
}

void *extractCalls(void *foo) 
{
    Config *config = (Config*)foo;
    bam_hdr_t *hdr = config->hdr;
    bam_mplp_t iter = NULL;
    int ret, tid, i, seqlen, type, rv;
    hts_pos_t pos;
    int n_plp;
    int strand, direction;
    uint32_t nmethyl = 0, nunmethyl = 0;
    const bam_pileup1_t **plp = NULL;
    char *seq = NULL, base = 'A';
    htsFile *fp = NULL;
    char* current_raw_chrom = NULL;
    char* current_std_chrom = NULL;
    char* current_fasta_chrom = NULL;
    aux_t **data = NULL;
    
    // Initialize Worker structure for this thread
    Worker worker;
    worker.thread_id = (int)(unsigned long)pthread_self();
    worker.current_chrom[0] = '\0';
    worker.start_pos = -1;
    worker.end_pos = -1;
    worker.hdf5_handles = -1;
    worker.buffer_pos = 0;
    worker.had_error = 0;
    
    // Register thread
    pthread_mutex_lock(&chrom_transition.mutex);
    chrom_transition.active_threads++;
    pthread_mutex_unlock(&chrom_transition.mutex);
    
    // Open a thread-local BAM file handle only once per thread
    fp = hts_open(config->BAMName, "rb");
    if (!fp) 
    {
        worker.had_error = 1;
        goto cleanup;
    }

    // Set thread-local BAM handle options
    hts_set_opt(fp, HTS_OPT_CACHE_SIZE, 50000000);
    hts_set_opt(fp, HTS_OPT_BLOCK_SIZE, 1 << 20);

    // Allocate data for the pileup engine only once
    data = calloc(1, sizeof(aux_t*));
    if (!data) 
    {
        worker.had_error = 1;
        goto cleanup;
    }
    data[0] = calloc(1, sizeof(aux_t));
    if (!data[0]) 
    {
        worker.had_error = 1;
        goto cleanup;
    }
    data[0]->fp = fp;
    
    plp = calloc(1, sizeof(bam_pileup1_t *));
    if (plp == NULL) 
    {
        worker.had_error = 1;
        goto cleanup;
    }
    
    while (1) 
    {
        WorkUnit unit;
        if (dequeue_work(&unit) != 0) 
        {
            if (current_std_chrom != NULL) 
                signal_chrom_done(current_std_chrom);
            break;
        }
        
        // Extract details from the work unit
        tid = unit.tid;
        hts_pos_t localPos = unit.start_pos;
        hts_pos_t localEnd = unit.end_pos;
        const char* std_chrom = unit.chrom;
        const char* raw_chrom = hdr->target_name[tid];
        
        // Update worker's current positions
        worker.start_pos = localPos;
        worker.end_pos = localEnd;
        
        // If we transitioned chromosomes OR if current chromosome names are not set
        if (worker.current_chrom[0] == '\0' || strcmp(worker.current_chrom, std_chrom) != 0) 
        {
            // No need to flush local buffer since we use shared buffer
            worker.buffer_pos = 0; // Reset position, though not used
            
            // Update current chromosome
            strncpy(worker.current_chrom, std_chrom, sizeof(worker.current_chrom) - 1);
            worker.current_chrom[sizeof(worker.current_chrom) - 1] = '\0';
            
            pthread_mutex_lock(&global_chrom_buffer.mutex);
            global_chrom_buffer.active_threads++;
            pthread_mutex_unlock(&global_chrom_buffer.mutex);
            
            // Free and update chromosome names
            if (current_raw_chrom) free(current_raw_chrom);
            if (current_std_chrom) free(current_std_chrom);
            if (current_fasta_chrom) free(current_fasta_chrom);
            current_raw_chrom = strdup(raw_chrom);
            current_std_chrom = strdup(std_chrom);
            current_fasta_chrom = NULL; // Will be set after checking FASTA file
            if (!current_raw_chrom || !current_std_chrom) 
            {
                worker.had_error = 1;
                continue;
            }
        }
        
        // Process the chunk
        hts_pos_t localPos2 = localPos > 1 ? localPos - 2 : 0;
        
        // Verify chromosome exists in FASTA index using standardized name
        if (!check_chromosome_in_fai(config->fai, current_std_chrom)) 
        {
            worker.had_error = 1;
            continue;
        }

        // Determine the correct chromosome name to use with FASTA file
        if (!current_fasta_chrom) 
        {
            int len = faidx_seq_len(config->fai, current_std_chrom);
            if (len >= 0) 
            {
                current_fasta_chrom = strdup(current_std_chrom);
            } 
            else 
            {
                char alt_chrom[256];
                if (strncmp(current_std_chrom, "chr", 3) == 0) 
                {
                    strncpy(alt_chrom, current_std_chrom + 3, sizeof(alt_chrom) - 1);
                    alt_chrom[sizeof(alt_chrom) - 1] = '\0';
                    len = faidx_seq_len(config->fai, alt_chrom);
                    if (len >= 0) 
                    {
                        current_fasta_chrom = strdup(alt_chrom);
                    }
                } 
                else 
                {
                    snprintf(alt_chrom, sizeof(alt_chrom), "chr%s", current_std_chrom);
                    len = faidx_seq_len(config->fai, alt_chrom);
                    if (len >= 0) 
                    {
                        current_fasta_chrom = strdup(alt_chrom);
                    }
                }
            }
            if (!current_fasta_chrom) 
            {
                worker.had_error = 1;
                continue;
            }
        }

        // Get the chromosome length from FASTA index to check bounds
        int chrom_len = faidx_seq_len(config->fai, current_fasta_chrom);
        if (chrom_len < 0) 
        {
            worker.had_error = 1;
            continue;
        }

        // Check if the requested coordinates are within bounds
        if (localEnd > chrom_len) 
            localEnd = chrom_len;

        // Use preloaded sequence if available, otherwise fetch it
        seq = NULL;
        seqlen = 0;
        for (int i = 0; i < num_chromosomes; i++) 
        {
            if (chrom_sequences && strcmp(current_std_chrom, chrom_sequences[i].chrom_name) == 0 && chrom_sequences[i].sequence) 
            {
                seq = chrom_sequences[i].sequence + localPos2;
                seqlen = chrom_sequences[i].seq_len - localPos2;
                if (seqlen > localEnd - localPos2) 
                    seqlen = localEnd - localPos2;
                break;
            }
        }
        if (!seq) 
        {
            // Fetch sequence only for the required range
            hts_pos_t fetch_end = localEnd;
            if (fetch_end > chrom_len) fetch_end = chrom_len;
            seq = faidx_fetch_seq(config->fai, current_fasta_chrom, localPos2, fetch_end - 1, &seqlen);
            if (seqlen < 0) 
            {
                worker.had_error = 1;
                pthread_mutex_lock(&chrom_transition.mutex);
                chrom_transition.threads_in_transition++;
                pthread_cond_broadcast(&chrom_transition.cond);
                pthread_mutex_unlock(&chrom_transition.mutex);
                continue;
            }
        }
        
        // Set up region for pileup using predefined chunks
        data[0]->tid = tid;
        data[0]->start = localPos;
        data[0]->end = localEnd;
        if (data[0]->iter) 
        {
            hts_itr_destroy(data[0]->iter);
            data[0]->iter = NULL;
        }
        data[0]->iter = sam_itr_queryi(config->bai, tid, localPos, localEnd);
        if (!data[0]->iter) 
        {
            worker.had_error = 1;
            pthread_mutex_lock(&chrom_transition.mutex);
            chrom_transition.threads_in_transition++;
            pthread_cond_broadcast(&chrom_transition.cond);
            pthread_mutex_unlock(&chrom_transition.mutex);
            continue;
        }

        // Initialize pileup engine only if not already initialized
        if (iter) 
        {
            bam_mplp_destroy(iter);
            iter = NULL;
        }
        iter = bam_mplp_init(1, plp_get_read, (void**)data);
        if (!iter) 
        {
            worker.had_error = 1;
            pthread_mutex_lock(&chrom_transition.mutex);
            chrom_transition.threads_in_transition++;
            pthread_cond_broadcast(&chrom_transition.cond);
            pthread_mutex_unlock(&chrom_transition.mutex);
            continue;
        }
        bam_mplp_set_maxcnt(iter, INT_MAX);
        plp_set_region(data[0], tid, localPos, localEnd);

        // Process pileup with a timeout to prevent stalls, but no iteration limit
        struct timespec pileup_start_time, current_time;
        clock_gettime(CLOCK_MONOTONIC, &pileup_start_time);
        long pileup_timeout_seconds = 300;
        while ((ret = bam_mplp64_auto(iter, &tid, &pos, &n_plp, plp)) > 0) 
        {          
            // Check for timeout
            clock_gettime(CLOCK_MONOTONIC, &current_time);
            long elapsed_seconds = current_time.tv_sec - pileup_start_time.tv_sec;
            if (elapsed_seconds > pileup_timeout_seconds) 
            {
                worker.had_error = 1;
                pthread_mutex_lock(&chrom_transition.mutex);
                chrom_transition.threads_in_transition++;
                pthread_cond_broadcast(&chrom_transition.cond);
                pthread_mutex_unlock(&chrom_transition.mutex);
                break;
            }
            
            if (pos < localPos || pos >= localEnd) continue;
            
            // Check methylation context
            if ((direction = isCpG(seq, pos-localPos2, seqlen))) 
            {
                if (!config->keepCpG) continue;
                type = 0;
            } 
            else if ((direction = isCHG(seq, pos-localPos2, seqlen))) 
            {
                if (!config->keepCHG) continue;
                type = 1;
            } 
            else if ((direction = isCHH(seq, pos-localPos2, seqlen))) 
            {
                if (!config->keepCHH) continue;
                type = 2;
            } 
            else 
                continue;
            
            // Count methylation
            nmethyl = nunmethyl = 0;
            base = *(seq + pos - localPos2);
            for (i = 0; i < n_plp; i++) 
            {
                if (plp[0][i].is_del || plp[0][i].is_refskip) continue;
                strand = getStrand(plp[0][i].b);
                if (strand & 1) 
                {
                    if (base != 'C' && base != 'c') continue;
                } 
                else 
                {
                    if (base != 'G' && base != 'g') continue;
                }
                rv = updateMetrics(config, &plp[0][i]);
                if (rv > 0) 
                    nmethyl++;
                else if (rv < 0) 
                    nunmethyl++;
            }

            // Check depth threshold
            uint32_t depth = nmethyl + nunmethyl;
            if (depth < config->minDepth || depth > config->maxDepth) continue;
            
            // Create methylation data record
            MethylData mdata;
            setMethylFlags(
                &mdata, 
                chromToNumber(hdr->target_name[tid]),
                encodeStrand((base == 'C' || base == 'c') ? '+' : '-'),
                type
            );
            mdata.pos = pos;
            mdata.nmethyl = nmethyl;
            mdata.nunmethyl = nunmethyl;
            mdata.tnc = encodeTrinucleotide(seq, pos - localPos2, seqlen, direction);
            
            // Add to shared chromosome buffer instead of local buffer or direct HDF5 write
            if (add_to_chrom_buffer(&global_chrom_buffer, &mdata) < 0) 
            {
                worker.had_error = 1;
                break;
            }
        }
        
        // Signal that this thread has finished processing a chunk for the current chromosome
        signal_chrom_done(current_std_chrom);

        // Signal completion for this thread on the current chromosome
        if (signal_thread_completion(&global_chrom_buffer, current_std_chrom, config) < 0) 
        {
            worker.had_error = 1;
        }

        // Check for chunk processing timeout before looping back
        struct timespec chunk_start_time, chunk_current_time;
        clock_gettime(CLOCK_MONOTONIC, &chunk_start_time);
        long chunk_timeout_seconds = 600;

        // Check for chunk processing timeout before looping back
        clock_gettime(CLOCK_MONOTONIC, &chunk_current_time);
        long chunk_elapsed_seconds = chunk_current_time.tv_sec - chunk_start_time.tv_sec;
        if (chunk_elapsed_seconds > chunk_timeout_seconds) 
        {
            //fprintf(stderr, "Thread %d: Exceeded timeout of %ld seconds for processing chunk of chromosome %s (range %"PRIi64"-%"PRIi64"), last processed pos: %"PRIi64", aborting to prevent stall\n",
            //        worker.thread_id, chunk_timeout_seconds, current_std_chrom, localPos, localEnd, last_processed_pos);
            worker.had_error = 1;
            pthread_mutex_lock(&chrom_transition.mutex);
            chrom_transition.threads_in_transition++;
            //fprintf(stderr, "Thread %d: Aborted processing chunk for chromosome %s due to timeout, marking as done (threads in transition: %d, working threads: %d)\n",
            //        worker.thread_id, current_std_chrom, chrom_transition.threads_in_transition, chrom_transition.working_threads);
            pthread_cond_broadcast(&chrom_transition.cond);
            pthread_mutex_unlock(&chrom_transition.mutex);
            break;
        }
    }

cleanup:
    // Ensure final buffer is flushed only if the file handle is still valid - not needed with shared buffer
    worker.buffer_pos = 0;
    
    // Clean up resources
    if (current_std_chrom != NULL) 
    {
        free(current_raw_chrom);
        free(current_std_chrom);
        free(current_fasta_chrom);
    }
    
    // No local buffer to free
    // if (worker.buffer) free(worker.buffer); // Removed as no local buffer is allocated
    
    if (iter) bam_mplp_destroy(iter);
    if (data && data[0]) 
    {
        if (data[0]->iter) hts_itr_destroy(data[0]->iter);
        free(data[0]);
    }
    free(data);
    free(plp);
    if (fp) hts_close(fp);
    
    // Deregister thread and update error state
    pthread_mutex_lock(&chrom_transition.mutex);
    chrom_transition.active_threads--;
    if (worker.had_error) 
    {
        chrom_transition.error_count++;
    }
    if (chrom_transition.threads_in_transition > 0) 
    {
        pthread_cond_broadcast(&chrom_transition.cond);
    }
    pthread_mutex_unlock(&chrom_transition.mutex);
    
    return (void*)(intptr_t)(worker.had_error ? 1 : 0);
}

int main(int argc, char *argv[]) 
{
    // Initialize HDF5 with thread safety
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    hbool_t is_threadsafe;
    herr_t thread_safe = H5is_library_threadsafe(&is_threadsafe);
    if (thread_safe < 0 || !is_threadsafe) 
    {
        //fprintf(stderr, "HDF5 library was not built with thread safety\n");
        return 1;
    }

    Config config;
    int c;
    
    config.keepCpG = 1;
    config.keepCHG = 0;  // Default to 0 (false)
    config.keepCHH = 0;  // Default to 0 (false)
    config.minMapq = 30; // Default minimum mapping quality
    config.minPhred = 20;// Default minimum Phred score
    config.minDepth = 4; // Default minimum depth
    config.maxDepth = 100; // Default maximum depth 
    config.nThreads = 8;
    config.chunkSize = 1000000;
    config.bufferSizePerThread = 100000;
    config.hdf5_output_dir = NULL;
    config.hdf5_compression_level = 9;
    config.hdf5_chunk_size = 1000;
    
    static struct option lopts[] = 
    {
        {"o", 1, NULL, 1},
        {"hdf5-compression", 1, NULL, 2},
        {"hdf5-chunk-size", 1, NULL, 3},
        {"@", 1, NULL, '@'},
        {"CHG", 0, NULL, 'g'},
        {"CHH", 0, NULL, 'h'},
        {"q", 1, NULL, 'q'},
        {"p", 1, NULL, 'p'},
        {"c", 1, NULL, 'c'},
        {"C", 1, NULL, 'C'},
        {"help", 0, NULL, 'H'},
        {0, 0, NULL, 0}
    };
    
    while ((c = getopt_long(argc, argv, "@ghqpcCH", lopts, NULL)) >= 0) 
    {
        switch (c) 
        {
            case 1: config.hdf5_output_dir = strdup(optarg); break;
            case 2:
                config.hdf5_compression_level = atoi(optarg);
                if (config.hdf5_compression_level < 0 || config.hdf5_compression_level > 9) 
                {
                    fprintf(stderr, "Compression level must be between 0 and 9\n");
                    return 1;
                }
                break;
            case 3:
                config.hdf5_chunk_size = atoi(optarg);
                if (config.hdf5_chunk_size < 1) 
                {
                    fprintf(stderr, "Chunk size must be at least 1\n");
                    return 1;
                }
                break;
            case '@': 
                config.nThreads = atoi(optarg); 
                break;
            case 'g':
                config.keepCHG = 1;
                break;
            case 'h':
                config.keepCHH = 1;
                break;
            case 'q':
                config.minMapq = atoi(optarg);
                break;
            case 'p':
                config.minPhred = atoi(optarg);
                break;
            case 'c':
                config.minDepth = atoi(optarg);
                break;
            case 'C':
                config.maxDepth = atoi(optarg);
                break;
            case 'H':
                fprintf(stderr, "Usage: %s [options] <ref.fa> <sorted_alignments.bam>\n", argv[0]);
                fprintf(stderr, "Options:\n");
                fprintf(stderr, "  --o DIR                  Output directory for HDF5 files\n");
                fprintf(stderr, "  --hdf5-compression INT   Compression level (0-9, default: 9)\n");
                fprintf(stderr, "  --hdf5-chunk-size INT    Chunk size for HDF5 datasets (default: 1000)\n");
                fprintf(stderr, "  --@ INT                  Number of threads to use (overrides auto-detection)\n");
                fprintf(stderr, "  --CHG                    Keep CHG context methylation data (default: false)\n");
                fprintf(stderr, "  --CHH                    Keep CHH context methylation data (default: false)\n");
                fprintf(stderr, "  --q INT                  Minimum mapping quality (default: 30)\n");
                fprintf(stderr, "  --p INT                  Minimum Phred score (default: 20)\n");
                fprintf(stderr, "  --c INT                  Minimum coverage (default: 4)\n");
                fprintf(stderr, "  --C INT                  Maximum coverage (default: 100)\n");
                return 0;
        }
    }
    
    if (argc - optind < 2) 
    {
        fprintf(stderr, "Usage: %s [options] <ref.fa> <sorted_alignments.bam>\n", argv[0]);
        return 1;
    }
    
    if (!config.hdf5_output_dir) 
    {
        fprintf(stderr, "Error: HDF5 output directory must be specified with --hdf5-output-dir\n");
        return 1;
    }
    
    if (access(argv[optind], R_OK) != 0) 
    { 
        fprintf(stderr, "Cannot read FASTA file %s: %s\n", argv[optind], strerror(errno));
        return 1;
    }
    if (access(argv[optind+1], R_OK) != 0) 
    {
        fprintf(stderr, "Cannot read BAM file %s: %s\n", argv[optind+1], strerror(errno));
        return 1;
    }
    
    if (ensure_directory_exists(config.hdf5_output_dir) != 0) return 1;
    
    config.FastaName = argv[optind];
    config.BAMName = argv[optind+1];
    
    // Print simple header with extraction details
    fprintf(stderr, "Methylation Extraction Pipeline\n");
    fprintf(stderr, "--------------------------------\n");
    fprintf(stderr, "Input BAM: %s\n", config.BAMName);
    fprintf(stderr, "Input FASTA: %s\n", config.FastaName);
    fprintf(stderr, "Output Directory: %s\n", config.hdf5_output_dir);
    fprintf(stderr, "Keep CpG: %d\n", config.keepCpG);
    fprintf(stderr, "Keep CHG: %d\n", config.keepCHG);
    fprintf(stderr, "Keep CHH: %d\n", config.keepCHH);
    fprintf(stderr, "Minimum mapping quality: %d\n", config.minMapq);
    fprintf(stderr, "Minimum Phred score: %d\n", config.minPhred);
    fprintf(stderr, "Minimum depth: %d\n", config.minDepth);
    fprintf(stderr, "Maximum depth: %d\n", config.maxDepth);
    fprintf(stderr, "Number of threads: %d\n", config.nThreads);
    fprintf(stderr, "Chunk size: %zu\n", config.chunkSize);
    fprintf(stderr, "Buffer size per thread: %zu\n", config.bufferSizePerThread);
    fprintf(stderr, "HDF5 compression level: %d\n", config.hdf5_compression_level);
    fprintf(stderr, "HDF5 chunk size: %d\n", config.hdf5_chunk_size);
 
    fprintf(stderr, "Progress: [");
    
    // Load indexes and header before creating threads
    htsFile *fp = hts_open(config.BAMName, "rb");
    if (!fp) 
    {
        fprintf(stderr, "Could not open BAM file %s\n", config.BAMName);
        return 1;
    }
    
    config.bai = sam_index_load(fp, config.BAMName);
    if (config.bai == NULL) 
    {
        fprintf(stderr, "Index not found for %s, attempting to create it...\n", config.BAMName);
        hts_close(fp);  // Close file before index creation
        
        if (sam_index_build(config.BAMName, 0) < 0) 
        {
            fprintf(stderr, "Failed to create index for %s\n", config.BAMName);
            return 1;
        }
        
        // Reopen file and load the newly created index
        fp = hts_open(config.BAMName, "rb");
        if (!fp) 
        {
            fprintf(stderr, "Could not reopen BAM file %s\n", config.BAMName);
            return 1;
        }
        
        config.bai = sam_index_load(fp, config.BAMName);
        if (config.bai == NULL) 
        {
            fprintf(stderr, "Failed to load newly created index for %s\n", config.BAMName);
            hts_close(fp);
            return 1;
        }
    }
    
    config.hdr = sam_hdr_read(fp);
    if (config.hdr == NULL) 
    {
        fprintf(stderr, "Could not read header from %s\n", config.BAMName);
        hts_idx_destroy(config.bai);
        hts_close(fp);
        return 1;
    }

    // After loading the BAM header
    bam_hdr_t *hdr = config.hdr; 
    hts_pos_t max_chrom_len = 0;
    for (int tid = 0; tid < hdr->n_targets; tid++) 
    {
        const char* chrom_name = hdr->target_name[tid]; // Name of the chromosome
        int chrom_num = chromToNumber(chrom_name);
        if (chrom_num) 
        {
            hts_pos_t chrom_len = hdr->target_len[tid];     // Length of the chromosome
            //fprintf(stderr, "Chromosome %s has length %" PRIi64 "\n", chrom_name, chrom_len);
            if (chrom_len > max_chrom_len) 
                max_chrom_len = chrom_len;
        }
    }
    //fprintf(stderr, "Maximum chromosome length: %" PRIi64 "\n", max_chrom_len);

    config.fai = fai_load(config.FastaName);
    if (config.fai == NULL) 
    {
        fprintf(stderr, "Index not found for %s, attempting to create it...\n", config.FastaName);
        if (fai_build(config.FastaName) != 0) 
        {
            fprintf(stderr, "Failed to create index for %s\n", config.FastaName);
            bam_hdr_destroy(config.hdr);
            hts_idx_destroy(config.bai);
            hts_close(fp);
            return 1;
        }
        
        config.fai = fai_load(config.FastaName);
        if (config.fai == NULL) 
        {
            fprintf(stderr, "Failed to load newly created index for %s\n", config.FastaName);
            bam_hdr_destroy(config.hdr);
            hts_idx_destroy(config.bai);
            hts_close(fp);
            return 1;
        }
    }

    hts_close(fp);  // Close the initial file handle
    
    // Calculate optimal thread count and buffer sizes based on system resources
    struct sysinfo si;
    if (sysinfo(&si) == 0) 
    {
        if (config.nThreads <= 1) 
        {
            config.nThreads = si.procs * 8 / 10;
            if (config.nThreads < 1) config.nThreads = 1;
        }
      
        config.chunkSize = (size_t)(max_chrom_len * sizeof(MethylData) / MAX(10, config.nThreads));
        if (config.chunkSize < 100000) config.chunkSize = 100000;
    }

    // Initialize the global chromosome buffer with a large capacity
    size_t initial_buffer_capacity = config.bufferSizePerThread * config.nThreads * 10; // Large enough for biggest chromosome
    initialize_global_chrom_buffer(initial_buffer_capacity);

    // Set total chromosomes for progress reporting
    total_chromosomes_global = config.hdr->n_targets;

    // Initialize work queue
    work_queue.units = malloc(1000 * sizeof(WorkUnit));
    if (!work_queue.units) 
    {
        return 1;
    }
    work_queue.capacity = 1000;
    work_queue.head = work_queue.tail = work_queue.current_count = 0;
    work_queue.current_chrom[0] = '\0';
    work_queue.is_current_chrom_done = 0;
    pthread_mutex_init(&work_queue.mutex, NULL);
    pthread_cond_init(&work_queue.not_empty, NULL);
    pthread_cond_init(&work_queue.not_full, NULL);
    pthread_cond_init(&work_queue.chrom_done, NULL);

    // Preload chromosome sequences into memory if possible
    typedef struct 
    {
        char* chrom_name;
        char* sequence;
        int seq_len;
    } ChromSequence;

    ChromSequence* chrom_sequences = malloc(hdr->n_targets * sizeof(ChromSequence));
    if (!chrom_sequences) 
    {
        return 1;
    }
    num_chromosomes = hdr->n_targets;
    for (int tid = 0; tid < hdr->n_targets; tid++) 
    {
        chrom_sequences[tid].chrom_name = strdup(hdr->target_name[tid]);
        chrom_sequences[tid].sequence = NULL;
        chrom_sequences[tid].seq_len = 0;
        
        // Check if chromosome exists in FASTA index
        if (check_chromosome_in_fai(config.fai, hdr->target_name[tid])) 
        {
            int len = faidx_seq_len(config.fai, hdr->target_name[tid]);
            if (len >= 0) 
            {
                chrom_sequences[tid].sequence = faidx_fetch_seq(config.fai, hdr->target_name[tid], 0, len, &chrom_sequences[tid].seq_len);
                if (chrom_sequences[tid].seq_len > 0) 
                {
                    //fprintf(stderr, "Preloaded sequence for chromosome %s, length: %d\n", hdr->target_name[tid], chrom_sequences[tid].seq_len);
                } 
                else 
                {
                    free(chrom_sequences[tid].sequence);
                    chrom_sequences[tid].sequence = NULL;
                }
            } 
            else 
            {
                // Try alternative name
                char alt_chrom[256];
                if (chromToNumber(hdr->target_name[tid]) != 0 && strncmp(hdr->target_name[tid], "chr", 3) == 0) 
                {
                    strncpy(alt_chrom, hdr->target_name[tid] + 3, sizeof(alt_chrom) - 1);
                    alt_chrom[sizeof(alt_chrom) - 1] = '\0';
                    len = faidx_seq_len(config.fai, alt_chrom);
                    if (len >= 0) 
                    {
                        chrom_sequences[tid].sequence = faidx_fetch_seq(config.fai, alt_chrom, 0, len, &chrom_sequences[tid].seq_len);
                        if (chrom_sequences[tid].seq_len > 0) 
                        {
                            //fprintf(stderr, "Preloaded sequence for chromosome %s using alt name %s, length: %d\n", hdr->target_name[tid], alt_chrom, chrom_sequences[tid].seq_len);
                        } 
                        else 
                        {
                            free(chrom_sequences[tid].sequence);
                            chrom_sequences[tid].sequence = NULL;
                        }
                    }
                } 
                else 
                {
                    snprintf(alt_chrom, sizeof(alt_chrom), "chr%s", hdr->target_name[tid]);
                    len = faidx_seq_len(config.fai, alt_chrom);
                    if (len >= 0) 
                    {
                        chrom_sequences[tid].sequence = faidx_fetch_seq(config.fai, alt_chrom, 0, len, &chrom_sequences[tid].seq_len);
                        if (chrom_sequences[tid].seq_len > 0) 
                        {
                            //fprintf(stderr, "Preloaded sequence for chromosome %s using alt name %s, length: %d\n", hdr->target_name[tid], alt_chrom, chrom_sequences[tid].seq_len);
                        } 
                        else 
                        {
                            free(chrom_sequences[tid].sequence);
                            chrom_sequences[tid].sequence = NULL;
                        }
                    }
                }
            }
        }
    }

    // Create work units for each chromosome
    for (int tid = 0; tid < config.hdr->n_targets; tid++) 
    {
        const char* raw_chrom = config.hdr->target_name[tid];
        int chrom_num = chromToNumber(raw_chrom);
        if (!chrom_num) continue;

        hts_pos_t chrom_len = config.hdr->target_len[tid];
        hts_pos_t pos = 0;
        
        while (pos < chrom_len) 
        {
            WorkUnit unit;
            strncpy(unit.chrom, raw_chrom, sizeof(unit.chrom) - 1);
            unit.chrom[sizeof(unit.chrom) - 1] = '\0';
            unit.tid = tid;
            unit.start_pos = pos;
            unit.end_pos = pos + config.chunkSize;
            if (unit.end_pos > chrom_len) unit.end_pos = chrom_len;
            
            // Add to work queue
            pthread_mutex_lock(&work_queue.mutex);
            while (work_queue.current_count >= work_queue.capacity) 
                pthread_cond_wait(&work_queue.not_full, &work_queue.mutex);
            work_queue.units[work_queue.tail] = unit;
            work_queue.tail = (work_queue.tail + 1) % work_queue.capacity;
            work_queue.current_count++;
            pthread_cond_signal(&work_queue.not_empty);
            pthread_mutex_unlock(&work_queue.mutex);
            
            pos = unit.end_pos;
        }
    }

    // Create and start worker threads
    pthread_t *threads = calloc(config.nThreads, sizeof(pthread_t));
    if (!threads) 
    {
        //fprintf(stderr, "Failed to allocate thread array\n");
        return 1;
    }

    for (int i = 0; i < config.nThreads; i++) 
    {
        if (pthread_create(&threads[i], NULL, extractCalls, &config) != 0) 
        {
            //fprintf(stderr, "Failed to create thread %d\n", i);
            return 1;
        }
    }

    // Wait for all threads to complete, no progress bar update here
    int total_errors = 0;
    for (int i = 0; i < config.nThreads; i++) 
    {
        void *result;
        pthread_join(threads[i], &result);
        if (result) total_errors++;
    }

    // Progress bar will be updated in signal_thread_completion, just print completion
    fprintf(stderr, "] Done\n");

    // Cleanup
    free(threads);
    free(work_queue.units);
    
    // Cleanup synchronization primitives
    pthread_mutex_destroy(&work_queue.mutex);
    pthread_cond_destroy(&work_queue.not_empty);
    pthread_cond_destroy(&work_queue.not_full);
    pthread_cond_destroy(&work_queue.chrom_done);
    pthread_mutex_destroy(&progress_mutex);
    
    // No need to destroy multiple chromosome buffer mutexes since we use a global buffer
    // The global buffer mutex is destroyed in free_global_chrom_buffer

    // Free the global chromosome buffer
    free_global_chrom_buffer();

    // Free preloaded chromosome sequences
    for (int tid = 0; tid < num_chromosomes; tid++) 
    {
        if (chrom_sequences[tid].chrom_name) free(chrom_sequences[tid].chrom_name);
        if (chrom_sequences[tid].sequence) free(chrom_sequences[tid].sequence);
    }
    free(chrom_sequences);
    chrom_sequences = NULL;
    num_chromosomes = 0;

    // Cleanup resources
    free(config.hdf5_output_dir);
    fai_destroy(config.fai);
    bam_hdr_destroy(config.hdr);
    hts_idx_destroy(config.bai);

    return total_errors ? 1 : 0;
}
