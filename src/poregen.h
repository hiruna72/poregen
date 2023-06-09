/* @file poregen.h
**
******************************************************************************/

#ifndef POREGEN_H
#define POREGEN_H

#include <stdint.h>
#include <slow5/slow5.h>
#include <string>
#include <vector>
#include <inttypes.h>

#define POREGEN_VERSION "0.1.0"

/*******************************************************
 * flags related to the user specified options (opt_t) *
 *******************************************************/

#define POREGEN_RNA 0x001 //if RNA or not
#define POREGEN_DTW 0x002 //if dtw-std
#define POREGEN_INV 0x004 //if set, reverse reference events instead of query events
#define POREGEN_SEC 0x008 //if secondaries are printed
#define POREGEN_REF 0x010 //map to the whole reference
#define POREGEN_END 0x020 //map the end of the query

#define WORK_STEAL 1 //simple work stealing enabled or not (no work stealing mean no load balancing)
#define STEAL_THRESH 1 //stealing threshold

#define TO_PICOAMPS(RAW_VAL,DIGITISATION,OFFSET,RANGE) (((RAW_VAL)+(OFFSET))*((RANGE)/(DIGITISATION)))
#define PROGRESS_BATCH_SIZE 10000
#define MAX_MIN_THRESHOLD 2.0
#define KMER_SAMPLE_LIMIT 100
#define KMERS_TO_DUMP_LIMIT 500
#define SIGNAL_PRINT_MARGIN 0
#define KMER_INDEX_START 1
#define NUM_DNA_BASES 4
#define DEFAULT_KMER_SIZE 9
#define MAXIMUM_MOVE_DURATION 70
#define MINIMUM_MOVE_DURATION 5
#define PA_MAX 180.0
#define PA_MIN 40.0
#define KMER_PICK_MARGIN 2

enum signal_scaling{ noscale, medmad_scale };

/* user specified options */
typedef struct {
    FILE *f_out;

    uint32_t flag;              //flags
    uint32_t flag_rna;          //RNA flag
    int32_t batch_size;         //max reads loaded at once: K
    int64_t batch_size_bytes;   //max bytes loaded at once: B

    int32_t num_thread; //t
    uint32_t kmer_size; //k
    uint32_t kmer_start_offset; //m
    uint32_t file_limit;
    uint32_t index_start;
    uint32_t index_end;
    int32_t debug_break;
    uint8_t use_paf_format;
    uint8_t delimit_files;
    char *arg_fname_out;

    //signal_related
    uint32_t max_dur; //max_dur
    uint32_t min_dur; //min_dur
    uint32_t sig_move_offset; //m
    uint32_t signal_print_margin;
    uint32_t sample_limit;
    enum signal_scaling signal_scale;
    double pa_max;
    double pa_min;

    uint32_t kmer_pick_margin;

} opt_t;


/* a batch of read data (dynamic data based on the reads) */
typedef struct {

    int32_t n_rec;
    int32_t capacity_rec;

    char **mem_records;
    size_t *mem_bytes;

    slow5_rec_t **slow5_rec;

    double *means;

    //stats
    int64_t sum_bytes;
    int64_t total_reads; //total number mapped entries in the bam file (after filtering based on flags, mapq etc)


} db_t;



/* core data structure (mostly static data throughout the program lifetime) */
typedef struct {

    //slow5
    slow5_file_t *sp;

    // options
    opt_t opt;

    //realtime0
    double realtime0;

    double load_db_time;
    double process_db_time;
    double output_time;

    //stats //set by output_db
    int64_t sum_bytes;
    int64_t total_reads; //total number mapped entries in the bam file (after filtering based on flags, mapq etc)

} core_t;


/* argument wrapper for the multithreaded framework used for data processing */
typedef struct {
    core_t* core;
    db_t* db;
    int32_t starti;
    int32_t endi;
    void (*func)(core_t*,db_t*,int);
    int32_t thread_index;
#ifdef WORK_STEAL
    void *all_pthread_args;
#endif
#ifdef HAVE_CUDA
    int32_t *ultra_long_reads; //reads that are assigned to the CPU due to the unsuitability to process on the GPU
    double ret1;    //return value
#endif
} pthread_arg_t;

/* return status by the load_db - used for termination when all the data is processed */
typedef struct {
    int32_t num_reads;
    int64_t num_bytes;
} ret_status_t;

/******************************************
 * function prototype for major functions *
 ******************************************/

/* initialise user specified options */
void init_opt(opt_t* opt);

/* initialise the core data structure */
core_t* init_core(char *slow5file, opt_t opt, double realtime0);

/* initialise a data batch */
db_t* init_db(core_t* core);

/* load a data batch from disk */
ret_status_t load_db(core_t* dg, db_t* db);

void work_per_single_read(core_t* core,db_t* db, int32_t i);
/* process all reads in the given batch db */
void work_db(core_t* core, db_t* db, void (*func)(core_t*,db_t*,int));

/* process a data batch */
void process_db(core_t* core, db_t* db);

/* align a single read specified by index i*/
void process_single(core_t* core, db_t* db, int32_t i);

/* write the output for a processed data batch */
void output_db(core_t* core, db_t* db);

/* partially free a data batch - only the read dependent allocations are freed */
void free_db_tmp(db_t* db);

/* completely free a data batch */
void free_db(db_t* db);

/* free the core data structure */
void free_core(core_t* core,opt_t opt);

// introduce new functions here
void generate_kmers(char set[], std::string prefix, int num_bases, int kmer_length, std::vector<std::string>& kmers);
#endif
