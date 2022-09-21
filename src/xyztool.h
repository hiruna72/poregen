/* @file xyztool.h
**
******************************************************************************/

#ifndef XYZTOOL_H
#define XYZTOOL_H

#include <stdint.h>
#include <slow5/slow5.h>

#define XYZTOOL_VERSION "0.1.0"

/*******************************************************
 * flags related to the user specified options (opt_t) *
 *******************************************************/

#define XYZTOOL_RNA 0x001 //if RNA or not
#define XYZTOOL_DTW 0x002 //if dtw-std
#define XYZTOOL_INV 0x004 //if set, reverse reference events instead of query events
#define XYZTOOL_SEC 0x008 //if secondaries are printed
#define XYZTOOL_REF 0x010 //map to the whole reference
#define XYZTOOL_END 0x020 //map the end of the query

#define WORK_STEAL 1 //simple work stealing enabled or not (no work stealing mean no load balancing)
#define STEAL_THRESH 1 //stealing threshold

/* user specified options */
typedef struct {

    uint32_t flag;              //flags
    int32_t batch_size;         //max reads loaded at once: K
    int64_t batch_size_bytes;   //max bytes loaded at once: B

    int32_t num_thread; //t
    int32_t debug_break;

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

#endif
