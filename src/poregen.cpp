/* @file poregen.c
**
** @@
******************************************************************************/
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "poregen.h"
#include "misc.h"
#include "error.h"

#include <slow5/slow5.h>

#include <sys/wait.h>
#include <unistd.h>

enum poregen_log_level_opt _log_level = LOG_INFO;

/* initialise the core data structure */
core_t* init_core(char *slow5file, opt_t opt,double realtime0) {

    core_t* core = (core_t*)malloc(sizeof(core_t));
    MALLOC_CHK(core);

    core->sp = slow5_open(slow5file,"r");
    if (core->sp == NULL) {
        VERBOSE("Error opening SLOW5 file %s\n",slow5file);
        exit(EXIT_FAILURE);
    }


    core->opt = opt;

    //realtime0
    core->realtime0=realtime0;

    core->load_db_time=0;
    core->process_db_time=0;
    core->output_time=0;

    core->sum_bytes=0;
    core->total_reads=0; //total number mapped entries in the bam file (after filtering based on flags, mapq etc)

    return core;
}

/* free the core data structure */
void free_core(core_t* core,opt_t opt) {
    slow5_close(core->sp);
    free(core);
}

/* initialise a data batch */
db_t* init_db(core_t* core) {
    db_t* db = (db_t*)(malloc(sizeof(db_t)));
    MALLOC_CHK(db);

    db->capacity_rec = core->opt.batch_size;
    db->n_rec = 0;

    db->mem_records = (char**)(calloc(db->capacity_rec,sizeof(char*)));
    MALLOC_CHK(db->mem_records);
    db->mem_bytes = (size_t*)(calloc(db->capacity_rec,sizeof(size_t)));
    MALLOC_CHK(db->mem_bytes);

    db->slow5_rec = (slow5_rec_t**)calloc(db->capacity_rec,sizeof(slow5_rec_t*));
    MALLOC_CHK(db->slow5_rec);

    db->means = (double*)calloc(db->capacity_rec,sizeof(double));
    MALLOC_CHK(db->means);

    db->total_reads=0;
    db->sum_bytes=0;


    return db;
}

/* load a data batch from disk */
ret_status_t load_db(core_t* core, db_t* db) {

    double load_start = realtime();

    db->n_rec = 0;
    db->sum_bytes = 0;
    db->total_reads = 0;

    ret_status_t status={0,0};
    int32_t i = 0;
    while (db->n_rec < db->capacity_rec && db->sum_bytes<core->opt.batch_size_bytes) {
        i=db->n_rec;

        if (slow5_get_next_bytes(&db->mem_records[i],&db->mem_bytes[i],core->sp) < 0){
            if (slow5_errno != SLOW5_ERR_EOF) {
                ERROR("Error reading from SLOW5 file %d", slow5_errno);
                exit(EXIT_FAILURE);
            } else {
                break;
            }
        } else {
            db->n_rec++;
            db->total_reads++; // candidate read
            db->sum_bytes += db->mem_bytes[i];
        }
    }

    status.num_reads=db->n_rec;
    status.num_bytes=db->sum_bytes;

    double load_end = realtime();
    core->load_db_time += (load_end-load_start);

    return status;
}


void parse_single(core_t* core,db_t* db, int32_t i){

    assert(db->mem_bytes[i]>0);
    assert(db->mem_records[i]!=NULL);
    //db->slow5_rec[i]=NULL;
    int ret=slow5_decode(&db->mem_records[i], &db->mem_bytes[i], &db->slow5_rec[i], core->sp);
    if(ret<0){
        ERROR("Error parsing the record %d",i);
        exit(EXIT_FAILURE);
    }

}

#define TO_PICOAMPS(RAW_VAL,DIGITISATION,OFFSET,RANGE) (((RAW_VAL)+(OFFSET))*((RANGE)/(DIGITISATION)))

void mean_single(core_t* core,db_t* db, int32_t i){

    slow5_rec_t* rec = db->slow5_rec[i];
    uint64_t len_raw_signal = rec->len_raw_signal;

    if(len_raw_signal>0){
        double sum = 0;
        for(uint64_t i=0;i<len_raw_signal;i++){
            double pA = TO_PICOAMPS(rec->raw_signal[i],rec->digitisation,rec->offset,rec->range);
            sum += pA;
        }
        double mean = sum/len_raw_signal;
        db->means[i]=mean;
    }

}

void work_per_single_read(core_t* core,db_t* db, int32_t i){
    parse_single(core,db,i);
    mean_single(core,db,i);
}

void process_db(core_t* core,db_t* db){
    double proc_start = realtime();
    work_db(core, db, work_per_single_read);
    double proc_end = realtime();
    core->process_db_time += (proc_end-proc_start);
}

/* write the output for a processed data batch */
void output_db(core_t* core, db_t* db) {

    double output_start = realtime();

    int32_t i = 0;
    for (i = 0; i < db->n_rec; i++) {
        if(db->slow5_rec[i]->len_raw_signal>0){
            printf("%s\t%f\n",db->slow5_rec[i]->read_id,db->means[i]);
        }
    }

    core->sum_bytes += db->sum_bytes;
    core->total_reads += db->total_reads;

    //core->read_index = core->read_index + db->n_rec;
    double output_end = realtime();
    core->output_time += (output_end-output_start);

}

/* partially free a data batch - only the read dependent allocations are freed */
void free_db_tmp(db_t* db) {
    int32_t i = 0;
    for (i = 0; i < db->n_rec; ++i) {
        free(db->mem_records[i]);
    }
}

/* completely free a data batch */
void free_db(db_t* db) {

    int32_t i = 0;
    for (i = 0; i < db->capacity_rec; ++i) {
        slow5_rec_free(db->slow5_rec[i]);
    }
    free(db->slow5_rec);
    free(db->mem_records);
    free(db->mem_bytes);
    free(db->means);
    free(db);
}

/* initialise user specified options */
void init_opt(opt_t* opt) {
    memset(opt, 0, sizeof(opt_t));
    opt->batch_size = 512;
    opt->batch_size_bytes = 20*1000*1000;
    opt->num_thread = 8;
    opt->kmer_size = DEFAULT_KMER_SIZE;
    opt->sig_move_offset = 0;
    opt->kmer_start_offset = 0;
    opt->signal_print_margin = SIGNAL_PRINT_MARGIN;
    opt->sample_limit = KMER_SAMPLE_LIMIT;
    opt->file_limit = KMERS_TO_DUMP_LIMIT;
    opt->index_start = KMER_INDEX_START;
    opt->index_end = opt->index_start + KMERS_TO_DUMP_LIMIT - 1;
    opt->signal_scale = medmad_scale;
    opt->debug_break=-1;
    opt->use_paf_format = 0;
    opt->delimit_files = 0;
    opt->arg_fname_out = NULL;
    opt->f_out = stdout;

    opt->max_dur = MAXIMUM_MOVE_DURATION;
    opt->min_dur = MINIMUM_MOVE_DURATION;
    opt->pa_max = PA_MAX;
    opt->pa_min = PA_MIN;

    opt->flag_rna = 0;
    opt->kmer_pick_margin = KMER_PICK_MARGIN;

}


enum poregen_log_level_opt get_log_level(){
    return _log_level;
}

void set_log_level(enum poregen_log_level_opt level){
    _log_level = level;
}

void generate_kmers(char set[], std::string prefix, int num_bases, int kmer_length, std::vector<std::string>& kmers){
    // Base case: k is 0,
    // print prefix
    if (kmer_length == 0){
        //        std::cout << (prefix) << std::endl;
        kmers.push_back(prefix);
        return;
    }
    // One by one add all characters
    // from set and recursively
    // call for k equals to k-1
    for (int i = 0; i < num_bases; i++){
        std::string newPrefix;
        // Next character of input added
        newPrefix = prefix + set[i];
        // k is decreased, because
        // we have added a new character
        generate_kmers(set, newPrefix, num_bases, kmer_length - 1, kmers);
    }
}


