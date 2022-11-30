/* @file  gmove.cpp
**
** Hiruna Samarakoon - hiruna72@gmail.com
******************************************************************************/
#include "poregen.h"
#include "error.h"
#include "misc.h"
#include <assert.h>
#include <getopt.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <vector>
#include <string>
#include <dirent.h>
#include <map>

#ifndef DISABLE_KSORT
#include "ksort.h"
KSORT_INIT_GENERIC(float)
KSORT_INIT_GENERIC(double)
KSORT_INIT_GENERIC(int16_t)
KSORT_INIT_GENERIC(int)
#endif

#define TO_PICOAMPS(RAW_VAL,DIGITISATION,OFFSET,RANGE) (((RAW_VAL)+(OFFSET))*((RANGE)/(DIGITISATION)))
#define PROGRESS_BATCH_SIZE 10000

static struct option long_options[] = {
    {"kmer_size", required_argument, 0, 'k'},      //0 kmer_size [6]
    {"scaling", required_argument, NULL, 0},       //1 scaling 1-medmad
    {"verbose", required_argument, 0, 'v'},        //3 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //4
    {"version", no_argument, 0, 'V'},              //5
    {"debug-break",required_argument, 0, 0},       //7 break after processing the first batch (used for debugging)
    {0, 0, 0, 0}};


static inline void print_help_msg(FILE *fp_help, opt_t opt){
    fprintf(fp_help,"Usage: poregen gmove reads.blow5 move_table output_dir\n");
    fprintf(fp_help,"\nbasic options:\n");
    fprintf(fp_help,"   -k INT                     kmer_size [%d]\n",opt.kmer_size);
    fprintf(fp_help,"   --scaling                  scaling [%d] (0-no scaling, 1-medmad)\n",opt.signal_scale);
    fprintf(fp_help,"   -h                         help\n");
    fprintf(fp_help,"   --verbose INT              verbosity level [%d]\n",(int)get_log_level());
    fprintf(fp_help,"   --version                  print version\n");
}

// from nanopolish
std::vector< std::string > list_directory(const std::string& file_name){
    std::vector< std::string > res;
    DIR* dir;
    struct dirent *ent;

    dir = opendir(file_name.c_str());
    if(not dir) {
    return res;
    }
    while((ent = readdir(dir))) {
    res.push_back(ent->d_name);
    }
    closedir(dir);
    return res;
}

//return 0; created dir
//return -1; dir exists
//return -2; could not create dir
int create_dir(const char *dir_name) {
    struct stat st = {0};
    if (stat(dir_name, &st) == -1) {
        int ret_mkdir = mkdir(dir_name, 0700);
        if(ret_mkdir == -1){
            return -2;
        }
    }else{
        std::vector< std::string > dir_list = list_directory(dir_name);
        if(dir_list.size()>2){
            return -1;
        }
    }
    return 0;
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

static inline double calc_median(const double* x, int n) {

    double *copy = (double *)malloc(n * sizeof(double));
    memcpy(copy, x, n * sizeof(double));
    double m = ks_ksmall_double(n, copy, n / 2);
    free(copy);
    return m;

}


/** Median Absolute Deviation of an array
 *
 *	@param x   An array to calculate the MAD of
 *	@param n   Length of array
 *	@param med Median of the array.	 If NAN then median is calculated.
 *
 *	@return MAD of array on success, NAN otherwise.
 **/
static double calc_madf(const double* x, size_t n, const double* med) {
    const double mad_scaling_factor = 1.4826;
    if (NULL == x) {
        return NAN;
    }
    if (1 == n) {
        return 0.0;
    }

    double* absdiff = (double*)malloc(n * sizeof(double));
    if (NULL == absdiff) {
        return NAN;
    }

    const double _med = (NULL == med) ? calc_median(x, n) : *med;

    for (size_t i = 0; i < n; i++) {
        absdiff[i] = fabs(x[i] - _med);
    }

    const double mad = calc_median(absdiff, n);
    free(absdiff);
    return mad * mad_scaling_factor;
}


void close_kmer_files(std::map<std::string, FILE *> *kmer_file_pointer_array_ptr, int num_kmers, std::vector<std::string>& kmers) {
    std::map<std::string, FILE *> kmer_file_pointer_array = *kmer_file_pointer_array_ptr;
    for(int i=0;i<num_kmers;i++){
        fclose(kmer_file_pointer_array[kmers[i]]);
    }
}


void delimit_kmer_files(std::map<std::string, FILE *> *kmer_file_pointer_array_ptr, int num_kmers, std::vector<std::string>& kmers) {
    std::map<std::string, FILE *> kmer_file_pointer_array = *kmer_file_pointer_array_ptr;
    for(int i=0;i<num_kmers;i++){
        fprintf(kmer_file_pointer_array[kmers[i]],":");
    }
}

int gmove(int argc, char* argv[]) {

//    fprintf(stderr,"%d", get_log_level());

    const char* optstring = "k:";

    int longindex = 0;
    int32_t c = -1;

    char *slow5file = NULL;
    char *move_table = NULL;
    char *output_dir = NULL;
    int signal_scale = 0;

    char dna_set[] = {'A', 'C', 'G', 'T'};
    int num_dna_bases = 4;
    int kmer_start_offset = 0;
    FILE *fp_help = stderr;

    opt_t opt;
    init_opt(&opt); //initialise options to defaults

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {

        if (c == 'k') {
            opt.kmer_size = atoi(optarg);
            if (opt.kmer_size < 1) {
                ERROR("Kmer length should larger than 0. You entered %d", opt.kmer_size);
                exit(EXIT_FAILURE);
            }
        } else if (c=='v'){
            int v = atoi(optarg);
            set_log_level((enum poregen_log_level_opt)v);
        } else if (c=='V'){
            fprintf(stdout,"gmove %s\n",POREGEN_VERSION);
            exit(EXIT_SUCCESS);
        } else if (c=='h'){
            fp_help = stdout;
            fp_help = stdout;
        } else if(c == 0 && longindex == 1){ //debug break
            signal_scale = atoi(optarg);
        }
        else if(c == 0 && longindex == 15){ //debug break
            opt.debug_break = atoi(optarg);
        }
    }

    // No arguments given
    if (argc - optind != 3 || fp_help == stdout) {
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }
    slow5file = argv[optind++];
    if (slow5file == NULL) {
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }
    move_table = argv[optind++];
    if (move_table == NULL) {
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }
    output_dir = argv[optind++];
    if (output_dir == NULL) {
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }


    int move_start_offset = opt.kmer_size - 1;
    fprintf(stderr,"slow5_file_path: %s\n", slow5file);
    fprintf(stderr,"guppy_sam_output_file: %s\n", move_table);
    fprintf(stderr,"kmer_output_dir: %s\n", output_dir);
    fprintf(stderr,"kmer_size: %d\n", opt.kmer_size);
    fprintf(stderr,"move_start_offset: %d\n", move_start_offset);
    if(signal_scale == 0){
        opt.signal_scale = noscale;
        fprintf(stderr,"scaling: %s\n", "no scale");
    }else if(signal_scale == 1){
        opt.signal_scale = medmad_scale;
        fprintf(stderr,"scaling: %s\n", "medmad scale");
    }else {
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    int ret_create_dir = create_dir(output_dir);
    if(ret_create_dir == -1){
        fprintf(stderr,"Output directory %s is not empty. Please remove it or specify another directory.", output_dir);
        return EXIT_FAILURE;
    }
    if(ret_create_dir == -2){
        fprintf(stderr,"Could not create the output dir %s.", output_dir);
        return EXIT_FAILURE;
    }
    std::string kmer_dump = std::string(std::string(output_dir) + "/kmer_dump");
    ret_create_dir = create_dir(kmer_dump.c_str());
    if(ret_create_dir == -1){
        fprintf(stderr,"Output directory %s is not empty. Please remove it or specify another directory.", kmer_dump.c_str());
        return EXIT_FAILURE;
    }
    if(ret_create_dir == -2){
        fprintf(stderr,"Could not create the output dir %s.", kmer_dump.c_str());
        return EXIT_FAILURE;
    }

    std::vector<std::string> kmers;
    generate_kmers(dna_set, "", num_dna_bases, opt.kmer_size, kmers);
    int num_kmers = kmers.size();
    fprintf(stderr,"num_kmers: %d\n", num_kmers);

    std::map<std::string,FILE*> kmer_file_pointer_array;
    for(int i=0;i<num_kmers;i++){
        FILE * kmer_output = NULL;
        std::string output_path = kmer_dump + "/" + kmers[i];
        kmer_output = fopen(output_path.c_str(), "w");
        if (kmer_output == NULL){
            fprintf(stderr,"Error in opening %dth kmer-file %s\n", i, output_path.c_str());
            exit(EXIT_FAILURE);
        }
        fprintf(kmer_output,"%.8f,", 0.0);
        kmer_file_pointer_array[kmers[i]] = kmer_output;
    }
    std::map<std::string,uint64_t> kmer_frequency_map;
    for(int i=0;i<num_kmers;i++){
        kmer_frequency_map[kmers[i]] = 0;
    }


    slow5_file_t *sp = slow5_open(slow5file,"r");
    if(sp==NULL){
        fprintf(stderr,"Error in opening file %s\n", slow5file);
        exit(EXIT_FAILURE);
    }
    slow5_rec_t *rec = NULL;
    int ret=0;

    ret = slow5_idx_load(sp);
    if(ret<0){
        fprintf(stderr,"Error in loading index\n");
        exit(EXIT_FAILURE);
    }

    FILE * basecall_fp = NULL;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    basecall_fp = fopen(move_table, "r");
    if (basecall_fp == NULL){
        fprintf(stderr,"Error in opening file %s\n", move_table);
        exit(EXIT_FAILURE);
    }
    int count_reads = 0;

    while ((read = getline(&line, &len, basecall_fp)) != -1) {
        int fastq_len;
        uint64_t signal_len;
        int trim_offset;
        int stride;
        size_t move_len;
        std::string read_id;
        std::string fastq_seq;
        std::string move_seq;

        read_id = strtok(line, "\t");
        fastq_len = atoi(strtok(NULL, "\t"));
        fastq_seq = strtok(NULL, "\t");
        stride = atoi(strtok(NULL, "\t"));
        move_seq = strtok(NULL, "\t");
        signal_len = std::stoull(strtok(NULL, "\t"));
        trim_offset = atoi(strtok(NULL, "\t"));
        move_len = move_seq.size();

        VERBOSE("%s\n",read_id.c_str());
        VERBOSE("%s\t%d\t%" PRIu64 "\t%d\t%d\n",read_id.c_str(),fastq_len,signal_len,trim_offset,stride);

        ret = slow5_get(read_id.c_str(), &rec, sp);
        if(ret < 0){
            fprintf(stderr,"Error in when fetching the read\n");
            exit(EXIT_FAILURE);
        }
        uint64_t len_raw_signal = rec->len_raw_signal;
        std::vector<double> raw_signal(len_raw_signal);
        assert(len_raw_signal==signal_len);
        assert((uint64_t)trim_offset < len_raw_signal);
        len_raw_signal = len_raw_signal - trim_offset;
        for(uint64_t i=0;i<len_raw_signal;i++){
            double pA = TO_PICOAMPS(rec->raw_signal[i+trim_offset],rec->digitisation,rec->offset,rec->range);
            raw_signal[i] = pA;
//                printf("%f ",pA);
        }
        //calculate medmad
        if(signal_scale == medmad_scale){
            double read_median = calc_median(raw_signal.data(), len_raw_signal);
            if(read_median == NAN){
                continue;
            }
            double read_mad = calc_madf(raw_signal.data(),len_raw_signal, &read_median);
            if(read_mad == NAN){
                continue;
            }
            read_mad = (read_mad > 1.0) ? read_mad : 1.0;

            for(uint64_t i=0;i<len_raw_signal;i++){
                raw_signal[i] = (raw_signal[i]-read_median) / read_mad;
            }
        }

        if(fastq_len < 10){
            fprintf(stderr,"\t%d", fastq_len);
            continue;
        }
        int move_count = 0;
        size_t move_idx=0;
        size_t start_move_idx = 0;
        while(move_count <= move_start_offset){
            if(move_seq[move_idx]=='1'){
                move_count++;
                start_move_idx = move_idx;
            }
            move_idx++;
        }
        move_idx = start_move_idx + 1;
        size_t seq_start = kmer_start_offset;
        for(;move_idx<=move_len;move_idx++){
            if(move_seq[move_idx]=='1'){
                std::string kmer = fastq_seq.substr(seq_start, opt.kmer_size);
                fprintf(stdout,"%s\t", kmer.c_str());
//                int raw_start_local = trim_offset + (start_move_idx*stride);
//                int raw_end_local = trim_offset + (move_idx*stride);
                int raw_start_local = start_move_idx*stride;
                int raw_end_local = move_idx*stride;

                fprintf(stdout,"%d\t%d\n", raw_start_local, raw_end_local);

                std::vector<double> raw_signal_local(raw_signal.begin()+raw_start_local,raw_signal.begin()+raw_end_local);

                size_t raw_signal_kmer_length = raw_signal_local.size();
                size_t k;
                for(k=0;k<raw_signal_kmer_length-1;k++){
                    fprintf(kmer_file_pointer_array[kmer],"%.8f,", raw_signal_local[k]);
                    fprintf(stdout, "%.8f,", raw_signal_local[k]);
                }
                fprintf(kmer_file_pointer_array[kmer],"%.8f;", raw_signal_local[k]);
                fprintf(stdout, "%.8f\n", raw_signal_local[k]);
                if(seq_start+opt.kmer_size == fastq_seq.size()){
                    break;
                }
                seq_start++;
                start_move_idx = move_idx;
                kmer_frequency_map[kmer] = kmer_frequency_map[kmer] + 1;
            }
        }
        delimit_kmer_files(&kmer_file_pointer_array, num_kmers, kmers);
        if(count_reads == PROGRESS_BATCH_SIZE){
            fprintf(stderr,"*");
            count_reads = 0;
        }
        count_reads++;
    }
    fclose(basecall_fp);
    if (line){
        free(line);
    }

    std::string kmer_freq_file_name = std::string(std::string(output_dir)+"/kmer_freq.txt");
    FILE* kmer_freq_file = fopen(kmer_freq_file_name.c_str(), "w");
    if (kmer_freq_file == NULL){
        fprintf(stderr,"Error in opening %s\n", kmer_freq_file_name.c_str());
        exit(EXIT_FAILURE);
    }
    for(int i=0;i<num_kmers;i++){
        fprintf(kmer_freq_file, "%s\t%" PRIu64 "\n", kmers[i].c_str(), kmer_frequency_map[kmers[i]]);
    }

    return 0;
}
