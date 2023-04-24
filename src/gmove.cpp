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
#include <algorithm>
#include <htslib/faidx.h>
#include <htslib/sam.h>

#ifndef DISABLE_KSORT
#include "ksort.h"
KSORT_INIT_GENERIC(float)
KSORT_INIT_GENERIC(double)
KSORT_INIT_GENERIC(int16_t)
KSORT_INIT_GENERIC(int)
#endif

#define MAX_LEN_KMER 2000
#define EXPECTED_STRIDE 5

typedef struct{
    char *rid;
    int32_t qlen;
    int32_t query_start;
    int32_t query_end;
    int8_t strand;
    char *tid;
    int32_t tlen;
    int32_t target_start;
    int32_t target_end;
    uint8_t mapq;
    char *ss;
}paf_rec_t;

static struct option long_options[] = {
    {"kmer_size", required_argument, 0, 'k'},           //0 kmer_size [6]
    {"sig_move_offset", required_argument, 0, 'm'},   //1 sig_move_offset [5]
    {"kmer_start_offset", required_argument, 0, 's'},   //2 kmer_start_offset [0]
    {"scaling", required_argument, NULL, 0},       //3 scaling 1-medmad
    {"margin", required_argument, NULL, 0},        //4
    {"sample_limit", required_argument, NULL, 0},  //5
    {"file_limit", required_argument, NULL, 0},    //6
    {"kmer_file", required_argument, NULL, 0},     //7
    {"index_start", required_argument, NULL, 0},   //8
    {"index_end", required_argument, NULL, 0},     //9
    {"fastq", required_argument, NULL, 0},     //10
    {"", no_argument, 0, 'd'},                 //11 delimit output files
    {"max_dur", required_argument, 0, 0},           //12
    {"pa_min", required_argument, 0, 0},           //13
    {"pa_max", required_argument, 0, 0},           //14
    {"verbose", required_argument, 0, 'v'},        //15 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //16
    {"version", no_argument, 0, 'V'},              //17
    {"debug-break",required_argument, 0, 0},       //18 break after processing the first batch (used for debugging)
    {0, 0, 0, 0}};

void process_move_table_file(char *move_table, std::map<std::string,FILE*> &kmer_file_pointer_array, slow5_file_t **sp_ptr, opt_t *opt_ptr, std::map<std::string,uint64_t> &kmer_frequency_map, std::vector<std::string> &kmers);
void process_move_table_bam(char *move_table, std::map<std::string,FILE*> &kmer_file_pointer_array, slow5_file_t **sp_ptr, opt_t *opt_ptr, std::map<std::string,uint64_t> &kmer_frequency_map, std::vector<std::string> &kmers);
void process_move_table_paf(char *move_table, std::map<std::string,FILE*> &kmer_file_pointer_array, slow5_file_t **sp_ptr, opt_t *opt_ptr, std::map<std::string,uint64_t> &kmer_frequency_map, std::vector<std::string> &kmers, char *input_fastq_file);
void free_paf_rec(paf_rec_t *paf);
paf_rec_t *parse_paf_rec(char *buffer);

static inline void print_help_msg(FILE *fp_help, opt_t opt){
    fprintf(fp_help,"Usage: poregen gmove reads.blow5 move_table output_dir\n");
    fprintf(fp_help,"\nbasic options:\n");
    fprintf(fp_help,"   -k INT                     kmer_size [%d]\n",opt.kmer_size);
    fprintf(fp_help,"   -m INT                     move start offset [%d]\n",opt.sig_move_offset);
    fprintf(fp_help,"   -s INT                     kmer start offset [%d]\n",opt.kmer_start_offset);
    fprintf(fp_help,"   --scaling INT              scaling [%d] (0-no scaling, 1-medmad)\n",opt.signal_scale);
    fprintf(fp_help,"   --margin INT               signal print margin on both sides of the sub signal[%u] \n",opt.signal_print_margin);
    fprintf(fp_help,"   --sample_limit INT         maximum number of instances to output for a kmer [%u] \n",opt.sample_limit);
    fprintf(fp_help,"   --file_limit INT           maximum number of kmer files to output [%u] \n",opt.file_limit);
    fprintf(fp_help,"   --kmer_file FILE           kmer file (optional) \n");
    fprintf(fp_help,"   --index_start INT          1-based closed interval index of start kmer [%u] \n",opt.file_limit);
    fprintf(fp_help,"   --index_end INT            1-based closed interval index of end kmer [%u] \n",opt.file_limit);
    fprintf(fp_help,"   --fastq FILE               fastq file (optional - should be provided with .paf) \n");
    fprintf(fp_help,"   -d                         delimit output files per read\n");
    fprintf(fp_help,"   --max_dur                  maximum move duration allowed for samples [%d]\n",opt.max_dur);
    fprintf(fp_help,"   --pa_min                   minimum pA level a sampling signal should have [%.3f]\n",opt.pa_min);
    fprintf(fp_help,"   --pa_madx                  maximum pA level a sampling signal should have [%.3f]\n",opt.pa_max);
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


void close_kmer_files(std::map<std::string, FILE *> *kmer_file_pointer_array_ptr, std::vector<std::string>& kmers, opt_t* opt) {
    std::map<std::string, FILE *> kmer_file_pointer_array = *kmer_file_pointer_array_ptr;
    for(uint32_t i=opt->index_start-1;i<opt->index_end;i++){
        if(kmer_file_pointer_array[kmers[i]]){
            fclose(kmer_file_pointer_array[kmers[i]]);
        }
    }
}

void delimit_kmer_files(std::map<std::string, FILE *> *kmer_file_pointer_array_ptr, std::vector<std::string>& kmers, opt_t* opt) {
    std::map<std::string, FILE *> kmer_file_pointer_array = *kmer_file_pointer_array_ptr;
    for(uint32_t i=opt->index_start-1;i<opt->index_end;i++){
        if(kmer_file_pointer_array[kmers[i]]){
            fprintf(kmer_file_pointer_array[kmers[i]],":");
        }
    }
}

int gmove(int argc, char* argv[]) {

//    fprintf(stderr,"%d", get_log_level());

    const char* optstring = "k:m:s:d";

    int longindex = 0;
    int32_t c = -1;
    int flag_kmer_file_avail = 0;
    int flag_fastq_file_avail = 0;

    char *slow5file = NULL;
    char *move_table = NULL;
    char *output_dir = NULL;
    char *input_kmer_file = NULL;
    char *input_fastq_file = NULL;
    int signal_scale = 0;

    char dna_set[] = {'A', 'C', 'G', 'T'};

    FILE *fp_help = stderr;

    opt_t opt;
    init_opt(&opt); //initialise options to defaults

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {

        if (c == 'k') {
            if (atoi(optarg) < 1) {
                ERROR("Kmer length should be larger than 0. You entered %d", atoi(optarg));
                exit(EXIT_FAILURE);
            }
            opt.kmer_size = atoi(optarg);
        } else if (c == 'm') {
            if (atoi(optarg) < 0) {
                ERROR("signal move offset value must not be less than zero. You entered %d", atoi(optarg));
                exit(EXIT_FAILURE);
            }
            opt.sig_move_offset = atoi(optarg);
        } else if (c == 's') {
            if (atoi(optarg) < 1) {
                ERROR("Kmer offset should be larger than 0. You entered %d", atoi(optarg));
                exit(EXIT_FAILURE);
            }
            opt.kmer_start_offset = atoi(optarg);
        } else if (c == 'd') {
            opt.delimit_files = 1;
        } else if (c=='v'){
            int v = atoi(optarg);
            set_log_level((enum poregen_log_level_opt)v);
        } else if (c=='V'){
            fprintf(stdout,"gmove %s\n",POREGEN_VERSION);
            exit(EXIT_SUCCESS);
        } else if (c=='h'){
            fp_help = stdout;
            fp_help = stdout;
        } else if(c == 0 && longindex == 3){ //scaling
            signal_scale = atoi(optarg);
        } else if (c == 0 && longindex == 4){ //margin
            if(atoi(optarg) < 0){
                ERROR("Signal print margin should be non negative. You entered %d", atoi(optarg));
                exit(EXIT_FAILURE);
            }
            opt.signal_print_margin = atoi(optarg);
        } else if (c == 0 && longindex == 5){ //kmer_limit
            if(atoi(optarg) < 0){
                ERROR("Maximum number of instances to output for a kmer should be non negative. You entered %d", atoi(optarg));
                exit(EXIT_FAILURE);
            }
            opt.sample_limit = atoi(optarg);
        } else if (c == 0 && longindex == 6){ //dump_limit
            if(atoi(optarg) < 0){
                ERROR("Maximum number of kmer files to output should be non negative. You entered %d", atoi(optarg));
                exit(EXIT_FAILURE);
            }
            opt.file_limit = atoi(optarg);
            opt.index_end = opt.index_start + opt.file_limit - 1;
        } else if (c == 0 && longindex == 7){
            input_kmer_file = optarg;
            flag_kmer_file_avail = 1;
        } else if (c == 0 && longindex == 8){
            if(atoi(optarg) < 1){
                ERROR("kmer index start should be a positive number. You entered %d", atoi(optarg));
                exit(EXIT_FAILURE);
            }
            opt.index_start = atoi(optarg);
            opt.file_limit = opt.index_end - opt.index_start + 1;
        } else if (c == 0 && longindex == 9){
            if(atoi(optarg) < 1){
                ERROR("kmer index end should be a positive number. You entered %d", atoi(optarg));
                exit(EXIT_FAILURE);
            }
            opt.index_end = atoi(optarg);
            opt.file_limit = opt.index_end - opt.index_start + 1;
        } else if (c == 0 && longindex == 10){
            input_fastq_file = optarg;
            flag_fastq_file_avail = 1;
        } else if (c == 0 && longindex == 12){
            opt.max_dur = atoi(optarg);
        } else if (c == 0 && longindex == 13){
            opt.pa_min = atof(optarg);
        } else if (c == 0 && longindex == 14){
            opt.pa_max = atof(optarg);
        } else if (c == 0 && longindex == 18){ //debug break
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

    if(opt.kmer_size < 1){
        ERROR("kmer length must be a positive integer%s", "")
        return -1;
    }
    if(opt.sig_move_offset < 0){
        ERROR("signal move offset value must not be less than zero%s", "")
        return -1;
    }
//    if(opt.kmer_size <= opt.sig_move_offset){
//        ERROR("signal move offset value must less than the kmer length%s", "")
//        return -1;
//    }

    int ret_create_dir = create_dir(output_dir);
    if(ret_create_dir == -1){
        fprintf(stderr,"Output directory %s is not empty. Please remove it or specify another directory.", output_dir);
        return EXIT_FAILURE;
    }
    if(ret_create_dir == -2){
        fprintf(stderr,"Could not create the output dir %s.", output_dir);
        return EXIT_FAILURE;
    }
    std::string kmer_dump = std::string(std::string(output_dir) + "/dump");
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
    if(flag_kmer_file_avail == 1){
        FILE* input_kmer_file_ptr = fopen(input_kmer_file, "r");
        if (input_kmer_file_ptr == NULL){
            fprintf(stderr,"Error in opening file %s\n", input_kmer_file);
            exit(EXIT_FAILURE);
        }

        char * line = NULL;
        size_t len = 0;
        ssize_t read;

        while ((read = getline(&line, &len, input_kmer_file_ptr)) != -1) {
            line[read-1] = '\0';
//            fprintf(stderr,"kmer_in %s $\n", line);
            if(read != opt.kmer_size+1){
                ERROR("The length of kmers in %s have a different value (%zu) than the kmer size %d.", input_kmer_file, read-1, opt.kmer_size);
                exit(EXIT_FAILURE);
            }
            kmers.push_back(line);
        }
        if (line){
            free(line);
        }
        fclose(input_kmer_file_ptr);

    }else{
        generate_kmers(dna_set, "", NUM_DNA_BASES, opt.kmer_size, kmers);
    }

    uint32_t num_kmers = kmers.size();
    fprintf(stderr,"num_kmers: %d\n", num_kmers);
    if(opt.file_limit < num_kmers){
        num_kmers = opt.file_limit;
        fprintf(stderr,"only dumping %d kmers in kmer interval [%d-%d]\n", num_kmers, opt.index_start, opt.index_end);
    } else if(opt.file_limit > num_kmers - opt.index_start + 1){
        if(opt.index_end > num_kmers){
            opt.file_limit = num_kmers - opt.index_start + 1;
            opt.index_end = opt.index_start + opt.file_limit - 1 ;
        }else{
            opt.file_limit = opt.index_end - opt.index_start + 1;
        }
    }

    fprintf(stderr,"slow5_file_path: %s\n", slow5file);
    fprintf(stderr,"guppy_sam_output_file: %s\n", move_table);
    fprintf(stderr,"kmer_output_dir: %s\n", output_dir);
    fprintf(stderr,"kmer_size: %d\n", opt.kmer_size);
    fprintf(stderr,"sig_move_offset: %d\n", opt.sig_move_offset);
    fprintf(stderr,"signal_print_margin: %d\n", opt.signal_print_margin);
    fprintf(stderr,"kmer index closed interval : [%d-%d]\n", opt.index_start, opt.index_end);
    fprintf(stderr,"no.of output files: %d\n", opt.file_limit);

    std::map<std::string,FILE*> kmer_file_pointer_array;
    uint32_t count_kmers = 1;
    for(uint32_t i=opt.index_start-1;i<opt.index_end;i++){
        FILE * kmer_output = NULL;
        std::string output_path = kmer_dump + "/" + kmers[i];
        kmer_output = fopen(output_path.c_str(), "w");
        if (kmer_output == NULL){
            fprintf(stderr,"Error in opening %dth kmer-file %s\n", count_kmers, output_path.c_str());
            exit(EXIT_FAILURE);
        }
//        fprintf(kmer_output,"%.8f,", 0.0);
        kmer_file_pointer_array[kmers[i]] = kmer_output;
        count_kmers++;
    }
    std::map<std::string,uint64_t> kmer_frequency_map;
    for(uint32_t i=opt.index_start-1;i<opt.index_end;i++){
        kmer_frequency_map[kmers[i]] = 0;
    }

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

    slow5_file_t *sp = slow5_open(slow5file,"r");
    if(sp==NULL){
        fprintf(stderr,"Error in opening file %s\n", slow5file);
        exit(EXIT_FAILURE);
    }

    int ret = slow5_idx_load(sp);
    if(ret<0){
        fprintf(stderr,"Error in loading index\n");
        exit(EXIT_FAILURE);
    }

    std::string move_table_str = std::string(move_table);
    std::string extension = move_table_str.substr(move_table_str.length()-4, move_table_str.length());

    if(extension == ".paf"){
        //paf parsing
        if (flag_fastq_file_avail == 0){
            fprintf(stderr,".paf input requires an additional .fastq file\n");
            exit(EXIT_FAILURE);
        }
        INFO("%s", "sig_move_offset has no effect when using paf format's ss tag");
        process_move_table_paf(move_table, kmer_file_pointer_array, &sp, &opt, kmer_frequency_map, kmers, input_fastq_file);
    } else if(extension == ".bam" || extension == ".sam") {
        //SAM parsing
        process_move_table_bam(move_table, kmer_file_pointer_array, &sp, &opt, kmer_frequency_map, kmers);
    } else{
        process_move_table_file(move_table, kmer_file_pointer_array, &sp, &opt, kmer_frequency_map, kmers);
    }

    close_kmer_files(&kmer_file_pointer_array, kmers, &opt);

    std::string kmer_freq_file_name = std::string(std::string(output_dir)+"/freq.txt");
    FILE* kmer_freq_file = fopen(kmer_freq_file_name.c_str(), "w");
    if (kmer_freq_file == NULL){
        fprintf(stderr,"Error in opening %s\n", kmer_freq_file_name.c_str());
        exit(EXIT_FAILURE);
    }
//    fprintf(stdout,"%d\t%d\n", opt.index_start-1, opt.index_end);
    for(uint32_t i=opt.index_start-1;i<opt.index_end;i++){
        fprintf(kmer_freq_file, "%s\t%" PRIu64 "\n", kmers[i].c_str(), kmer_frequency_map[kmers[i]]);
    }

    return 0;
}

void process_move_table_file(char *move_table, std::map<std::string,FILE*> &kmer_file_pointer_array, slow5_file_t **sp_ptr, opt_t *opt_ptr, std::map<std::string,uint64_t> &kmer_frequency_map, std::vector<std::string> &kmers){
    slow5_file_t *sp = *sp_ptr;
    slow5_rec_t *rec = NULL;
    int ret=0;

    opt_t opt =  *opt_ptr;
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
    size_t num_kmers_complete = 0;
    while ((read = getline(&line, &len, basecall_fp)) != -1) {
        if(num_kmers_complete == kmers.size()){
            break;
        }
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
            if(pA < opt.pa_min || pA > opt.pa_max){
                continue;
            }
            raw_signal[i] = pA;
        }
        //calculate medmad
        if(opt.signal_scale == medmad_scale){
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
        uint32_t move_count = 0;
        size_t move_idx=0;
        size_t start_move_idx = 0;
        while(move_count < opt.sig_move_offset + 1){
            if(move_seq[move_idx]=='1'){
                move_count++;
                start_move_idx = move_idx;
            }
            move_idx++;
        }
        move_idx = start_move_idx + 1;
        size_t seq_start = opt.kmer_start_offset;
        for(;move_idx<=move_len;move_idx++){
            if(move_seq[move_idx]=='1'){
                std::string kmer = fastq_seq.substr(seq_start, opt.kmer_size);
//                fprintf(stdout,"%s\n", kmer.c_str());
                uint32_t raw_start_local = start_move_idx*stride;
                uint32_t raw_end_local = move_idx*stride;
                start_move_idx = move_idx;
                seq_start++;

                if (raw_end_local - raw_start_local > opt.max_dur){
                    continue;
                }
                if (kmer_frequency_map.find(kmer) == kmer_frequency_map.end()) {
                    continue;
                }
                if(kmer_frequency_map[kmer] == opt.sample_limit){
                    continue;
                }
                if(raw_start_local-opt.signal_print_margin < 0){
                    raw_start_local = 0;
                } else{
                    raw_start_local = raw_start_local - opt.signal_print_margin;
                }
                if(raw_end_local+opt.signal_print_margin > len_raw_signal){
                    raw_end_local = len_raw_signal;
                } else{
                    raw_end_local = raw_end_local + opt.signal_print_margin;
                }
//                fprintf(stdout,"%d\t%d\n",raw_start_local, raw_end_local);
                std::vector<double> raw_signal_local(raw_signal.begin()+raw_start_local,raw_signal.begin()+raw_end_local);

//                double max = *std::max_element(raw_signal_local.begin(),raw_signal_local.end());
//                double min = *std::min_element(raw_signal_local.begin(),raw_signal_local.end());
//                if(max-min > MAX_MIN_THRESHOLD){
//                    continue;
//                }

                size_t raw_signal_kmer_length = raw_signal_local.size();
                size_t k;
                for(k=0;k<raw_signal_kmer_length-1;k++){
                    fprintf(kmer_file_pointer_array[kmer],"%.8f,", raw_signal_local[k]);
//                    fprintf(stdout, "%.8f,", raw_signal_local[k]);
                }
                fprintf(kmer_file_pointer_array[kmer],"%.8f;", raw_signal_local[k]);
//                fprintf(stdout, "%.8f\n", raw_signal_local[k]);
                kmer_frequency_map[kmer] = kmer_frequency_map[kmer] + 1;
                if(kmer_frequency_map[kmer] == opt.sample_limit){
                    fclose(kmer_file_pointer_array[kmer]);
                    kmer_file_pointer_array[kmer] = NULL;
                    num_kmers_complete++;
                }

                if(seq_start+opt.kmer_size > fastq_seq.size()){
                    break;
                }
            }
        }
        if(opt.delimit_files == 1){
            delimit_kmer_files(&kmer_file_pointer_array, kmers, &opt);
        }
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
    return;
}
void process_move_table_paf(char *move_table, std::map<std::string,FILE*> &kmer_file_pointer_array, slow5_file_t **sp_ptr, opt_t *opt_ptr, std::map<std::string,uint64_t> &kmer_frequency_map, std::vector<std::string> &kmers, char *input_fastq_file){
    slow5_file_t *sp = *sp_ptr;
    slow5_rec_t *rec = NULL;
    int ret=0;

    opt_t opt =  *opt_ptr;
    FILE * basecall_fp = NULL;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    basecall_fp = fopen(move_table, "r");
    if (basecall_fp == NULL){
        fprintf(stderr,"Error in opening file %s\n", move_table);
        exit(EXIT_FAILURE);
    }

    faidx_t * m_fai = fai_load(input_fastq_file);
    if (m_fai == NULL){
        fprintf(stderr,"Error in loading fastq index for %s\n", input_fastq_file);
        exit(EXIT_FAILURE);
    }

    int count_reads = 0;
    size_t num_kmers_complete = 0;
    while ((read = getline(&line, &len, basecall_fp)) != -1) {
        if(num_kmers_complete == kmers.size()){
            break;
        }
        paf_rec_t *paf = parse_paf_rec(line);
        std::string read_id(paf->rid);

        int fastq_len;
        char* seq;
        // this call is not threadsafe
        seq = fai_fetch(m_fai, read_id.c_str(), &fastq_len);
        if(seq == NULL) {
            fprintf(stderr,"Error in fetching the fastq sequence for read: %s\n", paf->rid);
            exit(EXIT_FAILURE);
        }

        std::string fastq_seq(seq);
        free(seq);
        int32_t trim_offset = paf->query_start;

        VERBOSE("%s\n",read_id.c_str());
        VERBOSE("%s\t%d\t%" PRIu32 "\n", read_id.c_str(), fastq_len, trim_offset);

        ret = slow5_get(read_id.c_str(), &rec, sp);
        if(ret < 0){
            fprintf(stderr,"Error in when fetching the read\n");
            exit(EXIT_FAILURE);
        }
        uint64_t len_raw_signal = rec->len_raw_signal;
        std::vector<double> raw_signal(len_raw_signal);
        assert((uint64_t)trim_offset < len_raw_signal);
        len_raw_signal = len_raw_signal - trim_offset;
        for(uint64_t i=0;i<len_raw_signal;i++){
            double pA = TO_PICOAMPS(rec->raw_signal[i+trim_offset],rec->digitisation,rec->offset,rec->range);
            if(pA < opt.pa_min || pA > opt.pa_max){
                continue;
            }
            raw_signal[i] = pA;
        }
        //calculate medmad
        if(opt.signal_scale == medmad_scale){
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

        char *ss=paf->ss;
        size_t start_raw=paf->query_start; size_t end_raw=paf->query_end; //int len_raw_signal=paf->qlen;
        size_t start_kmer=paf->target_start; size_t end_kmer=paf->target_end;
//        int len_kmer=paf->tlen;
        // Raw signal start index for the corresponding k-mer and Raw signal end index for the corresponding k-mer
        size_t cap = MAX_LEN_KMER;
        int *st_raw_idx = (int *)malloc(sizeof(int)*cap);
        MALLOC_CHK(st_raw_idx);
        int *end_raw_idx = (int *)malloc(sizeof(int)*cap);
        MALLOC_CHK(end_raw_idx);

        //intialise to -1
        for(size_t i=0; i<cap; i++){ st_raw_idx[i]=end_raw_idx[i]=-1; }

        size_t st_k = start_kmer; size_t end_k = end_kmer; //if DNA, start k-kmer index is start_kmer column in paf and end k-kmer index is end_kmer column in paf
        int8_t rna = start_kmer > end_kmer ? 1 : 0; //if RNA start_kmer>end_kmer in paf
        if (rna){
            fprintf(stderr,"Error: RNA support is not implemented yet\n");
            exit(EXIT_FAILURE);
        }
        if(rna){ st_k = end_kmer; end_k = start_kmer; } //if RNA, start k-kmer index is end_kmer column in paf and end k-kmer index is start_kmer column in paf

        size_t i_k = st_k; size_t i_raw = start_raw; //current k-mer index and current raw signal index

        //buffer for storing digits preceding each operation and its index
        char buff[11]; int i_buff=0;
        while(*ss){
            if(*ss==',' || *ss=='I' || *ss=='D'){
                if(i_buff <= 0){ fprintf(stderr,"Bad ss: Preceding digit missing\n"); exit(1); }//if nothing in buff

                buff[i_buff]=0; //null terminate buff
                int num = atoi(buff);
                if(num < 0){ fprintf(stderr,"Bad ss: Cannot have negative numbers\n"); exit(1); }
                i_buff=0; buff[0]=0; //reset buff

                if(*ss=='I'){ //if an insertion, current raw signal index is incremented by num
                    i_raw += num;
                } else if(*ss=='D'){ //if an deletion, current k-mer index is incremented by num
                    i_k += num;
                } else if (*ss==','){ //if a mapping, increment accordingly and set raw signal indices for the current k-mer
                    end_raw_idx[i_k] = i_raw; i_raw += num;
                    st_raw_idx[i_k] = i_raw; i_k++;
                }
                if(i_k >= cap){
                    st_raw_idx=(int *)realloc(st_raw_idx, sizeof(int)*cap*2);
                    MALLOC_CHK(st_raw_idx);
                    end_raw_idx=(int *)realloc(end_raw_idx, sizeof(int)*cap*2);
                    MALLOC_CHK(end_raw_idx);
                    for(size_t i=cap; i<cap*2; i++){ st_raw_idx[i]=end_raw_idx[i]=-1; }
                    cap *= 2;
                }
            } else {
                if(!isdigit(*ss)){ fprintf(stderr,"Bad ss: A non-digit found when expected a digit\n"); exit(1); }
                buff[i_buff++]=*ss;
            }
            ss++;
        }

        if(i_raw!=end_raw){ fprintf(stderr,"Bad ss: Signal end mismatch\n"); exit(1); } //current raw signal index should be equal to end_raw
        if(i_k!=end_k){ fprintf(stderr,"Bad ss: Kmer end mismatch\n"); exit(1); } //current k-mer index should be equal to end_k

        for(size_t i=st_k; i<end_k; i++){
            if(end_raw_idx[i]==-1){
                if(st_raw_idx[i] != -1) { fprintf(stderr,"Bad ss: This shoud not have happened\n"); exit(1); }//if st_raw_idx[i] is -1, then end_raw_idx[i] should also be -1
//                printf("%s\t%d\t.\t.\n", paf->rid, rna ? len_kmer-i-1 : i);
            }else {
//                printf("%s\t%d\t%d\t%d\n", paf->rid, rna ? len_kmer-i-1 : i, end_raw_idx[i], st_raw_idx[i]);
                std::string kmer = fastq_seq.substr(i, opt.kmer_size);
//                fprintf(stdout,"%s\n", kmer.c_str());
                uint32_t raw_start_local = end_raw_idx[i];
                uint32_t raw_end_local = st_raw_idx[i];
                if (raw_end_local - raw_start_local > opt.max_dur){
                    continue;
                }
                if (kmer_frequency_map.find(kmer) == kmer_frequency_map.end()) {
                    continue;
                }
                if(kmer_frequency_map[kmer] == opt.sample_limit){
                    continue;
                }
                if(raw_start_local-opt.signal_print_margin < 0){
                    raw_start_local = 0;
                } else{
                    raw_start_local = raw_start_local - opt.signal_print_margin;
                }
                if(raw_end_local+opt.signal_print_margin > len_raw_signal){
                    raw_end_local = len_raw_signal;
                } else{
                    raw_end_local = raw_end_local + opt.signal_print_margin;
                }
//                printf("%d\t%d\n", end_raw_idx[i], st_raw_idx[i]);
                std::vector<double> raw_signal_local(raw_signal.begin()+raw_start_local,raw_signal.begin()+raw_end_local);
                size_t raw_signal_kmer_length = raw_signal_local.size();
                size_t k;
                for(k=0;k<raw_signal_kmer_length-1;k++){
                    fprintf(kmer_file_pointer_array[kmer],"%.8f,", raw_signal_local[k]);
                }
                fprintf(kmer_file_pointer_array[kmer],"%.8f;", raw_signal_local[k]);
                kmer_frequency_map[kmer] = kmer_frequency_map[kmer] + 1;
                if(kmer_frequency_map[kmer] == opt.sample_limit){
                    fclose(kmer_file_pointer_array[kmer]);
                    kmer_file_pointer_array[kmer] = NULL;
                    num_kmers_complete++;
                }
                if(i+opt.kmer_size > fastq_seq.size()){
                    break;
                }
            }
        }

        free(st_raw_idx);
        free(end_raw_idx);

        if(opt.delimit_files == 1){
            delimit_kmer_files(&kmer_file_pointer_array, kmers, &opt);
        }
        if(count_reads == PROGRESS_BATCH_SIZE){
            fprintf(stderr,"*");
            count_reads = 0;
        }
        count_reads++;
        free_paf_rec(paf);
    }
    fclose(basecall_fp);
    if (line){
        free(line);
    }
    return;
}

paf_rec_t *parse_paf_rec(char *buffer){

    char *pch=NULL;

    paf_rec_t *paf = (paf_rec_t *)malloc(sizeof(paf_rec_t));
    MALLOC_CHK(paf);

    //read name
    pch = strtok (buffer,"\t\r\n"); assert(pch!=NULL);
    paf->rid = strdup(pch);

    //readlen
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->qlen = atoi(pch);

    //query start
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->query_start = atoi(pch);

    //query end
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->query_end= atoi(pch);

    //relative strand
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    if(strcmp(pch,"+")==0){
        paf->strand=0;
    }
    else if(strcmp(pch,"-")==0){
        paf->strand=1;
    }
    else{
        assert(0);
    }

    //targetname
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->tid = strdup(pch);

    //target len
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->tlen = atoi(pch);

    //target start
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->target_start = atoi(pch);

    //target end
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->target_end= atoi(pch);

    //num residue
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);

    //num block
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);

    //mapq
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->mapq = atoi(pch);

    paf->ss = NULL;
    while((pch = strtok(NULL,"\t\r\n"))){
        if(strncmp("ss:Z:",pch,5)==0){
            //fprintf(stderr,"ss:Z:%s\n",pch);
            paf->ss = strdup(pch+5);
        }
    }
    //error check
    if(paf->ss==NULL){
        ERROR("ss:Z: tag not found in paf record for %s", paf->rid);
        exit(EXIT_FAILURE);
    }

    return paf;
}

void free_paf_rec(paf_rec_t *paf){
    free(paf->rid);
    free(paf->tid);
    free(paf->ss);
    free(paf);
}

void process_move_table_bam(char *move_table, std::map<std::string,FILE*> &kmer_file_pointer_array, slow5_file_t **sp_ptr, opt_t *opt_ptr, std::map<std::string,uint64_t> &kmer_frequency_map, std::vector<std::string> &kmers){
    slow5_file_t *sp = *sp_ptr;
    slow5_rec_t *rec = NULL;
    int ret=0;

    opt_t opt =  *opt_ptr;
    htsFile* bam_fp = sam_open(move_table, "r"); //open bam file
    F_CHK(bam_fp, move_table);

    bam_hdr_t* bam_hdr = sam_hdr_read(bam_fp);
    NULL_CHK(bam_hdr);

    bam1_t *aln = bam_init1(); //initialize an alignment
    if (!aln) {
        ERROR("%s failed\n", "bam_init1)(");
        exit(EXIT_FAILURE);
    }
    int count_reads = 0;
    size_t num_kmers_complete = 0;
    while((ret = sam_read1(bam_fp, bam_hdr, aln)) >= 0){
        if(num_kmers_complete == kmers.size()){
            break;
        }
        const char tag_ns[] = {'n', 's'};
        const uint8_t * ns_ptr = bam_aux_get(aln, tag_ns);
        if(!ns_ptr){
            ERROR("tag 'ns' is not found. Please check your SAM/BAM file: %s", "");
            exit(EXIT_FAILURE);
        }
        uint64_t signal_len = bam_aux2i(ns_ptr);

        const char tag_ts[] = {'t', 's'};
        const uint8_t * ts_ptr = bam_aux_get(aln, tag_ts);
        if(!ts_ptr){
            ERROR("tag 'ts' is not found. Please check your SAM/BAM file: %s", "");
            exit(EXIT_FAILURE);
        }
        uint64_t trim_offset = bam_aux2i(ts_ptr);

        const char tag_mv[] = {'m', 'v'};
        const uint8_t * mv_array = bam_aux_get(aln, tag_mv);
        if(!mv_array){
            ERROR("NULL returned for tag mv: %s", "");
            exit(EXIT_FAILURE);
        }
        int stride;
        uint32_t move_len;
        if( (char)mv_array[0] == 'B' && (char)mv_array[1] == 'c') {
            move_len = bam_auxB_len(mv_array);
            if (move_len == 0) {
                ERROR("mv array length is 0: %s", "");
                exit(EXIT_FAILURE);
            }
            LOG_DEBUG("len_mv:%d\n", move_len);
            stride = bam_auxB2i(mv_array, 0);
            if ( stride != EXPECTED_STRIDE) {
                ERROR("expected stride of %d is missing.", EXPECTED_STRIDE);
                exit(EXIT_FAILURE);
            }
        } else{
            ERROR("tag 'mv' specification is incorrect%s", "");
            exit(EXIT_FAILURE);
        }

        std::string read_id(bam_get_qname(aln));
        int32_t fastq_len = aln->core.l_qseq;
        uint8_t*  bam_seq_ptr = bam_get_seq(aln);
        std::string alphabet = "NACNGNNNT";
        std::string fastq_seq = "";
        for(int i=0; i<fastq_len; i++){
//            fprintf(stderr,"%d", bam_seqi(bam_seq_ptr,i));
//            fprintf(stderr,"%c", alphabet[bam_seqi(bam_seq_ptr,i)]);
            fastq_seq.push_back(alphabet[bam_seqi(bam_seq_ptr,i)]);
        }

        VERBOSE("%s\n",read_id.c_str());
        VERBOSE("%s\t%" PRIu32 "\t%" PRIu64 "\t%" PRIu64 "\t%d\n",read_id.c_str(),fastq_len,signal_len,trim_offset,stride);

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
        int flag_skip_signal = 0;
        for(uint64_t i=0;i<len_raw_signal;i++){
            double pA = TO_PICOAMPS(rec->raw_signal[i+trim_offset],rec->digitisation,rec->offset,rec->range);
            if(pA < opt.pa_min || pA > opt.pa_max){
                flag_skip_signal = 1;
                break;
            }
            raw_signal[i] = pA;
        }
        if(flag_skip_signal == 1){
            continue;
        }
        //calculate medmad
        if(opt.signal_scale == medmad_scale){
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
        uint32_t move_count = 0;
        size_t move_idx=0;
        size_t start_move_idx = 0;
        while(move_count < opt.sig_move_offset + 1){
            int8_t value = bam_auxB2i(mv_array, move_idx+1);
            if(value == 1){
                move_count++;
                start_move_idx = move_idx;
            }
            move_idx++;
        }
        move_idx = start_move_idx + 1;
        size_t seq_start = opt.kmer_start_offset;
        for(;move_idx<=move_len;move_idx++){
            int8_t value = bam_auxB2i(mv_array, move_idx+1);
            if(value == 1){
                std::string kmer = fastq_seq.substr(seq_start, opt.kmer_size);
                uint32_t raw_start_local = start_move_idx*stride;
                uint32_t raw_end_local = move_idx*stride;
                start_move_idx = move_idx;
                seq_start++;

                if (raw_end_local - raw_start_local > opt.max_dur){
                    continue;
                }
                if (kmer_frequency_map.find(kmer) == kmer_frequency_map.end()) {
                    continue;
                }
                if(kmer_frequency_map[kmer] == opt.sample_limit){
                    continue;
                }
                if(raw_start_local-opt.signal_print_margin < 0){
                    raw_start_local = 0;
                } else{
                    raw_start_local = raw_start_local - opt.signal_print_margin;
                }
                if(raw_end_local+opt.signal_print_margin > len_raw_signal){
                    raw_end_local = len_raw_signal;
                } else{
                    raw_end_local = raw_end_local + opt.signal_print_margin;
                }
//                fprintf(stdout,"%d\t%d\n", raw_start_local, raw_end_local);
                std::vector<double> raw_signal_local(raw_signal.begin()+raw_start_local,raw_signal.begin()+raw_end_local);

//                double max = *std::max_element(raw_signal_local.begin(),raw_signal_local.end());
//                double min = *std::min_element(raw_signal_local.begin(),raw_signal_local.end());
//                if(max-min > MAX_MIN_THRESHOLD){
//                    continue;
//                }
                size_t raw_signal_kmer_length = raw_signal_local.size();
                size_t k;
                for(k=0;k<raw_signal_kmer_length-1;k++){
                    fprintf(kmer_file_pointer_array[kmer],"%.8f,", raw_signal_local[k]);
//                    fprintf(stdout, "%.8f,", raw_signal_local[k]);
                }
                fprintf(kmer_file_pointer_array[kmer],"%.8f;", raw_signal_local[k]);
//                fprintf(stdout, "%.8f\n", raw_signal_local[k]);
                kmer_frequency_map[kmer] = kmer_frequency_map[kmer] + 1;
                if(kmer_frequency_map[kmer] == opt.sample_limit){
                    fclose(kmer_file_pointer_array[kmer]);
                    kmer_file_pointer_array[kmer] = NULL;
                    num_kmers_complete++;
                }
                if(seq_start+opt.kmer_size > fastq_seq.size()){
                    break;
                }
            }
        }
        if(opt.delimit_files == 1){
            delimit_kmer_files(&kmer_file_pointer_array, kmers, &opt);
        }
        if(count_reads == PROGRESS_BATCH_SIZE){
            fprintf(stderr,"*");
            count_reads = 0;
        }
        count_reads++;
    }
    bam_destroy1(aln);
    bam_hdr_destroy(bam_hdr);
    sam_close(bam_fp);
    return;
}






























