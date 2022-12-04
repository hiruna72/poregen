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

#ifndef DISABLE_KSORT
#include "ksort.h"
KSORT_INIT_GENERIC(float)
KSORT_INIT_GENERIC(double)
KSORT_INIT_GENERIC(int16_t)
KSORT_INIT_GENERIC(int)
#endif

static struct option long_options[] = {
    {"kmer_size", required_argument, 0, 'k'},      //0 kmer_size [6]
    {"scaling", required_argument, NULL, 0},       //1 scaling 1-medmad
    {"margin", required_argument, NULL, 0},        //2
    {"sample_limit", required_argument, NULL, 0},  //3
    {"file_limit", required_argument, NULL, 0},    //4
    {"kmer_file", required_argument, NULL, 0},     //5
    {"index_start", required_argument, NULL, 0},   //6
    {"index_end", required_argument, NULL, 0},     //7
    {"verbose", required_argument, 0, 'v'},        //8 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //9
    {"version", no_argument, 0, 'V'},              //10
    {"debug-break",required_argument, 0, 0},       //11 break after processing the first batch (used for debugging)
    {0, 0, 0, 0}};


static inline void print_help_msg(FILE *fp_help, opt_t opt){
    fprintf(fp_help,"Usage: poregen gmove reads.blow5 move_table output_dir\n");
    fprintf(fp_help,"\nbasic options:\n");
    fprintf(fp_help,"   -k INT                     kmer_size [%d]\n",opt.kmer_size);
    fprintf(fp_help,"   --scaling INT              scaling [%d] (0-no scaling, 1-medmad)\n",opt.signal_scale);
    fprintf(fp_help,"   --margin INT               signal print margin on both sides of the sub signal[%u] \n",opt.signal_print_margin);
    fprintf(fp_help,"   --sample_limit INT         maximum number of instances to output for a kmer [%u] \n",opt.sample_limit);
    fprintf(fp_help,"   --file_limit INT           maximum number of kmer files to output [%u] \n",opt.file_limit);
    fprintf(fp_help,"   --kmer_file FILE           kmer file (optional) \n");
    fprintf(fp_help,"   --index_start INT          1-based closed interval index of start kmer [%u] \n",opt.file_limit);
    fprintf(fp_help,"   --index_end INT            1-based closed interval index of end kmer [%u] \n",opt.file_limit);
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
        fclose(kmer_file_pointer_array[kmers[i]]);
    }
}


void delimit_kmer_files(std::map<std::string, FILE *> *kmer_file_pointer_array_ptr, std::vector<std::string>& kmers, opt_t* opt) {
    std::map<std::string, FILE *> kmer_file_pointer_array = *kmer_file_pointer_array_ptr;
    for(uint32_t i=opt->index_start-1;i<opt->index_end;i++){
        fprintf(kmer_file_pointer_array[kmers[i]],":");
    }
}

int gmove(int argc, char* argv[]) {

//    fprintf(stderr,"%d", get_log_level());

    const char* optstring = "k:";

    int longindex = 0;
    int32_t c = -1;
    int flag_kmer_file_avail = 0;

    char *slow5file = NULL;
    char *move_table = NULL;
    char *output_dir = NULL;
    char *input_kmer_file = NULL;
    int signal_scale = 0;

    char dna_set[] = {'A', 'C', 'G', 'T'};

    int kmer_start_offset = 0;
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
        } else if (c=='v'){
            int v = atoi(optarg);
            set_log_level((enum poregen_log_level_opt)v);
        } else if (c=='V'){
            fprintf(stdout,"gmove %s\n",POREGEN_VERSION);
            exit(EXIT_SUCCESS);
        } else if (c=='h'){
            fp_help = stdout;
            fp_help = stdout;
        } else if(c == 0 && longindex == 1){ //scaling
            signal_scale = atoi(optarg);
        } else if (c == 0 && longindex == 2){ //margin
            if(atoi(optarg) < 0){
                ERROR("Signal print margin should be non negative. You entered %d", atoi(optarg));
                exit(EXIT_FAILURE);
            }
            opt.signal_print_margin = atoi(optarg);
        } else if (c == 0 && longindex == 3){ //kmer_limit
            if(atoi(optarg) < 0){
                ERROR("Maximum number of instances to output for a kmer should be non negative. You entered %d", atoi(optarg));
                exit(EXIT_FAILURE);
            }
            opt.sample_limit = atoi(optarg);
        } else if (c == 0 && longindex == 4){ //dump_limit
            if(atoi(optarg) < 0){
                ERROR("Maximum number of kmer files to output should be non negative. You entered %d", atoi(optarg));
                exit(EXIT_FAILURE);
            }
            opt.file_limit = atoi(optarg);
            opt.index_end = opt.index_start + opt.file_limit;
        } else if (c == 0 && longindex == 5){ //dump_limit
            input_kmer_file = optarg;
            flag_kmer_file_avail = 1;
        } else if (c == 0 && longindex == 6){ //dump_limit
            if(atoi(optarg) < 1){
                ERROR("kmer index start should be a positive number. You entered %d", atoi(optarg));
                exit(EXIT_FAILURE);
            }
            opt.index_start = atoi(optarg);
            opt.file_limit = opt.index_end - opt.index_start + 1;
        } else if (c == 0 && longindex == 7){ //dump_limit
            if(atoi(optarg) < 1){
                ERROR("kmer index end should be a positive number. You entered %d", atoi(optarg));
                exit(EXIT_FAILURE);
            }
            opt.index_end = atoi(optarg);
            opt.file_limit = opt.index_end - opt.index_start + 1;
        } else if (c == 0 && longindex == 11){ //debug break
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

    int move_start_offset = opt.kmer_size - 1;
    fprintf(stderr,"slow5_file_path: %s\n", slow5file);
    fprintf(stderr,"guppy_sam_output_file: %s\n", move_table);
    fprintf(stderr,"kmer_output_dir: %s\n", output_dir);
    fprintf(stderr,"kmer_size: %d\n", opt.kmer_size);
    fprintf(stderr,"move_start_offset: %d\n", move_start_offset);

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
        fprintf(kmer_output,"%.8f,", 0.0);
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
//                fprintf(stdout,"%s\n", kmer.c_str());
                int raw_start_local = start_move_idx*stride;
                int raw_end_local = move_idx*stride;
                start_move_idx = move_idx;
                seq_start++;

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
                if(seq_start+opt.kmer_size > fastq_seq.size()){
                    break;
                }
            }
        }
        delimit_kmer_files(&kmer_file_pointer_array, kmers, &opt);
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

    close_kmer_files(&kmer_file_pointer_array, kmers, &opt);

    std::string kmer_freq_file_name = std::string(std::string(output_dir)+"/freq.txt");
    FILE* kmer_freq_file = fopen(kmer_freq_file_name.c_str(), "w");
    if (kmer_freq_file == NULL){
        fprintf(stderr,"Error in opening %s\n", kmer_freq_file_name.c_str());
        exit(EXIT_FAILURE);
    }
    for(uint32_t i=opt.index_start-1;i<opt.index_end;i++){
        fprintf(kmer_freq_file, "%s\t%" PRIu64 "\n", kmers[i].c_str(), kmer_frequency_map[kmers[i]]);
    }

    return 0;
}
