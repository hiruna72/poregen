/* @file  kmer_freq.c
**
** @@ Hiruna Samarakoon - hiruna72@gmail.com
** ./poregen kmer_freq 9 pass.fastq | grep  -e $'.*\t0' | wc
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
#include <map>
#include <algorithm>

static struct option long_options[] = {
    {"sort", required_argument, 0, 0},               //0
    {"print_absent_kmers", required_argument, 0, 0}, //1
    {"help", no_argument, 0, 'h'},                   //2
    {"version", no_argument, 0, 'V'},                //3
    {"output",required_argument, 0, 'o'},            //4 output to a file [stdout]
    {"debug-break",required_argument, 0, 0},         //5 break after processing the first batch (used for debugging)
    {0, 0, 0, 0}};


static inline void print_help_msg(FILE *fp_help, opt_t opt){
    fprintf(fp_help,"Usage: poregen kmer_freq kmer_size reads.fastq\n");
    fprintf(fp_help,"\nbasic options:\n");
    fprintf(fp_help,"   --sort INT                 sort based on frequency (0-no sorting, 1-ascend, 2-descend) [0] \n");
    fprintf(fp_help,"   --print_absent_kmers INT   print kmers with 0 frequency (0-do not print, 1-print) [1] \n");
    fprintf(fp_help,"   -o FILE                    output to file [stdout]\n");
    fprintf(fp_help,"   --verbose INT              verbosity level [%d]\n",(int)get_log_level());
    fprintf(fp_help,"   --version                  print version\n");

    fprintf(fp_help,"\nadvanced options:\n");
    fprintf(fp_help,"   --debug-break INT          break after processing the specified no. of batches\n");
}

//parse yes or no arguments : taken from minimap2
static inline void yes_or_no(opt_t* opt, uint64_t flag, int long_idx,
                             const char* arg,
                             int yes_to_set)
{
    if (yes_to_set) {
        if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) {
            opt->flag |= flag;
        } else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) {
            opt->flag &= ~flag;
        } else {
            WARNING("option '--%s' only accepts 'yes' or 'no'.",
                    long_options[long_idx].name);
        }
    } else {
        if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) {
            opt->flag &= ~flag;
        } else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) {
            opt->flag |= flag;
        } else {
            WARNING("option '--%s' only accepts 'yes' or 'no'.",
                    long_options[long_idx].name);
        }
    }
}

int kmer_freq(int argc, char* argv[]) {

    //signal(SIGSEGV, sig_handler);

    const char* optstring = "v:o:hV";

    int longindex = 0;
    int32_t c = -1;

    char *fastq_file = NULL;
    int flag_sort_freq = 0;
    int flag_print_absent_kmers = 1;
    char dna_set[] = {'A', 'C', 'G', 'T'};

    FILE *fp_help = stderr;

    opt_t opt;
    init_opt(&opt); //initialise options to defaults

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {
        if (c=='v'){
            int v = atoi(optarg);
            set_log_level((enum poregen_log_level_opt)v);
        } else if (c=='V'){
            fprintf(stdout,"subtool0 %s\n",POREGEN_VERSION);
            exit(EXIT_SUCCESS);
        } else if (c=='h'){
            fp_help = stdout;
            fp_help = stdout;
        } else if(c == 0 && longindex == 0){ //debug break
            if(atoi(optarg) != 0 && atoi(optarg) != 1 && atoi(optarg) != 2){
                ERROR("sort argument must be 0,1 or 2 You entered %d", atoi(optarg));
                exit(EXIT_FAILURE);
            }
            flag_sort_freq = atoi(optarg);
        } else if(c == 0 && longindex == 1){ //debug break
            if(atoi(optarg) != 0 && atoi(optarg) != 1){
                ERROR("print_absent_kmers flag must be 0 or 1 You entered %d", atoi(optarg));
                exit(EXIT_FAILURE);
            }
            flag_print_absent_kmers = atoi(optarg);
        } else if(c == 0 && longindex == 4){ //debug break
            opt.debug_break = atoi(optarg);
        }
    }

    // No arguments given
    if (argc - optind != 2 || fp_help == stdout) {
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }
    opt.kmer_size = atoi(argv[optind++]);
    fastq_file = argv[optind];
    fprintf(stderr,"kmer_size: %d\n", opt.kmer_size);

    if (fastq_file == NULL) {
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }


    std::vector<std::string> kmers;
    generate_kmers(dna_set, "", NUM_DNA_BASES, opt.kmer_size, kmers);
    std::map<std::string, uint64_t> kmer_frequency_map;
    uint32_t num_kmers = kmers.size();
    fprintf(stderr,"num_kmers: %d\n", num_kmers);

    for(uint32_t i=0;i<num_kmers;i++){
        kmer_frequency_map[kmers[i]] = 0;
    }

    FILE* input_fastq_file_ptr = fopen(fastq_file, "r");
    if (input_fastq_file_ptr == NULL){
        fprintf(stderr,"Error in opening file %s\n", fastq_file);
        exit(EXIT_FAILURE);
    }
    char * line = NULL;
    size_t len = 0;
    ssize_t read_len;

    size_t count_line = 0;
    while ((read_len = getline(&line, &len, input_fastq_file_ptr)) != -1) {
        line[read_len-1] = '\0';
        if(count_line%4 == 1){
//            fprintf(stderr,"kmer_in %s$\n", line);
            std::string fastq_seq = std::string(line);
            for(uint32_t i=0; i<read_len-opt.kmer_size; i++){
                std::string kmer = fastq_seq.substr(i, opt.kmer_size);
                kmer_frequency_map[kmer] = kmer_frequency_map[kmer] + 1;
            }
        }
        count_line++;
    }
    if (line){
        free(line);
    }

    fclose(input_fastq_file_ptr);

    // create an empty vector of pairs
    std::vector<std::pair<std::string, uint64_t>> vec;

    // copy key-value pairs from the map to the vector
    std::copy(kmer_frequency_map.begin(),
              kmer_frequency_map.end(),
              std::back_inserter<std::vector<std::pair<std::string, uint64_t>>>(vec));

    if(flag_sort_freq == 1){
        // sort the vector by increasing the order of its pair's second value
        // if the second value is equal, order by the pair's first value
        std::sort(vec.begin(), vec.end(), [](const std::pair<std::string, uint64_t> &l, const std::pair<std::string, uint64_t> &r){
            if (l.second != r.second) {
                return l.second < r.second;
            }
            return l.first < r.first;
        });
    } else if(flag_sort_freq == 2){
        std::sort(vec.begin(), vec.end(), [](const std::pair<std::string, uint64_t> &l, const std::pair<std::string, uint64_t> &r){
            if (l.second != r.second) {
                return l.second > r.second;
            }
            return l.first > r.first;
        });
    }

    // print the vector
    for (auto const &pair: vec) {
        if(flag_print_absent_kmers == 0 && pair.second == 0){
            continue;
        }
        fprintf(stdout, "%s\t%" PRIu64 "\n", pair.first.c_str(), pair.second);
    }

    return 0;
}
