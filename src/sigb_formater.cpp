/* @file  sigb_formater.c
** properly format move table
** @@ Hiruna Samarakoon
 * https://gist.github.com/PoisonAlien/350677acc03b2fbf98aa
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
#include <htslib/sam.h>

#define EXPECTED_STRIDE 5

static struct option long_options[] = {
    {"kmer", required_argument, 0, 'k'},           //0 kmer size [9]
    {"", no_argument, 0, 'c'},                     //1 optional paf format
    {"threads", required_argument, 0, 't'},        //2 number of threads [8]
    {"batchsize", required_argument, 0, 'K'},      //3 batchsize - number of reads loaded at once [512]
    {"max-bytes", required_argument, 0, 'B'},      //4 batchsize - number of bytes loaded at once
    {"verbose", required_argument, 0, 'v'},        //5 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //6
    {"version", no_argument, 0, 'V'},              //7
    {"output",required_argument, 0, 'o'},          //8 output to a file [stdout]
    {"debug-break",required_argument, 0, 0},       //9 break after processing the first batch (used for debugging)
    {0, 0, 0, 0}};


static inline void print_help_msg(FILE *fp_help, opt_t opt){
    fprintf(fp_help,"Usage: poregen sigbformater basecalled.SAM/BAM\n");
    fprintf(fp_help,"\nbasic options:\n");
    fprintf(fp_help,"   -k, --kmer                 kmer size [%d]\n",opt.kmer_size);
    fprintf(fp_help,"   -c                         write move table in paf format\n");
    fprintf(fp_help,"   -t INT                     number of processing threads [%d]\n",opt.num_thread);
    fprintf(fp_help,"   -K INT                     batch size (max number of reads loaded at once) [%d]\n",opt.batch_size);
    fprintf(fp_help,"   -B FLOAT[K/M/G]            max number of bytes loaded at once [%.1fM]\n",opt.batch_size_bytes/(float)(1000*1000));
    fprintf(fp_help,"   -h                         help\n");
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

int sigb_formater(int argc, char* argv[]) {

    //signal(SIGSEGV, sig_handler);

    const char* optstring = "ct:B:K:v:o:hV";

    int longindex = 0;
    int32_t c = -1;

    char *bam_file_name = NULL;

    FILE *fp_help = stderr;

    opt_t opt;
    init_opt(&opt); //initialise options to defaults

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {

        if (c == 'k') {
            opt.kmer_size = atoi(optarg);
        } else if (c == 'c') {
            opt.use_paf_format = 1;
        } else if (c == 'B') {
            opt.batch_size_bytes = mm_parse_num(optarg);
            if(opt.batch_size_bytes<=0){
                ERROR("%s","Maximum number of bytes should be larger than 0.");
                exit(EXIT_FAILURE);
            }
        } else if (c == 'K') {
            opt.batch_size = atoi(optarg);
            if (opt.batch_size < 1) {
                ERROR("Batch size should larger than 0. You entered %d",opt.batch_size);
                exit(EXIT_FAILURE);
            }
        } else if (c == 't') {
            opt.num_thread = atoi(optarg);
            if (opt.num_thread < 1) {
                ERROR("Number of threads should larger than 0. You entered %d", opt.num_thread);
                exit(EXIT_FAILURE);
            }
        } else if (c=='v'){
            int v = atoi(optarg);
            set_log_level((enum poregen_log_level_opt)v);
        } else if (c=='V'){
            fprintf(stdout,"subtool0 %s\n",POREGEN_VERSION);
            exit(EXIT_SUCCESS);
        } else if (c=='h'){
            fp_help = stdout;
        } else if(c == 0 && longindex == 8){ //debug break
            opt.debug_break = atoi(optarg);
        }
    }

    // No arguments given
    if (argc - optind != 1 || fp_help == stdout) {
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }
    bam_file_name = argv[optind];

    if (bam_file_name == NULL) {
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "bam_file : %s\n", bam_file_name);

    htsFile* bam_fp = sam_open(bam_file_name, "r"); //open bam file
    bam_hdr_t* bam_hdr = sam_hdr_read(bam_fp);
    bam1_t *aln = bam_init1(); //initialize an alignment

    while(sam_read1(bam_fp, bam_hdr, aln) > 0){
        LOG_DEBUG("%s\n", aln->data);
        uint32_t len_seq = aln->core.l_qseq; //length of the read.

        const char tag_ns[] = {'n', 's'};
        const uint8_t * ns_ptr = bam_aux_get(aln, tag_ns);
        if(!ns_ptr){
            ERROR("NULL returned for tag ns: %s", strerror(errno));
            return -1;
        }
        int64_t ns;
        if( (char)ns_ptr[0] == 'i'){
            ns = bam_aux2i(ns_ptr);
        }
        LOG_DEBUG("ns: %" PRIu64 "\n", ns);

        const char tag_ts[] = {'t', 's'};
        const uint8_t * ts_ptr = bam_aux_get(aln, tag_ts);
        if(!ts_ptr){
            ERROR("NULL returned for tag ts: %s", strerror(errno));
            return -1;
        }
        uint64_t ts;
        if( (char)ns_ptr[0] == 'i'){
            ts = bam_aux2i(ts_ptr);
        }
        LOG_DEBUG("ts: %" PRIu64 "\n", ts);

        const char tag_mv[] = {'m', 'v'};
        const uint8_t * mv_array = bam_aux_get(aln, tag_mv);
        if(!mv_array){
            ERROR("NULL returned for tag mv: %s", strerror(errno));
            return -1;
        }
        if( (char)mv_array[0] == 'B' && (char)mv_array[1] == 'c'){
            uint32_t len_mv = bam_auxB_len(mv_array);
            if(len_mv == 0){
                ERROR("mv array length is 0: %s", strerror(errno));
                return -1;
            }
            LOG_DEBUG("len_mv:%d\n", len_mv);
            if(bam_auxB2i(mv_array, 0) != EXPECTED_STRIDE){
                ERROR("expected stride of %d is missing.", EXPECTED_STRIDE);
                return -1;
            }
            if(opt.use_paf_format == 0){
                uint32_t move_count = 0;
                uint32_t i = 1;
                while(move_count < opt.kmer_size){
                    int8_t value = bam_auxB2i(mv_array, i);
                    if(value == 1){
                        move_count++;
                    }
                    i++;
                }
                uint64_t end_idx = ts + (i-1) * EXPECTED_STRIDE;
                uint64_t start_idx = end_idx - EXPECTED_STRIDE;
                uint32_t kmer_idx = 0;
                for(; i<len_mv; i++){
                    int8_t value = bam_auxB2i(mv_array, i);
                    LOG_DEBUG("%d", value);

                    if(value == 1 || i == len_mv-1){
                        fprintf(stdout, "%s\t", aln->data);
                        fprintf(stdout, "%" PRIu32 "\t", kmer_idx);
                        fprintf(stdout, "%" PRId8 "\t", i-1);
                        fprintf(stdout, "%" PRIu64 "\t", start_idx);
                        if( i == len_mv-1){
                            fprintf(stdout, "%" PRIu64 "\n", end_idx+EXPECTED_STRIDE);
                        } else{
                            fprintf(stdout, "%" PRIu64 "\n", end_idx);
                        }
                        start_idx = end_idx;
                        kmer_idx++;
                        len_seq--;
                    }
                    end_idx = end_idx + EXPECTED_STRIDE;
                }

            } else{
                fprintf(stderr,"paf is not yet implemented");
            }

        } else{
            ERROR("mv tag specification is incorrect%s", "");
            return -1;
        }

    }

    bam_destroy1(aln);
    bam_hdr_destroy(bam_hdr);
    sam_close(bam_fp);
    return 0;
}
