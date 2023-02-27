/* @file  reform.c
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
    {"kmer_length", required_argument, 0, 'k'},         //0 kmer length [9]
    {"sig_move_offset", required_argument, 0, 'm'},     //1 signal move offset  [9]
    {"", no_argument, 0, 'c'},                          //2 optional paf format
    {"threads", required_argument, 0, 't'},             //3 number of threads [8]
    {"batchsize", required_argument, 0, 'K'},           //4 batchsize - number of reads loaded at once [512]
    {"max-bytes", required_argument, 0, 'B'},           //5 batchsize - number of bytes loaded at once
    {"verbose", required_argument, 0, 'v'},             //6 verbosity level [1]
    {"help", no_argument, 0, 'h'},                      //7
    {"version", no_argument, 0, 'V'},                   //8
    {"output",required_argument, 0, 'o'},               //9 output to a file [stdout]
    {"debug-break",required_argument, 0, 0},            //10 break after processing the first batch (used for debugging)
    {0, 0, 0, 0}};


static inline void print_help_msg(FILE *fp_help, opt_t opt){
    fprintf(fp_help,"Usage: poregen reform basecalled.SAM/BAM\n");
    fprintf(fp_help,"\nbasic options:\n");
    fprintf(fp_help,"   -k, --kmer_length          kmer length [%d]\n",opt.kmer_size);
    fprintf(fp_help,"   -m, --sig_move_offset      signal move offset [%d]\n",opt.move_start_offset);
    fprintf(fp_help,"   -c                         write move table in paf format\n");
//    fprintf(fp_help,"   -t INT                     number of processing threads [%d]\n",opt.num_thread);
//    fprintf(fp_help,"   -K INT                     batch size (max number of reads loaded at once) [%d]\n",opt.batch_size);
//    fprintf(fp_help,"   -B FLOAT[K/M/G]            max number of bytes loaded at once [%.1fM]\n",opt.batch_size_bytes/(float)(1000*1000));
    fprintf(fp_help,"   -h                         help\n");
    fprintf(fp_help,"   -o FILE                    output to file [stdout]\n");
    fprintf(fp_help,"   --verbose INT              verbosity level [%d]\n",(int)get_log_level());
    fprintf(fp_help,"   --version                  print version\n");

//    fprintf(fp_help,"\nadvanced options:\n");
//    fprintf(fp_help,"   --debug-break INT          break after processing the specified no. of batches\n");
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

int reform(int argc, char* argv[]) {

    //signal(SIGSEGV, sig_handler);

    const char* optstring = "k:m:ct:B:K:v:o:hV";

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
        } else if (c == 'm') {
            opt.move_start_offset = atoi(optarg);
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
        } else if (c=='o'){
            opt.arg_fname_out = optarg;
        } else if (c=='V'){
            fprintf(stdout,"subtool0 %s\n",POREGEN_VERSION);
            exit(EXIT_SUCCESS);
        } else if (c=='h'){
            fp_help = stdout;
        } else if(c == 0 && longindex == 8){ //debug break
            opt.debug_break = atoi(optarg);
        }
    }


    if(opt.kmer_size < 1){
        ERROR("kmer length must be non zero%s", "")
        return -1;
    }
    if(opt.move_start_offset < 1){
        ERROR("signal move offset value must be non zero%s", "")
        return -1;
    }
    if(opt.kmer_size < opt.move_start_offset){
        ERROR("signal move offset value must not be larger than the kmer length%s", "")
        return -1;
    }

    // No arguments given
    if (argc - optind != 1 || fp_help == stdout) {
        ERROR("%s", "not enough arguments")
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        return EXIT_FAILURE;
    }
    bam_file_name = argv[optind];

    if (bam_file_name == NULL) {
        ERROR("%s", "input file name cannot be null")
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        return EXIT_FAILURE;

    }

    fprintf(stderr, "bam_file : %s\n", bam_file_name);
    fprintf(stderr, "kmer length : %" PRIu32 "\n", opt.kmer_size);
    fprintf(stderr, "signal move offset : %" PRIu32 "\n", opt.move_start_offset);
    if(opt.use_paf_format){
        fprintf(stderr, "output format : %s\n", "paf");
    } else{
        fprintf(stderr, "output format : %s\n", "tsv");
    }

    // Parse output argument
    if (opt.arg_fname_out != NULL) {
        LOG_DEBUG("opening output file%s","");
        // Create new file or
        // Truncate existing file
        FILE *new_file;
        new_file = fopen(opt.arg_fname_out, "w");
        F_CHK(new_file, opt.arg_fname_out);
        opt.f_out = new_file;
    }

    htsFile* bam_fp = sam_open(bam_file_name, "r"); //open bam file
    F_CHK(bam_fp, bam_file_name);

    bam_hdr_t* bam_hdr = sam_hdr_read(bam_fp);
    NULL_CHK(bam_hdr);

    bam1_t *aln = bam_init1(); //initialize an alignment
    if (!aln) {
        ERROR("%s failed\n", "bam_init1)(");
        return EXIT_FAILURE;
    }

    while(sam_read1(bam_fp, bam_hdr, aln) > 0){
        LOG_DEBUG("%s\n", aln->data);
        uint32_t len_seq = aln->core.l_qseq; //length of the read.
        len_seq = len_seq - opt.kmer_size + 1; //to get the number of kmers

        const char tag_ns[] = {'n', 's'};
        const uint8_t * ns_ptr = bam_aux_get(aln, tag_ns);
        if(!ns_ptr){
            ERROR("NULL returned for tag ns: %s", "");
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
            ERROR("NULL returned for tag ts: %s", "");
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
            ERROR("NULL returned for tag mv: %s", "");
            return -1;
        }
        if( (char)mv_array[0] == 'B' && (char)mv_array[1] == 'c'){
            uint32_t len_mv = bam_auxB_len(mv_array);
            if(len_mv == 0){
                ERROR("mv array length is 0: %s", "");
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
                while(move_count < opt.move_start_offset){
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

                    if(len_seq > 0 && (value == 1 || i == len_mv-1)){
                        fprintf(opt.f_out, "%s\t", aln->data);
                        fprintf(opt.f_out, "%" PRIu32 "\t", kmer_idx);
                        fprintf(opt.f_out, "%" PRId8 "\t", i-1);
                        fprintf(opt.f_out, "%" PRIu64 "\t", start_idx);
                        if( i == len_mv-1){
                            uint32_t l_duration =  end_idx+EXPECTED_STRIDE+(ns-i*EXPECTED_STRIDE);
                            fprintf(opt.f_out, "%" PRIu32 "", l_duration); //ss
                        } else{
                            fprintf(opt.f_out, "%" PRIu64 "\n", end_idx);
                        }
                        start_idx = end_idx;
                        kmer_idx++;
                        len_seq--;
                    }
                    end_idx = end_idx + EXPECTED_STRIDE;
                }

            } else{
                fprintf(opt.f_out, "%s\t", aln->data); //1
                fprintf(opt.f_out, "%" PRIu64 "\t", ns); //2
                uint32_t move_count = 0;
                uint32_t i = 1;
                uint32_t start_idx;
                uint32_t kmer_idx = 0;

                while(move_count < opt.move_start_offset){
                    int8_t value = bam_auxB2i(mv_array, i);
                    if(value == 1){
                        move_count++;
                        start_idx = i;
                    }
                    i++;
                }
                fprintf(opt.f_out, "%" PRIu64 "\t", ts + (i-2) * EXPECTED_STRIDE); //3

                uint32_t j = 1;
                uint64_t l_end_raw = 0;
                uint32_t len_seq_1 = len_seq+opt.move_start_offset;
                uint32_t end_idx = j + 1;
                for(; j<len_mv; j++){
                    int8_t value = bam_auxB2i(mv_array, j);
                    LOG_DEBUG("%d", value);
                    if(len_seq_1 > 0 && value == 1){
                        len_seq_1--;
                        end_idx = j;
                    }
                }
                if(len_seq_1 > 0 && j == len_mv){
                    l_end_raw =  (j-1)*EXPECTED_STRIDE+EXPECTED_STRIDE+(ns-j*EXPECTED_STRIDE);
                } else{
                    l_end_raw =  (end_idx-1)*EXPECTED_STRIDE;
                }

                fprintf(opt.f_out, "%" PRIu64 "\t", l_end_raw); //4
                fprintf(opt.f_out, "%s\t", "+"); //5
                fprintf(opt.f_out, "%s\t", aln->data); //6
                fprintf(opt.f_out, "%" PRIu32 "\t", len_seq + opt.kmer_size - 1); //7
                fprintf(opt.f_out, "%" PRIu32 "\t", kmer_idx); //8
                fprintf(opt.f_out, "%" PRIu32 "\t", len_seq - 1); //9
                fprintf(opt.f_out, "%" PRIu32 "\t", len_seq - kmer_idx); //10
                fprintf(opt.f_out, "%" PRIu32 "\t", len_seq + opt.kmer_size - 1); //11
                fprintf(opt.f_out, "%s\t", "255"); //12
                fprintf(opt.f_out, "%s", "ss:Z:"); //ss

                for(; i<len_mv; i++){
                    int8_t value = bam_auxB2i(mv_array, i);
                    LOG_DEBUG("%d", value);
                    if(len_seq > 0 && value == 1){
                        fprintf(opt.f_out, "%" PRIu32 ",", (i-start_idx)*EXPECTED_STRIDE); //ss
                        start_idx = i;
                        len_seq--;
                    } else if(len_seq > 0 && i == len_mv-1){
                        uint32_t l_duration =  (i-start_idx)*EXPECTED_STRIDE+EXPECTED_STRIDE+(ns-i*EXPECTED_STRIDE);
                        fprintf(opt.f_out, "%" PRIu32 ",", l_duration); //ss
                    }
                }
                fprintf(opt.f_out, "%s", "\n"); //newline

            }
        } else{
            ERROR("mv tag specification is incorrect%s", "");
            return -1;
        }
    }

    bam_destroy1(aln);
    bam_hdr_destroy(bam_hdr);
    sam_close(bam_fp);

    if (opt.arg_fname_out != NULL) {
        LOG_DEBUG("closing output file%s","");
        fclose(opt.f_out);
    }

    return 0;
}
