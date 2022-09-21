/* @file main.c
**
******************************************************************************/
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "error.h"
#include "misc.h"
#include "xyztool.h"

#ifdef HAVE_EXECINFO_H
    #include <execinfo.h>
#endif

extern enum xyztool_log_level_opt xyztool_log_level;

//make the segmentation faults a bit cool
void sig_handler(int sig) {
#ifdef HAVE_EXECINFO_H
    void* array[100];
    size_t size = backtrace(array, 100);
    ERROR("An unfortunate event of a segmentation fault has occurred. But at least "
          "it is better than a wrong answer%s",
          ".");
    fprintf(stderr,
            "[%s::DEBUG]\033[1;35m Here is the backtrace in case it is of any "
            "use:\n",
            __func__);
    backtrace_symbols_fd(&array[2], size - 1, STDERR_FILENO);
    fprintf(stderr, "\033[0m\n");
#else
    ERROR("An unfortunate event of a segmentation fault has occurred. But at least "
          "it is better than a wrong answer%s",
          ".");
#endif
    exit(EXIT_FAILURE);
}

int subtool1_main(int argc, char* argv[]);

int print_usage(FILE *fp_help){

    fprintf(fp_help,"Usage: xyztool <command> [options]\n\n");
    fprintf(fp_help,"command:\n");
    fprintf(fp_help,"         subtool1      do something\n");
    fprintf(fp_help,"         subtool2      do something\n");
    if(fp_help==stderr){
        return(EXIT_FAILURE);
    }
    else if(fp_help==stdout){
        return(EXIT_SUCCESS);
    } else {
        return(EXIT_FAILURE);
    }

}

int main(int argc, char* argv[]){

    double realtime0 = realtime();
    signal(SIGSEGV, sig_handler);

    int ret=1;

    if(argc<2){
        return print_usage(stderr);
    } else if (strcmp(argv[1],"subtool1")==0){
        ret=subtool1_main(argc-1, argv+1);
    } else if (strcmp(argv[1],"subtool2")==0){
        ret=subtool1_main(argc-1, argv+1);
    } else if(strcmp(argv[1],"--version")==0 || strcmp(argv[1],"-V")==0){
        fprintf(stdout,"xyztool %s\n",XYZTOOL_VERSION);
        exit(EXIT_SUCCESS);
    } else if(strcmp(argv[1],"--help")==0 || strcmp(argv[1],"-h")==0){
        return print_usage(stdout);
    } else{
        fprintf(stderr,"[xyztool] Unrecognised command %s\n",argv[1]);
        return print_usage(stderr);
    }

    fprintf(stderr,"[%s] Version: %s\n", __func__,XYZTOOL_VERSION);
    fprintf(stderr, "[%s] CMD:", __func__);
    for (int i = 0; i < argc; ++i) {
        fprintf(stderr, " %s", argv[i]);
    }

    fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU time: %.3f sec; Peak RAM: %.3f GB\n\n",
            __func__, realtime() - realtime0, cputime(),peakrss() / 1024.0 / 1024.0 / 1024.0);

    return ret;
}
