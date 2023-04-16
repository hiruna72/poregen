#!/bin/bash

# steps

RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo ; echo -e "${GREEN}$1${NC}" >&2 ; }

#redirect
verbose=1
exec 3>&1
exec 4>&2
if ((verbose)); then
  echo "verbose=1"
else
  echo "verbose=0"
  exec 1>/dev/null
  exec 2>/dev/null
fi
#echo "this should be seen if verbose"
#echo "this should always be seen" 1>&3 2>&4

# Relative path to "slow5/tests/"
REL_PATH="$(dirname $0)/"
#...directories files tools arguments commands clean
OUTPUT_DIR="${REL_PATH}/data/out/kmer_freq"
test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"
#commands ...

if [ "$1" = 'mem' ]; then
    mem=1
else
    mem=0
fi

ex() {
    if [ $mem -eq 1 ]; then
        valgrind --leak-check=full --error-exitcode=1 ./poregen kmer_freq "$@"
    else
        ./poregen kmer_freq "$@"
    fi
}

RAW_DIR="${REL_PATH}/data/raw/kmer_freq"
EXP_DIR="${REL_PATH}/data/exp/kmer_freq"

TESTCASE=1
info "testcase:$TESTCASE - help"
ex && die "testcase:$TESTCASE failed"

TESTCASE=2
info "testcase:$TESTCASE - basic"
ex 6 "${RAW_DIR}/read_0.fastq" -o ${OUTPUT_DIR}/${TESTCASE}.txt || die "testcase:$TESTCASE failed"

TESTCASE=3
info "testcase:$TESTCASE - print_absent_kmers:0"
ex 6 "${RAW_DIR}/read_0.fastq" -o ${OUTPUT_DIR}/${TESTCASE}.txt --print_absent_kmers 0 || die "testcase:$TESTCASE failed"

TESTCASE=4
info "testcase:$TESTCASE - print_absent_kmers:0 sort:1"
ex 6 "${RAW_DIR}/read_0.fastq" -o ${OUTPUT_DIR}/${TESTCASE}.txt --print_absent_kmers 0 --sort 1 || die "testcase:$TESTCASE failed"

TESTCASE=5
info "testcase:$TESTCASE - print_absent_kmers:1 sort:1"
ex 6 "${RAW_DIR}/read_0.fastq" -o ${OUTPUT_DIR}/${TESTCASE}.txt --print_absent_kmers 1 --sort 1 || die "testcase:$TESTCASE failed"

TESTCASE=6
info "testcase:$TESTCASE - print_absent_kmers:0 sort:2"
ex 6 "${RAW_DIR}/read_0.fastq" -o ${OUTPUT_DIR}/${TESTCASE}.txt --print_absent_kmers 0 --sort 1 || die "testcase:$TESTCASE failed"

#info "all $TESTCASE testcases passed"
#rm -r "$OUTPUT_DIR" || die "could not delete $OUTPUT_DIR"
exit 0
# If you want to log to the same file: command1 >> log_file 2>&1
# If you want different files: command1 >> log_file 2>> err_file
# use ANSI syntax format to view stdout/stderr on SublimeText
# use bash -n [script] and shellcheck [script] to check syntax
































