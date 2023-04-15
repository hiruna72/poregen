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
OUTPUT_DIR="${REL_PATH}/data/out/gmove"
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
        valgrind --leak-check=full --error-exitcode=1 ./poregen gmove "$@"
    else
        ./poregen gmove "$@"
    fi
}

RAW_DIR="${REL_PATH}/data/raw/gmove"
EXP_DIR="${REL_PATH}/data/exp/gmove"

TESTCASE=1
info "testcase:$TESTCASE - help"
ex && die "testcase:$TESTCASE failed"

TESTCASE=2
info "testcase:$TESTCASE - input:table"
ex "${RAW_DIR}/single_read/reads.slow5" "${RAW_DIR}/single_read/guppy_move" --file_limit 50 ${OUTPUT_DIR}/${TESTCASE} || die "testcase:$TESTCASE failed"

TESTCASE=3
info "testcase:$TESTCASE - input:table kmer_file:given"
ex "${RAW_DIR}/single_read/reads.slow5" "${RAW_DIR}/single_read/guppy_move" --kmer_file "${RAW_DIR}/single_read/kmer_file.txt" ${OUTPUT_DIR}/${TESTCASE} && die "testcase:$TESTCASE failed"

TESTCASE=4
info "testcase:$TESTCASE - input:table kmer_file:given"
ex -k 6 -m 0 "${RAW_DIR}/single_read/reads.slow5" "${RAW_DIR}/single_read/guppy_move" --kmer_file "${RAW_DIR}/single_read/kmer_file.txt" ${OUTPUT_DIR}/${TESTCASE} || die "testcase:$TESTCASE failed"

TESTCASE=5
info "testcase:$TESTCASE - input:paf"
ex "${RAW_DIR}/single_read/reads.slow5" "${RAW_DIR}/single_read/guppy_move.paf" --file_limit 50 ${OUTPUT_DIR}/${TESTCASE} && die "testcase:$TESTCASE failed"

TESTCASE=6
info "testcase:$TESTCASE - input:paf"
ex "${RAW_DIR}/single_read/reads.slow5" "${RAW_DIR}/single_read/guppy_move.paf" --file_limit 50 ${OUTPUT_DIR}/${TESTCASE} --fastq "${RAW_DIR}/single_read/read_0.fastq" || die "testcase:$TESTCASE failed"

TESTCASE=7
info "testcase:$TESTCASE - input:paf kmer_file:given"
ex -k 6 "${RAW_DIR}/single_read/reads.slow5" "${RAW_DIR}/single_read/guppy_move.paf" ${OUTPUT_DIR}/${TESTCASE} --fastq "${RAW_DIR}/single_read/read_0.fastq" --kmer_file "${RAW_DIR}/single_read/kmer_file.txt" || die "testcase:$TESTCASE failed"

diff ${OUTPUT_DIR}/4/freq.txt ${OUTPUT_DIR}/7/freq.txt
diff ${OUTPUT_DIR}/4/dump ${OUTPUT_DIR}/7/dump

info "all $TESTCASE testcases passed"
#rm -r "$OUTPUT_DIR" || die "could not delete $OUTPUT_DIR"
exit 0




























