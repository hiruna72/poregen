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
OUTPUT_DIR="${REL_PATH}/data/out/reform"
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
        valgrind --leak-check=full --error-exitcode=1 ./poregen reform "$@"
    else
        ./poregen reform "$@"
    fi
}

RAW_DIR="${REL_PATH}/data/raw/reform"
EXP_DIR="${REL_PATH}/data/exp/reform"

TESTCASE=1
info "testcase:$TESTCASE - help"
ex && die "testcase:$TESTCASE failed"

TESTCASE=2
info "testcase:$TESTCASE - read:1,kmer:0,move:1,output:paf"
ex -k0 -m1 -c "${RAW_DIR}/guppy_one_read.bam" > ${OUTPUT_DIR}/out.paf && die "testcase:$TESTCASE failed"

TESTCASE=3
info "testcase:$TESTCASE - read:1,kmer:9,move:-1,output:tsv"
ex -k0 -m0 -c "${RAW_DIR}/guppy_one_read.bam" > ${OUTPUT_DIR}/out.paf && die "testcase:$TESTCASE failed"

TESTCASE=4
info "testcase:$TESTCASE - read:1,kmer:9,move:10,output:paf"
ex -k9 -m10 -c "${RAW_DIR}/guppy_one_read.bam" > ${OUTPUT_DIR}/out.paf && die "testcase:$TESTCASE failed"

TESTCASE=5
info "testcase:$TESTCASE - read:1,kmer:1,move:0,output:paf"
ex -k1 -m0 -c "${RAW_DIR}/guppy_one_read.bam" > "${OUTPUT_DIR}/r1k1m0.paf" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/r1k1m0.paf" "${OUTPUT_DIR}/r1k1m0.paf" || die "testcase:${TESTCASE} diff failed"

TESTCASE=6
info "testcase:$TESTCASE - read:1,kmer:9,move:0,output:paf"
ex -k9 -m0 -c "${RAW_DIR}/guppy_one_read.bam" > "${OUTPUT_DIR}/r1k1m0.paf" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/r1k9m0.paf" "${OUTPUT_DIR}/r1k1m0.paf" || die "testcase:${TESTCASE} diff failed"

TESTCASE=7
info "testcase:$TESTCASE - read:1,kmer:9,move:1,output:paf"
ex -k9 -m1 -c "${RAW_DIR}/guppy_one_read.bam" > "${OUTPUT_DIR}/r1k9m1.paf" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/r1k9m1.paf" "${OUTPUT_DIR}/r1k9m1.paf" || die "testcase:${TESTCASE} diff failed"

TESTCASE=8
info "testcase:$TESTCASE - read:1,kmer:9,move:6,output:paf"
ex -k9 -m6 -c "${RAW_DIR}/guppy_one_read.bam" > "${OUTPUT_DIR}/r1k9m6.paf" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/r1k9m6.paf" "${OUTPUT_DIR}/r1k9m6.paf" || die "testcase:${TESTCASE} diff failed"

TESTCASE=9
info "testcase:$TESTCASE - read:1,kmer:9,move:8,output:paf"
ex -k9 -m8 -c "${RAW_DIR}/guppy_one_read.bam" > "${OUTPUT_DIR}/r1k9m8.paf" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/r1k9m8.paf" "${OUTPUT_DIR}/r1k9m8.paf" || die "testcase:${TESTCASE} diff failed"

TESTCASE=10
info "testcase:$TESTCASE - read:1,kmer:1,move:0,output:tsv"
ex -k1 -m0 "${RAW_DIR}/guppy_one_read.bam" > "${OUTPUT_DIR}/r1k1m0.tsv" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/r1k1m0.tsv" "${OUTPUT_DIR}/r1k1m0.tsv" || die "testcase:${TESTCASE} diff failed"

TESTCASE=11
info "testcase:$TESTCASE - read:1,kmer:9,move:0,output:tsv"
ex -k9 -m0 "${RAW_DIR}/guppy_one_read.bam" > "${OUTPUT_DIR}/r1k9m0.tsv" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/r1k9m0.tsv" "${OUTPUT_DIR}/r1k9m0.tsv" || die "testcase:${TESTCASE} diff failed"

TESTCASE=12
info "testcase:$TESTCASE - read:1,kmer:9,move:1,output:tsv"
ex -k9 -m1 "${RAW_DIR}/guppy_one_read.bam" > "${OUTPUT_DIR}/r1k9m1.tsv" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/r1k9m1.tsv" "${OUTPUT_DIR}/r1k9m1.tsv" || die "testcase:${TESTCASE} diff failed"

TESTCASE=13
info "testcase:$TESTCASE - read:1,kmer:9,move:6,output:tsv"
ex -k9 -m6 "${RAW_DIR}/guppy_one_read.bam" > "${OUTPUT_DIR}/r1k9m6.tsv" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/r1k9m6.tsv" "${OUTPUT_DIR}/r1k9m6.tsv" || die "testcase:${TESTCASE} diff failed"

TESTCASE=14
info "testcase:$TESTCASE - read:1,kmer:9,move:8,output:tsv"
ex -k9 -m8 "${RAW_DIR}/guppy_one_read.bam" > "${OUTPUT_DIR}/r1k9m8.tsv" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/r1k9m8.tsv" "${OUTPUT_DIR}/r1k9m8.tsv" || die "testcase:${TESTCASE} diff failed"

# dorado output
TESTCASE=15
info "testcase:$TESTCASE - read:2,kmer:9,move:8,output:tsv"
ex -k9 -m8 "${RAW_DIR}/slow5-dorado.sam" > "${OUTPUT_DIR}/dr2k9m8.tsv" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/dr2k9m8.tsv" "${OUTPUT_DIR}/dr2k9m8.tsv" || die "testcase:${TESTCASE} diff failed"

TESTCASE=16
info "testcase:$TESTCASE - read:2,kmer:9,move:8,output:paf"
ex -k9 -m8 -c "${RAW_DIR}/slow5-dorado.sam" > "${OUTPUT_DIR}/dr2k9m8.paf" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/dr2k9m8.paf" "${OUTPUT_DIR}/dr2k9m8.paf" || die "testcase:${TESTCASE} diff failed"

TESTCASE=17
info "testcase:$TESTCASE - read:2,kmer:1,move:0,output:tsv"
ex -k1 -m0 "${RAW_DIR}/slow5-dorado.sam" > "${OUTPUT_DIR}/dr2k1m0.tsv" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/dr2k1m0.tsv" "${OUTPUT_DIR}/dr2k1m0.tsv" || die "testcase:${TESTCASE} diff failed"

TESTCASE=18
info "testcase:$TESTCASE - read:2,kmer:1,move:0,output:paf"
ex -k1 -m0 -c "${RAW_DIR}/slow5-dorado.sam" > "${OUTPUT_DIR}/dr2k1m0.paf" || die "testcase:$TESTCASE failed"
diff "${EXP_DIR}/dr2k1m0.paf" "${OUTPUT_DIR}/dr2k1m0.paf" || die "testcase:${TESTCASE} diff failed"

info "all $TESTCASE testcases passed"
rm -r "$OUTPUT_DIR" || die "could not delete $OUTPUT_DIR"
exit 0
# If you want to log to the same file: command1 >> log_file 2>&1
# If you want different files: command1 >> log_file 2>> err_file
# use ANSI syntax format to view stdout/stderr on SublimeText
# use bash -n [script] and shellcheck [script] to check syntax
