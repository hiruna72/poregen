#!/bin/bash

# steps
# cat slow5 files
# cat blow5 files
# convert catenated blow5 to slow5
# compare slow5s against the expected
# additional error catching testcases

RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo ; echo -e "${GREEN}$1${NC}" >&2 ; }
# Relative path to "slow5/tests/"
REL_PATH="$(dirname $0)/"

#...directories files tools arguments commands clean
OUTPUT_DIR="$REL_PATH/data/out/sib_formater"
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
        valgrind --leak-check=full --error-exitcode=1 ./poregen sigb_formater
    else
        ./poregen sigb_formater
    fi
}

#redirect
verbose=0
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

TESTCASE=1
info "testcase:$TESTCASE - help"
ex && die "testcase:$TESTCASE slow5tools cat failed"
exit 0

TESTCASE=2
info "testcase:$TESTCASE - cat two slow5s. output-stdout"
$SLOW5TOOLS cat "$RAW_DIR/slow5s/" > "$OUTPUT_DIR/output.slow5" || die "testcase:$TESTCASE slow5tools cat failed"
diff $EXP_SLOW5_FILE "$OUTPUT_DIR/output.slow5" || die "testcase:$TESTCASE diff failed"
slow5tools_quickcheck $OUTPUT_DIR

TESTCASE=3
info "testcase:$TESTCASE - cat two blow5s. output-stdout"
$SLOW5TOOLS cat "$RAW_DIR/blow5s/" > "$OUTPUT_DIR/output.blow5" || die "testcase:$TESTCASE slow5tools cat failed"
$SLOW5TOOLS view "$OUTPUT_DIR/output.blow5" > "$OUTPUT_DIR/output.slow5" || die "testcase:$TESTCASE slow5tools view failed"
diff $EXP_SLOW5_FILE "$OUTPUT_DIR/output.slow5" || die "testcase:$TESTCASE diff failed"
slow5tools_quickcheck $OUTPUT_DIR

TESTCASE=4
info "testcase:$TESTCASE - cat two slow5s. output-file"
$SLOW5TOOLS cat "$RAW_DIR/slow5s/" -o "$OUTPUT_DIR/output.slow5" || die "testcase:$TESTCASE slow5tools cat failed"
diff $EXP_SLOW5_FILE "$OUTPUT_DIR/output.slow5" || die "testcase:$TESTCASE diff failed"
slow5tools_quickcheck $OUTPUT_DIR

TESTCASE=5
info "testcase:$TESTCASE - cat two blow5s. output-file"
$SLOW5TOOLS cat "$RAW_DIR/blow5s/" -o "$OUTPUT_DIR/output.blow5" || die "testcase:$TESTCASE slow5tools cat failed"
$SLOW5TOOLS view "$OUTPUT_DIR/output.blow5" > "$OUTPUT_DIR/output.slow5" || die "testcase:$TESTCASE slow5tools view failed"
diff $EXP_SLOW5_FILE "$OUTPUT_DIR/output.slow5" || die "testcase:$TESTCASE diff failed"
slow5tools_quickcheck $OUTPUT_DIR

TESTCASE=6
info "testcase:$TESTCASE - cat two slow5s. output-file. wrong file extension. $SLOW5TOOLS_ERROR"
$SLOW5TOOLS cat "$RAW_DIR/slow5s/" -o "$OUTPUT_DIR/output.blow5" && die "testcase:$TESTCASE slow5tools cat failed"

TESTCASE=7
info "testcase:$TESTCASE - cat two blow5s. output-file. wrong file extension. $SLOW5TOOLS_ERROR"
$SLOW5TOOLS cat "$RAW_DIR/blow5s/" -o "$OUTPUT_DIR/output.slow5" && die "testcase:$TESTCASE slow5tools cat failed"

TESTCASE=8
info "testcase:$TESTCASE - cat different format files. $SLOW5TOOLS_ERROR"
$SLOW5TOOLS cat "$RAW_DIR/mixed_format/" > "$OUTPUT_DIR/output.slow5" && die "testcase:$TESTCASE slow5tools cat failed"

TESTCASE=9
info "testcase:$TESTCASE - cat different compression types files. $SLOW5TOOLS_ERROR"
$SLOW5TOOLS cat "$RAW_DIR/mixed_compression/" > "$OUTPUT_DIR/output.slow5" && die "testcase:$TESTCASE slow5tools cat failed"

TESTCASE=10
info "testcase:$TESTCASE - cat different auxiliary attribute order. $SLOW5TOOLS_ERROR"
$SLOW5TOOLS cat "$RAW_DIR/different_aux_order/" > "$OUTPUT_DIR/output.slow5" && die "testcase:$TESTCASE slow5tools cat failed"

info "all $TESTCASE cat testcases passed"
rm -r "$OUTPUT_DIR" || die "could not delete $OUTPUT_DIR"
exit 0
# If you want to log to the same file: command1 >> log_file 2>&1
# If you want different files: command1 >> log_file 2>> err_file
# use ANSI syntax format to view stdout/stderr on SublimeText
# use bash -n [script] and shellcheck [script] to check syntax
