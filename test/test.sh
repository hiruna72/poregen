#!/bin/bash

RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}${1}${NC}" >&2 ; echo ; exit 1 ; } # terminate script
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

FIRST_FAILED_SET_OF_TESTCASES="NOT SET"
FLAG_FIRST_FAILED_SET_OF_TESTCASES=0
ret=0

if [ "$1" = 'mem' ]; then
    mem=1
else
    mem=0
fi

ex() {
    if [ $mem -eq 1 ]; then
        valgrind --leak-check=full --error-exitcode=1 "$@"
    else
        "$@"
    fi
}

fail() {
    echo 'FAILURE'
    ret=1
    if [ $FLAG_FIRST_FAILED_SET_OF_TESTCASES -eq 0 ]; then
        FLAG_FIRST_FAILED_SET_OF_TESTCASES=1
        FIRST_FAILED_SET_OF_TESTCASES=$1
    fi
}


#echo "DNA sref"
#ex  ./xyztool subtool1 test/example.blow5 > test/tmp.txt  || die "Running the tool failed"
#diff -q test/example.exp test/tmp.txt || die "diff failed"

TESTCASE_NAME="main_help"
info ${TESTCASE_NAME}
ex ./poregen  && die "Running the tool failed"


TESTCASE_NAME="sib_formater"
info $TESTCASE_NAME
if [ $mem -eq 1 ]; then
    if ! ./test/test_sigb_formater.sh mem; then
        fail "$TESTCASE_NAME"
    fi
else
    if ! ./test/test_sigb_formater.sh; then
        fail "$TESTCASE_NAME"
    fi
fi
if [ $ret -eq 1 ]; then
  info ">>>>>One or more test cases have failed. The first failed set of testcases is $FIRST_FAILED_SET_OF_TESTCASES<<<<<"
else
  info ">>>>>All tests passed<<<<<"
fi

exit $ret
