#!/bin/bash -lex

RUN_NAME="FLD02"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

rm -f POST_* &> /dev/null

for DELP_Z in -3.0 -2.0 -1.0 0.0 1.0 2.0 3.0; do
  rm -f ${RUN_NAME}* &> /dev/null
  time -p ${MFIX} inputs DELP_Z=${DELP_Z}
done

post_dats=../FLD02-x/AUTOTEST/POST*.dat

for test_post_file in ${post_dats}; do
    numdiff -a 0.000001 -r 0.05 ${test_post_file} \
      $(basename ${test_post_file})
done
