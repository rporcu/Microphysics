#!/bin/bash -lex

RUN_NAME="DEM04"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

rm -f POST_* &> /dev/null

DES_IM=ADAMS_BASHFORTH
for DES_MEW in 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0; do
  rm -f ${RUN_NAME}* &> /dev/null
  time -p ${MFIX} inputs DES_INTG_METHOD=\"${DES_IM}\" \
    MEW=${DES_MEW} MEW_W=${DES_MEW}
done

post_dats=../DEM04-y/AUTOTEST/POST*.dat

for test_post_file in ${post_dats}; do
    numdiff -a 0.000001 -r 0.05  -X 1:4 -X 2:4 ${test_post_file} $(basename ${test_post_file})
done
