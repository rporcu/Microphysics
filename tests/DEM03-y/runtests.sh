#!/bin/bash -lex

RUN_NAME="DEM03"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

rm -f POST_* &> /dev/null

rm -f ${RUN_NAME}* &> /dev/null
time -p ${MFIX} inputs

post_dats=AUTOTEST/POST*.dat

for test_post_file in ${post_dats}; do
    diff ${test_post_file} $(basename ${test_post_file})
done
