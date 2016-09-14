#!/bin/bash -l

RUN_NAME="DEM02"

MFIX=./mfix
if [ -n "$1" ]; then
   MFIX=$1
fi

rm -f ${RUN_NAME}* &> /dev/null

time -p ${MFIX}

#post_dats=AUTOTEST/POST*.dat
#
#for test_post_file in ${post_dats}; do
#    numdiff -a 0.000001 -r 0.05 ${test_post_file} $(basename ${test_post_file})
#done
