#!/bin/bash -lx

set -euo pipefail

RUN_NAME="DEM07"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

INPUTS=inputs_single
if [ -n "$3" ]; then
    INPUTS=$3
fi
echo "Using INPUTS file ${INPUTS}"

rm -rf POST_* ${RUN_NAME}* &> /dev/null
time -p ${MFIX} ${INPUTS}

if [ ${INPUTS} = 'inputs_single' ]; then
   post_dats=POST*.dat
   for post_dats in ${post_dats}; do
       numdiff "AUTOTEST/${post_dats}" "${post_dats}"
   done
fi
