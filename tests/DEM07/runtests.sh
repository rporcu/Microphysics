#!/bin/bash -lx

RUN_NAME="DEM07"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

rm -rf POST_* ${RUN_NAME}* &> /dev/null
time -p ${MFIX} inputs

# post_dats=POST*.dat
# for post_dats in ${post_dats}; do
#     numdiff AUTOTEST/${post_dats} ${post_dats}
# done
