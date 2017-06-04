#!/bin/bash -lx

RUN_NAME="DEM07"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

if [ "$MULTIGRID" -eq "1" ]; then
    if [ "$TILED" -eq "1" ]; then
        exit -1 # unsupported
    else
        INPUTS=inputs_multiple
    fi
else
    if [ "$TILED" -eq "1" ]; then
        INPUTS=inputs_tiled
    else
        INPUTS=inputs_single
    fi
fi

rm -rf POST_* ${RUN_NAME}* &> /dev/null
time -p ${MFIX} "${INPUTS}"

post_dats=POST*.dat
for post_dats in ${post_dats}; do
    numdiff AUTOTEST/${post_dats} ${post_dats}
done
