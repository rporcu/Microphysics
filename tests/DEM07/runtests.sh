#!/bin/bash -lx

set -euo pipefail

RUN_NAME="DEM07"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

GRID=${GRID:-"single multiple tiled"}

for grid_type in $GRID; do
    INPUTS=inputs_${grid_type}
    rm -rf POST_* ${RUN_NAME}* &> /dev/null
    time -p ${MFIX} "${INPUTS}"

    post_dats=POST*.dat
    for post_dats in ${post_dats}; do
        numdiff "AUTOTEST/${post_dats}" "${post_dats}"
    done

done
