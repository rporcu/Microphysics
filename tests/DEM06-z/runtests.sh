#!/bin/bash -lx

RUN_NAME="DEM06"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

if [ -n "$2" ]; then
    FCOMPARE=$2/plt_compare_diff_grids
fi

GRID=${GRID:-"single multiple tiled"}
MFIX_BENCHMARKS_HOME=${MFIX_BENCHMARKS_HOME:-}
FCOMPARE=${FCOMPARE:-}

for grid_type in $GRID; do
    INPUTS=inputs_${grid_type}
    rm -rf POST_* ${RUN_NAME}* &> /dev/null
    time -p ${MFIX} "${INPUTS}" DES_ONEWAY_COUPLED=.F.

    if ! [ -z "${MFIX_BENCHMARKS_HOME}" ] && ! [ -z "${FCOMPARE}" ]; then
        ${FCOMPARE} --infile1 "${MFIX_BENCHMARKS_HOME}/DEM06-z_plt00350" --infile2 DEM06_plt00350/
    fi

    post_dats=POST*.dat
    for result in ${post_dats}; do
        numdiff -a 0.0 "AUTOTEST/${result}" "${result}"
    done
done
