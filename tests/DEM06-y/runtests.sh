#!/bin/bash -lx

RUN_NAME="DEM06"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

if [ -n "$2" ]; then
    FCOMPARE=$2/plt_compare_diff_grids
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
time -p ${MFIX} "${INPUTS}" DES_ONEWAY_COUPLED=.F.

if ! [ -z "${MFIX_BENCHMARKS_HOME}" ] && ! [ -z "${FCOMPARE}" ]; then
  ${FCOMPARE} --infile1 "${MFIX_BENCHMARKS_HOME}/DEM06-y_plt00350" --infile2 DEM0600350/
fi

post_dats=POST*.dat
for result in ${post_dats}; do
    numdiff -a 0.0 "AUTOTEST/${result}" "${result}"
done
