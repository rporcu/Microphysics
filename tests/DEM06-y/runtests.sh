#!/bin/bash -lx

set -euo pipefail

RUN_NAME="DEM06"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

if [ -n "$2" ]; then
    FCOMPARE=$2/plt_compare_diff_grids
fi

INPUTS=inputs_single
if [ -n "$3" ]; then
    INPUTS=$3
fi
echo "Using INPUTS file ${INPUTS}"

if [ "$ENABLE_MPI" -eq "1" ]; then
    MPIRUN="mpirun -np 4"
else
    MPIRUN=""
fi

MFIX_BENCHMARKS_HOME=${MFIX_BENCHMARKS_HOME:-}
FCOMPARE=${FCOMPARE:-}

rm -rf POST_* ${RUN_NAME}* &> /dev/null
time -p ${MPIRUN} ${MFIX} "${INPUTS}" DES_ONEWAY_COUPLED=.F.

post_dats=POST*.dat
for result in ${post_dats}; do
    diff -w -B "AUTOTEST/${result}" "${result}"
done
