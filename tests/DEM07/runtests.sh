#!/bin/bash -lx

set -euo pipefail

RUN_NAME="DEM07"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

if [ -n "$2" ]; then
    FCOMPARE=$2/plt_compare_diff_grids
    FJOIN_PAR=$2/fjoin_par
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

rm -rf post_* ${RUN_NAME}* &> /dev/null
time -p ${MPIRUN} ${MFIX} ${INPUTS}

${FJOIN_PAR} -f DEM07_par --end 25 --var 100 --format 8 &> POST_GRAN_TEMP.NEW

post_dats=POST*.NEW
for result in ${post_dats}; do
    diff -w -B "AUTOTEST/${result}" "${result}"
done
