#!/bin/bash -lex

set -euo pipefail

RUN_NAME="FLD02"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

if [ -n "$2" ]; then
    FCOMPARE=$2/plt_compare_diff_grids
    FEXTRACT=$2/fextract
fi
if [ -z "${FEXTRACT}" ]; then
    echo "FEXTRACT is not set. Aborting test."
    exit 1
fi

INPUTS=inputs_single
if [ -n "$3" ]; then
    INPUTS=$3
fi
echo "Using INPUTS file ${INPUTS}"

if [ "$ENABLE_MPI" -eq "1" ]; then
    MPIRUN="mpirun -np 2"
else
    MPIRUN=""
fi

FCOMPARE=${FCOMPARE:-}

rm -rf POST_* ${RUN_NAME}* &> /dev/null
time -p ${MPIRUN} "${MFIX}" "${INPUTS}"

${FEXTRACT} -p FLD0200001/ -d 1 -t 1.0e-10 -v w_g -s POST_VG.dat
${FEXTRACT} -p FLD0200001/ -d 3 -t 1.0e-10 -v p_g -s POST_PG.dat

#post_dats=POST*.dat
#for result in ${post_dats}; do
#    diff -u -I '#.*' "../FLD02-y/AUTOTEST/${result}" "${result}"
#done
