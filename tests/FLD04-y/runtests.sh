#!/bin/bash -exl

set -euo pipefail

# set case directory
RUN_NAME="FLD04"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

if [ -n "$2" ]; then
    FEXTRACT=$2/fextract
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

rm -rf ${RUN_NAME}* POST_* &> /dev/null
time -p ${MPIRUN} "${MFIX}" "${INPUTS}"

if ! [ -z "${FEXTRACT}" ]; then
    ${FEXTRACT} -p FLD0400000/ -d 3 -v v_g -f 8 -s POST_UG.dat
    ${FEXTRACT} -p FLD0400000/ -d 2 -v w_g -f 8 -s POST_VG.dat

    post_dats=POST*.dat
    for result in ${post_dats}; do
        diff "AUTOTEST/${result}" "${result}"
    done
fi
