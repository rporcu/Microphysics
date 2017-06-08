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

if [ "$ENABLE_MPI" -eq "1" ]; then
    MPIRUN="mpirun -np 4"
else
    MPIRUN=""
fi

GRID=${GRID:-"single multiple tiled"}

for grid_type in $GRID; do
    INPUTS=inputs_${grid_type}
    rm -rf ${RUN_NAME}* POST_* &> /dev/null
    time -p ${MPIRUN} "${MFIX}" "${INPUTS}"

    if ! [ -z "${FEXTRACT}" ]; then
    ${FEXTRACT} -p FLD0100000/ -d 1 -v w_g -s POST_UG.dat
    ${FEXTRACT} -p FLD0100000/ -d 3 -v u_g -s POST_VG.dat

    post_dats=POST*.dat
    for result in ${post_dats}; do
        numdiff -a 0.0 "AUTOTEST/${result}" "${result}"
    done
    fi
done
