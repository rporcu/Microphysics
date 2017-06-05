#!/bin/bash -exl

set -euo pipefail

# set case directory
RUN_NAME="FLD01"

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

if [ "$ENABLE_MPI" -eq "1" ]; then
    MPIRUN="mpirun -np 4"
    REL_ERR="-r 0.0"
else
    MPIRUN=""
    REL_ERR=""
fi

GRID=${GRID:-"single multiple tiled"}
MFIX_BENCHMARKS_HOME=${MFIX_BENCHMARKS_HOME:-}
FCOMPARE=${FCOMPARE:-}

for grid_type in $GRID; do
    INPUTS=inputs_${grid_type}
    rm -rf POST_* ${RUN_NAME}* &> /dev/null
    time -p ${MPIRUN} "${MFIX}" "${INPUTS}"

    ${FEXTRACT} -p FLD0100000/ -d 2 -v u_g -s POST_UG.dat
    ${FEXTRACT} -p FLD0100000/ -d 1 -v p_g -s POST_PG.dat

    post_dats=POST*.dat
    for result in ${post_dats}; do
        numdiff -a 0.0 ${REL_ERR} "AUTOTEST/${result}" "${result}"
    done

    if ! [ -z "${MFIX_BENCHMARKS_HOME}" ] && ! [ -z "${FCOMPARE}" ]; then
        ${FCOMPARE} --infile1 "${MFIX_BENCHMARKS_HOME}/FLD01-x_FLD01-x_plt00000" --infile2 FLD0100000/
    fi

done
