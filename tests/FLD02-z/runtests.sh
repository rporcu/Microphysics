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

if [ "$ENABLE_MPI" -eq "1" ]; then
    MPIRUN="mpirun -np 2"
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

    ${FEXTRACT} -p FLD0200000/ -d 1 -v w_g -s POST_WG.dat
    ${FEXTRACT} -p FLD0200000/ -d 3 -v p_g -s POST_PG.dat

    post_dats=POST*.dat
    for result in ${post_dats}; do
        diff ${REL_ERR} "AUTOTEST/${result}" "${result}"
    done

    if ! [ -z "${MFIX_BENCHMARKS_HOME}" ] && ! [ -z "${FCOMPARE}" ]; then
        ${FCOMPARE} --infile1 "${MFIX_BENCHMARKS_HOME}/FLD02-z_FLD02-z_plt00000" --infile2 FLD0200000/
    fi
done
