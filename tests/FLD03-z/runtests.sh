#!/bin/bash -exl

set -euo pipefail

# set case directory
RUN_NAME="FLD03"

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
    if [ "$ENABLE_OMP" -eq "1" ]; then
	MPIRUN="mpirun -np 2"
    else 
	MPIRUN="mpirun -np 4"
    fi
else
    MPIRUN=""
fi

MFIX_BENCHMARKS_HOME=${MFIX_BENCHMARKS_HOME:-}
FCOMPARE=${FCOMPARE:-}

rm -rf POST_* ${RUN_NAME}* &> /dev/null
time -p ${MPIRUN} "${MFIX}" "${INPUTS}"

${FEXTRACT} -p FLD0300001/ -d 1 -v w_g -s POST_VG.dat
${FEXTRACT} -p FLD0300001/ -d 3 -v p_g -s POST_PG.dat

post_dats=POST*.dat
for result in ${post_dats}; do
    diff -u -I '#.*' "../FLD03-y/AUTOTEST/${result}" "${result}"
done

if ! [ -z "${MFIX_BENCHMARKS_HOME}" ] && ! [ -z "${FCOMPARE}" ]; then
    ${FCOMPARE} --infile1 "${MFIX_BENCHMARKS_HOME}/FLD03-z_FLD03-z_plt00001" --infile2 FLD0300001/
fi
