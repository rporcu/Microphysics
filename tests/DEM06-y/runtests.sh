#!/bin/bash -lx

set -euo pipefail

RUN_NAME="DEM06"

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
    if [ "$ENABLE_OMP" -eq "1" ]; then
	MPIRUN="mpirun -np 2"
    else
	MPIRUN="mpirun -np 4"
    fi
else
    MPIRUN=""
fi

FCOMPARE=${FCOMPARE:-}

rm -rf POST_* ${RUN_NAME}* &> /dev/null
time -p ${MPIRUN} ${MFIX} "${INPUTS}" TSTOP=0.150
time -p ${MPIRUN} ${MFIX} "${INPUTS}" "amr.restart=DEM06_chk00150"

${FJOIN_PAR} -f DEM06_par --end 350 --var  2 --format 4 --dt 0.001 &> POST_POS.NEW
${FJOIN_PAR} -f DEM06_par --end 350 --var 10 --format 4 --dt 0.001 &> POST_VEL.NEW

post_dats=POST*.NEW
for result in ${post_dats}; do
    diff -w -B "AUTOTEST/${result}" "${result}"
done
