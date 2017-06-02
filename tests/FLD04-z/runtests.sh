#!/bin/bash -exl

set -euo pipefail

# set case directory
RUN_NAME="FLD01"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

if [ -n "$2" ]; then
    FEXTRACT=$2/fextract
fi

if [ "$ENABLE_MPI" -eq "1" ]; then
    MPIRUN="mpirun -np 4"
    POST_U=POST_UG_MPI.dat
    POST_V=POST_VG_MPI.dat
else
    MPIRUN=""
    POST_U=POST_UG.dat
    POST_V=POST_VG.dat
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

rm -rf ${RUN_NAME}* POST_* &> /dev/null
time -p ${MPIRUN} "${MFIX}" "${INPUTS}"

if ! [ -z "${FEXTRACT}" ]; then
  ${FEXTRACT} -p FLD0100000/ -d 1 -v w_g -s "${POST_U}"
  ${FEXTRACT} -p FLD0100000/ -d 3 -v u_g -s "${POST_V}"

  post_dats=POST*.dat
  for result in ${post_dats}; do
    numdiff -a 0.0 "AUTOTEST/${result}" "${result}"
  done
fi
