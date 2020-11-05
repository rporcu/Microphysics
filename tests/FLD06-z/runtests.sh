#!/bin/bash -exl

set -euo pipefail

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

if [ -n "$2" ]; then
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

if [ "$MFIX_MPI" -eq "1" ]; then
    if [ "$MFIX_OMP" -eq "1" ]; then
  MPIRUN="mpirun -np 2"
    else
  MPIRUN="mpirun -np 4"
    fi
else
    MPIRUN=""
fi

rm -rf POST_* const_plt* &> /dev/null
time -p ${MPIRUN} "${MFIX}" "${INPUTS}"

${FEXTRACT} -p FLD0600100/ -d 1 -v w_g -f 8 -s POST_VG.dat
${FEXTRACT} -p FLD0600100/ -d 3 -v p_g -f 8 -s POST_PG.dat

post_dats=POST*.dat
for result in ${post_dats}; do
    diff -b -u -I '#.*' "../FLD06-y/AUTOTEST/${result}" "${result}"
done
