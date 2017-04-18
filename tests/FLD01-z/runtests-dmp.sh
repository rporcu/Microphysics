#!/bin/bash -exl

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

rm -rf POST_* ${RUN_NAME}* &> /dev/null
time -p  mpirun -np 4 ${MFIX} inputs

if ! [ -z "${MFIX_BENCHMARKS_HOME}" ] && ! [ -z "${FCOMPARE}" ]; then
  ${FCOMPARE} --infile1 ${MFIX_BENCHMARKS_HOME}/FLD01-z_FLD01-z_plt00000 --infile2 FLD0100000/
fi

if ! [ -z "${FEXTRACT}" ]; then
  ${FEXTRACT} -p FLD0100000/ -d 1 -v w_g && mv FLD0100000.slice POST_WG_MPI.dat
  ${FEXTRACT} -p FLD0100000/ -d 3 -v p_g && mv FLD0100000.slice POST_PG_MPI.dat

  post_dats=POST*.dat
  for result in ${post_dats}; do
    numdiff -a 0.0 -r 0.0 AUTOTEST/${result} ${result}
  done
fi
