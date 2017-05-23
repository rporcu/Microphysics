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

rm -rf ${RUN_NAME}* POST_* &> /dev/null
time -p ${MFIX} inputs

if ! [ -z "${FEXTRACT}" ]; then
  ${FEXTRACT} -p FLD0100000/ -d 1 -v v_g && mv FLD0100000.slice POST_UG.dat
  ${FEXTRACT} -p FLD0100000/ -d 2 -v u_g && mv FLD0100000.slice POST_VG.dat

  post_dats=POST*.dat
  for result in ${post_dats}; do
    numdiff -a 0.0 AUTOTEST/${result} ${result}
  done
fi
